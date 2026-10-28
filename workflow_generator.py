#!/usr/bin/env python3

"""Pegasus workflow generator for gw-workflow.

Independent reproduction of the GWTC-1 binary black hole detections from
publicly released GWOSC strain data. Built to SPEC.md; sourcing rules in
CLAUDE.md.

Pipeline (SPEC.md §7.1):
  fetch_strain      GWOSC 4096 s / 4 kHz segment per event, per detector
  check_quality     DQ + injection bitmasks and NaN gate for the analysis span
  condition_strain  decimate, PSD, inverse-PSD truncation, gating, off-source spans
  bank_psd          harmonic-mean placement PSD over every conditioned segment
  bank_row          one constant-mass-ratio row of the bank, per q (fan-out)
  merge_bank        concatenate the rows, apply the B6 adequacy gate
  matched_filter    per (event, detector, bank chunk); also drives the null test
  coincidence       exact-template, 15 ms coincidence per event / per off-source span
  estimate_params   detector-frame chirp mass, merger time, network SNR
  spectrogram       whitened Q-transform around each event
  validate_event    compare against GWTC-1 with explicit pass/fail criteria
  control_offsource V5a null test (hard gate)
  control_injection V5b injection recovery (hard gate)
  aggregate_report  one report across every event and control

Usage:
    ./workflow_generator.py                          # Tier A, H1+L1
    ./workflow_generator.py --tier all
    ./workflow_generator.py --smoke-test
    ./workflow_generator.py --events GW150914 --detectors H1 L1
"""

import argparse
import json
import logging
import math
import os
import sys
from pathlib import Path

from Pegasus.api import *

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# SPEC.md §1.2 -- tiering is a property of the events, not of this pipeline.
TIER_A = ["GW150914", "GW151226", "GW170104", "GW170608",
          "GW170809", "GW170814", "GW170823"]
TIER_B = ["GW151012", "GW170729", "GW170818"]
# SPEC.md §1.2 / N11: GW170817 is a BNS and is out of scope.
EXCLUDED = {"GW170817"}

# SPEC.md §2 / W3: every job takes one core and <= 1.5 GB so four pack into an
# 8 GB node with headroom. The bank build and the report may take 2 GB.
TOOL_CONFIGS = {
    "fetch_strain":      {"memory": "1 GB",   "cores": 1},
    "check_quality":     {"memory": "1 GB",   "cores": 1},
    "condition_strain":  {"memory": "1500 MB", "cores": 1},
    "bank_psd":          {"memory": "1500 MB", "cores": 1},
    "bank_row":          {"memory": "1500 MB", "cores": 1},
    "merge_bank":        {"memory": "2 GB",   "cores": 1},
    "matched_filter":    {"memory": "1500 MB", "cores": 1},
    "coincidence":       {"memory": "1 GB",   "cores": 1},
    "estimate_params":   {"memory": "1 GB",   "cores": 1},
    "spectrogram":       {"memory": "1500 MB", "cores": 1},
    "validate_event":    {"memory": "1 GB",   "cores": 1},
    "control_offsource": {"memory": "1 GB",   "cores": 1},
    "control_injection": {"memory": "1500 MB", "cores": 1},
    "aggregate_report":  {"memory": "2 GB",   "cores": 1},
}

# Support files: staged as DATA into each job's working directory and imported
# with os.getcwd() (PEGASUS.md, "Finding Support Files"). gwlib.py is the shared
# signal-processing core; matched_filter.py is also a support file because
# control_injection.py imports filter_segment() from it, so the null test and the
# injection test cannot drift apart from the on-source filter.
SUPPORT_FILES = ["gwlib.py", "matched_filter.py"]
REFERENCE_JSON = "gwtc1_reference.json"


def split_list(values, cast=str, flag=""):
    """Accept both space- and comma-separated list arguments.

    The Pegasus AI Studio GUI passes a list argument as ONE comma-separated
    string (--detectors H1,L1), while a shell passes it space-separated
    (--detectors H1 L1). argparse's nargs="+" only understands the latter, and
    with type=float a comma-joined string fails outright:

        error: argument --injection-chirp-masses: invalid float value: '8,15,25,40,55'

    So list arguments are taken as strings and split here, accepting either form
    (and any mixture).
    """
    out = []
    for raw in values or []:
        for part in str(raw).replace(",", " ").split():
            try:
                out.append(cast(part))
            except (TypeError, ValueError):
                raise argparse.ArgumentTypeError(
                    f"{flag or 'list argument'}: {part!r} is not a valid "
                    f"{getattr(cast, '__name__', cast)}"
                )
    return out


class GwWorkflow:
    """GWTC-1 BBH reproduction from GWOSC open data."""

    wf = sc = tc = rc = props = None
    wf_name = "gw_workflow"

    def __init__(self, dagfile="workflow.yml"):
        self.dagfile = dagfile
        self.wf_dir = str(Path(__file__).parent.resolve())
        self.shared_scratch_dir = os.path.join(self.wf_dir, "scratch")
        self.local_storage_dir = os.path.join(self.wf_dir, "output")
        self.support = {}          # lfn -> File
        self.reference = None

    # ------------------------------------------------------------------
    def write(self):
        if self.sc is not None:
            self.sc.write()
        self.props.write()
        self.rc.write()
        self.tc.write()
        self.wf.write(file=self.dagfile)

    def create_pegasus_properties(self):
        self.props = Properties()
        self.props["pegasus.transfer.threads"] = "16"
        # CLAUDE.md §4 / SPEC.md W4: GWOSC is a shared public service. Cap the
        # concurrent downloads via a DAGMan category, not politeness alone.
        self.props["dagman.download.maxjobs"] = "4"

    def create_sites_catalog(self, exec_site_name="condorpool"):
        self.sc = SiteCatalog()
        local = Site("local").add_directories(
            Directory(Directory.SHARED_SCRATCH, self.shared_scratch_dir)
            .add_file_servers(
                FileServer("file://" + self.shared_scratch_dir, Operation.ALL)),
            Directory(Directory.LOCAL_STORAGE, self.local_storage_dir)
            .add_file_servers(
                FileServer("file://" + self.local_storage_dir, Operation.ALL)),
        )
        exec_site = (
            Site(exec_site_name)
            .add_condor_profile(universe="vanilla")
            .add_pegasus_profile(style="condor")
        )
        self.sc.add_sites(local, exec_site)

    # ------------------------------------------------------------------
    def create_transformation_catalog(self, exec_site_name="condorpool",
                                      container_image=None):
        self.tc = TransformationCatalog()

        # SPEC.md W7: one container for every job. image_site is passed by
        # KEYWORD -- the 4th positional argument is `arguments`, and passing a
        # site name there leaves image_site unset and fails at plan time with an
        # unresolvable CONTAINER_SITE placeholder (PEGASUS.md).
        sif = container_image or os.path.join(self.wf_dir, "GW_Container.sif")
        container = Container(
            "gw_container",
            container_type=Container.SINGULARITY,
            image="file://" + os.path.abspath(sif),
            image_site="local",
        )

        transformations = []
        for tool, cfg in TOOL_CONFIGS.items():
            tx = Transformation(
                tool,
                site=exec_site_name,
                pfn=os.path.join(self.wf_dir, f"bin/{tool}.py"),
                is_stageable=True,
                container=container,
            ).add_pegasus_profile(memory=cfg["memory"], cores=cfg["cores"])
            if tool == "fetch_strain":
                # SPEC.md W4: rate-limit GWOSC downloads.
                tx.add_profiles(Namespace.DAGMAN, key="category", value="download")
            transformations.append(tx)

        self.tc.add_containers(container)
        self.tc.add_transformations(*transformations)

    # ------------------------------------------------------------------
    def create_replica_catalog(self):
        """Support files and the reference values.

        No strain files are registered: the first job fetches them from GWOSC at
        runtime (SPEC.md D1), so there is nothing to stage in.
        """
        self.rc = ReplicaCatalog()
        for lfn in SUPPORT_FILES:
            path = os.path.join(self.wf_dir, "bin", lfn)
            if not os.path.exists(path):
                raise FileNotFoundError(f"support file missing: {path}")
            self.rc.add_replica("local", lfn, "file://" + path)
            self.support[lfn] = File(lfn)

        ref = os.path.join(self.wf_dir, "references", REFERENCE_JSON)
        if not os.path.exists(ref):
            raise FileNotFoundError(f"reference values missing: {ref}")
        self.rc.add_replica("local", REFERENCE_JSON, "file://" + ref)
        self.reference = File(REFERENCE_JSON)

    # ------------------------------------------------------------------
    def _job(self, tool, _id, label=None):
        """A job with the shared support files already attached as inputs."""
        j = Job(tool, _id=_id, node_label=_id)
        j.add_inputs(*self.support.values())
        if label:
            j.add_pegasus_profiles(label=label)
        return j

    def create_workflow(self, args):
        self.wf = Workflow(self.wf_name, infer_dependencies=True)

        conditioned = {}      # (event, det) -> File
        cond_reports = []
        spec_reports = []
        validations = []
        controls = []

        # ---- per event, per detector: fetch -> quality -> condition -----
        for event in args.events:
            for det in args.detectors:
                tag = f"{event}_{det}"
                strain = File(f"{tag}_strain.hdf5")
                meta = File(f"{tag}_fetch.json")
                quality = File(f"{tag}_quality.json")
                cond = File(f"{tag}_conditioned.hdf5")
                creport = File(f"{tag}_conditioning.json")

                fetch = (
                    self._job("fetch_strain", f"fetch_{tag}", label=event)
                    .add_args("--event", event, "--detector", det,
                              "--duration", str(args.segment_duration),
                              "--sample-rate", str(args.source_sample_rate),
                              "--output", strain, "--meta-output", meta)
                    .add_outputs(strain, stage_out=False, register_replica=False)
                    .add_outputs(meta, stage_out=False, register_replica=False)
                )

                qual = (
                    self._job("check_quality", f"quality_{tag}", label=event)
                    .add_args("--strain", strain, "--meta", meta,
                              "--analysis-length", str(args.analysis_length),
                              "--output", quality)
                    .add_inputs(strain, meta)
                    .add_outputs(quality, stage_out=True, register_replica=False)
                )

                condition = (
                    self._job("condition_strain", f"condition_{tag}", label=event)
                    .add_args("--strain", strain, "--meta", meta,
                              "--quality", quality,
                              "--sample-rate", str(args.sample_rate),
                              "--f-low", str(args.f_low),
                              "--f-high", str(args.f_high),
                              "--analysis-length", str(args.analysis_length),
                              "--psd-length", str(args.psd_length),
                              "--offsource-length", str(args.analysis_length),
                              "--n-offsource", str(args.n_offsource),
                              "--highpass", str(args.highpass),
                              "--gate-threshold", str(args.gate_threshold),
                              "--output", cond, "--report", creport)
                    .add_inputs(strain, meta, quality)
                    .add_outputs(cond, stage_out=False, register_replica=False)
                    .add_outputs(creport, stage_out=True, register_replica=False)
                )
                self.wf.add_jobs(fetch, qual, condition)
                conditioned[(event, det)] = cond
                cond_reports.append(creport)

                if not args.skip_spectrograms:
                    png = File(f"{tag}_spectrogram.png")
                    srep = File(f"{tag}_spectrogram.json")
                    self.wf.add_jobs(
                        self._job("spectrogram", f"spectrogram_{tag}", label=event)
                        .add_args("--conditioned", cond,
                                  "--f-low", str(args.f_low),
                                  "--f-high", str(args.f_high),
                                  "--output", png, "--report", srep)
                        .add_inputs(cond)
                        .add_outputs(png, stage_out=True, register_replica=False)
                        .add_outputs(srep, stage_out=True, register_replica=False)
                    )
                    spec_reports.append(srep)

        # ---- bank: one, shared, from the harmonic-mean PSD --------------
        # SPEC.md B7 + Usman §3.1. This is a genuine fan-in: exact-match
        # coincidence requires every detector filtered against the SAME bank,
        # and the placement PSD is the harmonic mean over these segments.
        #
        # Built in three stages rather than one job. Placement at fixed mass
        # ratio never consults another row, so the rows fan out; measured
        # serially the default bank (3453 templates) takes ~29 min on one core,
        # which would sit in front of every filter job. The PSD is split out
        # first so the rows stage one ~4 MB spectrum instead of every
        # conditioned segment (~30 MB x N) N_q times over.
        bank_psd = File("bank_psd.hdf5")
        bank_psd_report = File("bank_psd.json")
        psd_job = (
            self._job("bank_psd", "bank_psd")
            .add_args("--output", bank_psd, "--report", bank_psd_report)
            .add_outputs(bank_psd, stage_out=False, register_replica=False)
            .add_outputs(bank_psd_report, stage_out=True, register_replica=False)
        )
        for cond in conditioned.values():
            psd_job.add_args("--psd", cond)
            psd_job.add_inputs(cond)
        self.wf.add_jobs(psd_job)

        # The q grid is computed here and each value passed explicitly, so the
        # rows cannot disagree with one another about the grid they came from.
        # Log-spaced from 1 to q_max, matching numpy.geomspace, computed in
        # plain Python so the submit host needs no numpy.
        q_grid = ([1.0] if args.n_q < 2 else
                  [math.exp(math.log(args.q_max) * i / (args.n_q - 1))
                   for i in range(args.n_q)])
        row_files = []
        for i, q in enumerate(q_grid):
            row = File(f"bank_row_{i:02d}.json")
            self.wf.add_jobs(
                self._job("bank_row", f"bank_row_{i:02d}")
                .add_args("--psd", bank_psd, "--q", repr(float(q)),
                          "--row", str(i),
                          "--mc-min", str(args.mc_min), "--mc-max", str(args.mc_max),
                          "--minimal-match", str(args.minimal_match),
                          "--f-low", str(args.f_low), "--f-high", str(args.f_high),
                          "--approximant", args.approximant,
                          "--output", row)
                .add_inputs(bank_psd)
                .add_outputs(row, stage_out=False, register_replica=False)
            )
            row_files.append(row)

        bank = File("bank.hdf5")
        bank_manifest = File("bank_manifest.json")
        merge_job = (
            self._job("merge_bank", "merge_bank")
            .add_args("--psd", bank_psd, "--psd-report", bank_psd_report,
                      "--reference", self.reference,
                      "--tier-a", *[e for e in args.events if e in TIER_A],
                      "--output", bank, "--manifest", bank_manifest)
            .add_inputs(bank_psd, bank_psd_report, self.reference)
            .add_outputs(bank, stage_out=False, register_replica=False)
            .add_outputs(bank_manifest, stage_out=True, register_replica=False)
        )
        for row in row_files:
            merge_job.add_args("--row", row)
            merge_job.add_inputs(row)
        self.wf.add_jobs(merge_job)

        # ---- per event: filter -> coincidence -> estimate -> validate ---
        for event in args.events:
            tier = "A" if event in TIER_A else "B"
            trig_files = []
            for det in args.detectors:
                cond = conditioned[(event, det)]
                for k in range(args.bank_chunks):
                    trig = File(f"{event}_{det}_triggers_{k:03d}.json")
                    self.wf.add_jobs(
                        self._job("matched_filter", f"filter_{event}_{det}_{k:03d}",
                                  label=event)
                        .add_args("--conditioned", cond, "--bank", bank,
                                  "--dataset", "analysis",
                                  "--chunk", str(k),
                                  "--n-chunks", str(args.bank_chunks),
                                  "--snr-threshold", str(args.snr_threshold),
                                  "--output", trig)
                        .add_inputs(cond, bank)
                        .add_outputs(trig, stage_out=False, register_replica=False)
                    )
                    trig_files.append(trig)

            coinc = File(f"{event}_coincidence.json")
            cj = (
                self._job("coincidence", f"coincidence_{event}", label=event)
                .add_args("--event", event,
                          "--coincidence-window-ms", str(args.coincidence_window_ms),
                          "--output", coinc)
                .add_inputs(*trig_files)
                .add_outputs(coinc, stage_out=True, register_replica=False)
            )
            for t in trig_files:
                cj.add_args("--triggers", t)

            estimates = File(f"{event}_estimates.json")
            ej = (
                self._job("estimate_params", f"estimate_{event}", label=event)
                .add_args("--coincidence", coinc, "--reference", self.reference,
                          "--event", event, "--output", estimates)
                .add_inputs(coinc, self.reference)
                .add_outputs(estimates, stage_out=True, register_replica=False)
            )

            validation = File(f"{event}_validation.json")
            vj = (
                self._job("validate_event", f"validate_{event}", label=event)
                .add_args("--estimates", estimates, "--coincidence", coinc,
                          "--bank-manifest", bank_manifest,
                          "--reference", self.reference,
                          "--event", event, "--tier", tier,
                          "--output", validation)
                .add_inputs(estimates, coinc, bank_manifest, self.reference)
                .add_outputs(validation, stage_out=True, register_replica=False)
            )
            self.wf.add_jobs(cj, ej, vj)
            validations.append(validation)
            if event == args.control_event:
                self._control_onsource = coinc

        # ---- controls (SPEC.md V5) --------------------------------------
        if not args.skip_controls:
            controls += self._add_offsource_control(args, conditioned, bank)
            controls += self._add_injection_control(args, conditioned, bank)

        # ---- fan-in report ---------------------------------------------
        report_md = File("report.md")
        report_json = File("report.json")
        rj = (
            self._job("aggregate_report", "aggregate_report")
            .add_args("--bank-manifest", bank_manifest,
                      "--output-md", report_md, "--output-json", report_json)
            .add_inputs(bank_manifest, *validations, *controls,
                        *cond_reports, *spec_reports)
            .add_outputs(report_md, stage_out=True, register_replica=False)
            .add_outputs(report_json, stage_out=True, register_replica=False)
        )
        for v in validations:
            rj.add_args("--validation", v)
        for c in controls:
            rj.add_args("--control", c)
        for c in cond_reports:
            rj.add_args("--conditioning", c)
        for s in spec_reports:
            rj.add_args("--spectrogram", s)
        for e in args.events:
            # Evidence census: the report must see a usable validation for
            # every event this DAG was generated for, or verdict INCOMPLETE.
            rj.add_args("--expected-event", e)
        if args.skip_controls:
            # SPEC.md §7.3: skipping controls forces INCOMPLETE, never PASS.
            rj.add_args("--controls-skipped")
        self.wf.add_jobs(rj)

    # ------------------------------------------------------------------
    def _add_offsource_control(self, args, conditioned, bank):
        """SPEC.md V5a: filter off-source spans through the identical path."""
        ev = args.control_event
        dets = [d for d in args.detectors if (ev, d) in conditioned]
        if len(dets) < 2:
            logger.warning("off-source null test needs two detectors; skipping")
            return []

        span_coincs = []
        for s in range(args.n_offsource):
            trigs = []
            for det in dets:
                for k in range(args.control_chunks):
                    t = File(f"off_{ev}_{det}_s{s:02d}_{k:03d}.json")
                    self.wf.add_jobs(
                        self._job("matched_filter",
                                  f"offfilter_{ev}_{det}_s{s:02d}_{k:03d}",
                                  label="control")
                        .add_args("--conditioned", conditioned[(ev, det)],
                                  "--bank", bank,
                                  "--dataset", f"offsource/span_{s:02d}",
                                  "--chunk", str(k),
                                  "--n-chunks", str(args.control_chunks),
                                  "--snr-threshold", str(args.snr_threshold),
                                  "--output", t)
                        .add_inputs(conditioned[(ev, det)], bank)
                        .add_outputs(t, stage_out=False, register_replica=False)
                    )
                    trigs.append(t)
            sc = File(f"off_{ev}_s{s:02d}_coincidence.json")
            cj = (
                self._job("coincidence", f"offcoinc_{ev}_s{s:02d}", label="control")
                .add_args("--event", f"{ev}_offsource_span{s:02d}",
                          "--coincidence-window-ms", str(args.coincidence_window_ms),
                          "--output", sc)
                .add_inputs(*trigs)
                .add_outputs(sc, stage_out=False, register_replica=False)
            )
            for t in trigs:
                cj.add_args("--triggers", t)
            self.wf.add_jobs(cj)
            span_coincs.append(sc)

        out = File("control_offsource.json")
        j = (
            self._job("control_offsource", "control_offsource", label="control")
            .add_args("--onsource", self._control_onsource, "--event", ev,
                      "--output", out)
            .add_inputs(self._control_onsource, *span_coincs)
            .add_outputs(out, stage_out=True, register_replica=False)
        )
        for sc in span_coincs:
            j.add_args("--offsource", sc)
        self.wf.add_jobs(j)
        return [out]

    def _add_injection_control(self, args, conditioned, bank):
        """SPEC.md V5b: recover known-parameter injections."""
        ev = args.control_event
        det = args.detectors[0]
        if (ev, det) not in conditioned:
            logger.warning("injection control needs a conditioned segment; skipping")
            return []
        outs = []
        for i, mc in enumerate(args.injection_chirp_masses):
            out = File(f"control_injection_{i:02d}.json")
            self.wf.add_jobs(
                self._job("control_injection", f"control_injection_{i:02d}",
                          label="control")
                .add_args("--conditioned", conditioned[(ev, det)], "--bank", bank,
                          "--span", f"offsource/span_{i % args.n_offsource:02d}",
                          "--chirp-mass", str(mc),
                          "--target-snr", str(args.injection_snr),
                          "--crop-length", str(args.injection_crop_length),
                          "--snr-threshold", str(args.snr_threshold),
                          "--output", out)
                .add_inputs(conditioned[(ev, det)], bank)
                .add_outputs(out, stage_out=True, register_replica=False)
            )
            outs.append(out)
        return outs


# ======================================================================
def main():
    p = argparse.ArgumentParser(
        description="GWTC-1 BBH reproduction from GWOSC open data (see SPEC.md)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                   Tier A events, H1+L1
  %(prog)s --tier all                        Tier A + Tier B
  %(prog)s --smoke-test                      one event, coarse bank, no controls
  %(prog)s --events GW150914 --bank-chunks 4
""")

    # --- standard Pegasus arguments ---
    p.add_argument("-s", "--skip-sites-catalog", action="store_true")
    p.add_argument("-e", "--execution-site-name", default="condorpool")
    p.add_argument("-o", "--output", default="workflow.yml")
    p.add_argument("--container-image", default=None,
                   help="Path to GW_Container.sif (default: ./GW_Container.sif)")

    # --- event selection (SPEC.md §1.2) ---
    p.add_argument("--events", nargs="+", default=None,
                   help="Explicit event list, space- or comma-separated "
                        "(overrides --tier)")
    p.add_argument("--tier", choices=["A", "B", "all"], default="A",
                   help="A = the 7 primary events (default); all = A + B")
    p.add_argument("--detectors", nargs="+", default=["H1", "L1"],
                   help="Space- or comma-separated. SPEC.md §1.3: V1 is opt-in "
                        "and contributes little")

    # --- data / conditioning (SPEC.md D2-D5) ---
    p.add_argument("--segment-duration", type=int, default=4096)
    p.add_argument("--source-sample-rate", type=int, default=4096)
    p.add_argument("--sample-rate", type=int, default=2048)
    p.add_argument("--f-low", type=float, default=20.0)
    p.add_argument("--f-high", type=float, default=512.0)
    p.add_argument("--analysis-length", type=float, default=256.0)
    p.add_argument("--psd-length", type=float, default=1024.0)
    p.add_argument("--highpass", type=float, default=15.0,
                   help="High-pass cutoff below the band (SoftwareX §2.2)")
    p.add_argument("--gate-threshold", type=float, default=50.0)

    # --- bank (SPEC.md B1-B5) ---
    p.add_argument("--mc-min", type=float, default=5.0)
    p.add_argument("--mc-max", type=float, default=60.0)
    p.add_argument("--q-max", type=float, default=6.0)
    p.add_argument("--n-q", type=int, default=6)
    p.add_argument("--minimal-match", type=float, default=0.97)
    p.add_argument("--approximant", default="IMRPhenomD",
                   help="SPEC.md M2. SEOBNRv4_ROM is selectable but "
                        "needs its LAL data files (see bin/gwlib.py)")

    # --- search ---
    p.add_argument("--bank-chunks", type=int, default=8,
                   help="Filter fan-out width per (event, detector)")
    p.add_argument("--snr-threshold", type=float, default=5.5)
    p.add_argument("--coincidence-window-ms", type=float, default=15.0)

    # --- controls (SPEC.md V5) ---
    p.add_argument("--n-offsource", type=int, default=10)
    p.add_argument("--control-chunks", type=int, default=4)
    p.add_argument("--control-event", default=None,
                   help="Event whose data the controls use (default: first)")
    # Deliberately NOT type=float: the GUI sends "8,15,25,40,55" as one token,
    # which argparse would reject before split_list() ever sees it.
    p.add_argument("--injection-chirp-masses", nargs="+",
                   default=[8.0, 15.0, 25.0, 40.0, 55.0],
                   help="SPEC.md V5b: 5 masses spanning the bank, space- or "
                        "comma-separated")
    p.add_argument("--injection-snr", type=float, default=15.0)
    p.add_argument("--injection-crop-length", type=float, default=64.0)

    # --- iteration flags (SPEC.md §7.3) ---
    p.add_argument("--skip-controls", action="store_true",
                   help="For iteration only. Forces the report to INCOMPLETE.")
    p.add_argument("--skip-spectrograms", action="store_true",
                   help="For iteration only. Forces the report to INCOMPLETE.")
    p.add_argument("--smoke-test", action="store_true",
                   help="SPEC.md W8: one event, coarse bank, no controls")

    args = p.parse_args()

    # --- normalise list arguments (accept GUI comma form) ---
    try:
        if args.events is not None:
            args.events = split_list(args.events, str, "--events")
        args.detectors = split_list(args.detectors, str, "--detectors")
        args.injection_chirp_masses = split_list(
            args.injection_chirp_masses, float, "--injection-chirp-masses")
    except argparse.ArgumentTypeError as exc:
        p.error(str(exc))

    if not args.detectors:
        p.error("--detectors is empty")
    if not args.injection_chirp_masses:
        p.error("--injection-chirp-masses is empty")

    # --- resolve the event list ---
    if args.events is None:
        args.events = {"A": TIER_A, "B": TIER_B, "all": TIER_A + TIER_B}[args.tier]
    bad = EXCLUDED.intersection(args.events)
    if bad:
        p.error(f"{', '.join(sorted(bad))} is out of scope (SPEC.md N11: BNS)")
    unknown = [e for e in args.events if e not in TIER_A + TIER_B]
    if unknown:
        p.error(f"unknown event(s): {', '.join(unknown)}. "
                f"SPEC.md N12: only O1/O2 GWTC-1 BBHs are in scope.")

    if args.smoke_test:
        logger.info("Smoke test: one event, coarse bank, controls off")
        args.events = args.events[:1]
        args.bank_chunks = min(args.bank_chunks, 2)
        args.n_q = 2
        args.mc_min, args.mc_max = 25.0, 40.0
        args.minimal_match = 0.90
        args.analysis_length = 64.0
        args.psd_length = 256.0
        args.n_offsource = 2
        args.skip_controls = True

    if args.control_event is None:
        args.control_event = args.events[0]
    elif args.control_event not in args.events:
        p.error(f"--control-event {args.control_event} is not in the event list")

    # --- pre-submission validation (before any Pegasus API call) ---
    if len(args.detectors) < 2 and not args.skip_controls:
        p.error("the off-source null test needs at least two detectors; "
                "add a detector or pass --skip-controls")
    if args.f_high > args.sample_rate / 2:
        p.error(f"--f-high {args.f_high} exceeds Nyquist for "
                f"--sample-rate {args.sample_rate}")
    if args.mc_min >= args.mc_max:
        p.error("--mc-min must be below --mc-max")
    if args.highpass and args.highpass >= args.f_low:
        p.error(f"--highpass {args.highpass} must sit below --f-low {args.f_low}")

    wf_dir = Path(__file__).parent.resolve()
    for rel in ["bin/gwlib.py", f"references/{REFERENCE_JSON}"]:
        if not (wf_dir / rel).exists():
            p.error(f"required file missing: {rel}")
    for tool in TOOL_CONFIGS:
        if not (wf_dir / "bin" / f"{tool}.py").exists():
            p.error(f"wrapper missing: bin/{tool}.py")

    logger.info("Events (tier %s): %s", args.tier, ", ".join(args.events))
    logger.info("Detectors: %s | bank chunks: %d | controls: %s",
                ", ".join(args.detectors), args.bank_chunks,
                "off" if args.skip_controls else "on")

    wf = GwWorkflow(dagfile=args.output)
    wf.create_pegasus_properties()
    if not args.skip_sites_catalog:
        wf.create_sites_catalog(args.execution_site_name)
    wf.create_transformation_catalog(args.execution_site_name, args.container_image)
    wf.create_replica_catalog()
    wf.create_workflow(args)
    wf.write()

    n_jobs = len(wf.wf.jobs)
    logger.info("Wrote %s: %d jobs", args.output, n_jobs)
    if args.skip_controls or args.skip_spectrograms:
        logger.warning("controls and/or spectrograms skipped -- the report will "
                       "read INCOMPLETE, never PASS (SPEC.md §7.3)")


if __name__ == "__main__":
    main()
