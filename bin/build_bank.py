#!/usr/bin/env python3
"""Build the non-spinning template bank, once, for the whole workflow.

SPEC.md B1-B7. Usman §3.1: template density depends on the detectors' noise PSD,
and the pipeline places one bank from a PSD averaged over time AND over
detectors via the harmonic mean (eq. 7), so a single bank serves every detector
-- which is what makes the exact-match coincidence of §3.5 possible.

NOTE -- deviation from SPEC.md §7.2: the spec's DAG sketch has build_bank
independent of conditioning. Faithfulness to Usman §3.1 requires the real noise
curve, so this job consumes every conditioned PSD and is a fan-in over
condition_strain. Placement against an assumed analytic PSD would make the
measured minimal match meaningless.

Placement: for each mass ratio, walk chirp mass upward and choose each next
template by bisection so the MEASURED overlap with its neighbour equals the
minimal match. The bank therefore reports a measured minimal match rather than
asserting one (SPEC.md B4, V2d).

Placement and the B6 gate live in gwlib, shared with the parallel path
(bank_psd -> bank_row x n_q -> merge_bank), so the two cannot diverge. This
serial entry point remains for run_manual.sh and hand runs; the workflow uses
the fan-out, which turns a ~29 min serial build into ~5 min wall clock.
"""

import argparse
import json
import logging
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.getcwd())
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser(description="Build the template bank (SPEC.md B1-B7)")
    p.add_argument("--psd", action="append", required=True,
                   help="Conditioned HDF5 to draw a PSD from (repeatable)")
    p.add_argument("--reference", required=True,
                   help="gwtc1_reference.json, for the B6 adequacy gate")
    p.add_argument("--tier-a", nargs="+", required=True,
                   help="Tier A event names the bank must cover, space- or "
                        "comma-separated")
    p.add_argument("--mc-min", type=float, default=5.0, help="SPEC B1")
    p.add_argument("--mc-max", type=float, default=60.0, help="SPEC B1")
    p.add_argument("--q-max", type=float, default=6.0, help="SPEC B2")
    p.add_argument("--n-q", type=int, default=6, help="Mass-ratio rows (log spaced)")
    p.add_argument("--minimal-match", type=float, default=0.97, help="SPEC B4")
    p.add_argument("--f-low", type=float, default=20.0)
    p.add_argument("--f-high", type=float, default=512.0)
    p.add_argument("--approximant", default="IMRPhenomD", help="SPEC M2")
    p.add_argument("--adequacy-tol", type=float, default=0.01,
                   help="SPEC B6: template within this fractional Mc of each event")
    p.add_argument("--allow-fallback-waveform", action="store_true",
                   help="Permit the restricted-PN fallback when the requested "
                        "approximant is unavailable. FOR SMOKE TESTING ONLY: "
                        "downstream validation still refuses to report PASS "
                        "from a fallback bank (SPEC.md M2).")
    p.add_argument("--max-templates", type=int, default=20000,
                   help="Hard stop so a bad spacing cannot run away")
    p.add_argument("--output", required=True, help="Output bank HDF5")
    p.add_argument("--manifest", required=True, help="Output manifest JSON")
    args = p.parse_args()

    # Accept the GUI's comma-separated list form as well as space-separated.
    args.tier_a = [e for raw in args.tier_a
                   for e in str(raw).replace(",", " ").split()]
    if not args.tier_a:
        p.error("--tier-a is empty")

    # Validate inputs before touching h5py, so a typo'd path gives an
    # actionable message instead of a raw traceback.
    missing = [q for q in list(args.psd) + [args.reference] if not os.path.exists(q)]
    if missing:
        p.error("input file(s) not found: " + ", ".join(missing))

    for path in (args.output, args.manifest):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    # --- harmonic-mean PSD over every conditioned detector-segment ---------
    # Analysis segments may differ in length -- check_quality shortens a window
    # that cannot fit inside a clean science-mode run -- so PSDs are interpolated
    # onto a common grid rather than required to share a resolution.
    raw, sources = [], []
    for path in args.psd:
        with h5py.File(path, "r") as f:
            raw.append((f["psd/freq"][:], f["psd/psd_truncated"][:]))
            sources.append({"file": os.path.basename(path),
                            "event": f["meta"].attrs["event"],
                            "detector": f["meta"].attrs["detector"],
                            "psd_bins": int(len(raw[-1][0]))})
    # Finest grid wins, so no segment's resolution is thrown away.
    freqs = max((fr for fr, _ in raw), key=len)
    psds = [ps if len(fr) == len(freqs)
            else gwlib.interpolate_psd(fr, ps, freqs) for fr, ps in raw]
    if len({s["psd_bins"] for s in sources}) > 1:
        logger.info("PSD bin counts differ across segments (%s); interpolated "
                    "onto the finest grid (%d bins)",
                    sorted({s["psd_bins"] for s in sources}), len(freqs))
    psd = gwlib.harmonic_mean_psd(psds)
    df = float(freqs[1] - freqs[0])
    mask = gwlib.band_mask(freqs, args.f_low, args.f_high)
    logger.info("Harmonic-mean PSD from %d segment(s), df=%.6f Hz, %d band bins",
                len(psds), df, int(mask.sum()))

    # Probe the approximant ONCE before placing thousands of templates. Without
    # this, a missing SEOBNRv4_ROM data file degrades every template to the PN
    # toy and the failure only surfaces as a verdict at the very end of the run.
    _probe_h, _probe_fb = gwlib.fd_waveform(38.6, 33.3, args.f_low, df,
                                            args.f_high, args.approximant)
    if _probe_fb and not args.allow_fallback_waveform:
        logger.error(
            "approximant %r is unavailable, so every template would come from "
            "the restricted-PN fallback and no result could be validated "
            "(SPEC.md M2). If this is SEOBNRv4_ROM, LALSuite >= 7.25 needs its "
            "interpolation grids in $LAL_DATA_PATH -- see "
            "https://zenodo.org/records/14999310. Otherwise install PyCBC, or "
            "pass --allow-fallback-waveform to build a smoke-test bank whose "
            "results can never be validated.",
            args.approximant)
        sys.exit(1)
    if _probe_fb:
        logger.warning("approximant %r unavailable; building a FALLBACK bank "
                       "because --allow-fallback-waveform was given. No result "
                       "from this bank can be validated (SPEC.md M2).",
                       args.approximant)
    else:
        logger.info("Approximant %s available (%d bins probed)",
                    args.approximant, len(_probe_h))

    # --- placement (shared with the parallel bank_row path via gwlib) ------
    q_grid = np.geomspace(1.0, args.q_max, args.n_q)
    templates, adjacent, used_fallback, truncated = [], [], False, False
    for q in q_grid:
        rt, ra, rfb, rtr = gwlib.place_bank_row(
            q, args.mc_min, args.mc_max, args.minimal_match, psd, freqs, df,
            mask, args.f_low, args.f_high, args.approximant,
            args.max_templates - len(templates),
        )
        templates.extend(rt)
        adjacent.extend(ra)
        used_fallback = used_fallback or rfb
        truncated = truncated or rtr
        logger.info("q=%.2f: %d templates", q, len(rt))
        if len(templates) >= args.max_templates:
            truncated = True
            break

    if truncated:
        logger.warning("hit --max-templates=%d; bank is truncated (coverage "
                       "incomplete, reported in the manifest)", args.max_templates)

    mcs = np.array([t[0] for t in templates])
    logger.info("Bank: %d templates, Mc in [%.2f, %.2f]", len(templates),
                mcs.min(), mcs.max())

    # --- SPEC.md B6 adequacy gate -----------------------------------------
    with open(args.reference) as fh:
        ref = json.load(fh)
    adequacy, gate_failures = gwlib.bank_adequacy(mcs, ref["events"],
                                                  args.tier_a, args.adequacy_tol)

    manifest = {
        "n_templates": len(templates),
        "chirp_mass_detector_range": [float(mcs.min()), float(mcs.max())],
        "requested_mc_range": [args.mc_min, args.mc_max],
        "q_grid": [round(float(q), 4) for q in q_grid],
        "spins": "zero (SPEC.md M3, non-spinning bank)",
        "approximant": args.approximant,
        "waveform_fallback_used": bool(used_fallback),
        "minimal_match_requested": args.minimal_match,
        "minimal_match_measured_min": (round(float(np.min(adjacent)), 5)
                                       if adjacent else None),
        "minimal_match_measured_median": (round(float(np.median(adjacent)), 5)
                                          if adjacent else None),
        "band_hz": [args.f_low, args.f_high],
        "df_hz": df,
        "psd_sources": sources,
        "psd_average": "harmonic mean, Usman 2016 eq. (7)",
        "truncated": bool(truncated),
        "built_by": "build_bank (serial)",
        "adequacy_gate": {"tolerance": args.adequacy_tol,
                          "events": adequacy,
                          "failures": gate_failures,
                          "pass": not gate_failures},
    }
    if used_fallback:
        manifest["WARNING"] = (
            "Templates came from the restricted-PN fallback, not the requested "
            "approximant. Adequate for a smoke test only; validation must not "
            "report PASS from this bank (SPEC.md M2)."
        )

    with h5py.File(args.output, "w") as f:
        arr = np.array(templates, dtype=np.float64)
        d = f.create_dataset("templates", data=arr)
        d.attrs["columns"] = "chirp_mass_detector,q,mass1,mass2"
        g = f.create_group("psd")
        g.create_dataset("freq", data=freqs)
        g.create_dataset("psd", data=psd)
        f["/"].attrs["approximant"] = args.approximant
        f["/"].attrs["f_low"] = args.f_low
        f["/"].attrs["f_high"] = args.f_high
        f["/"].attrs["minimal_match"] = args.minimal_match
        f["/"].attrs["waveform_fallback_used"] = bool(used_fallback)

    with open(args.manifest, "w") as fh:
        json.dump(manifest, fh, indent=2)

    for msg in gate_failures:
        logger.error("B6 adequacy failure: %s", msg)
    logger.info("Bank: %s (%d templates), manifest: %s",
                args.output, len(templates), args.manifest)

    # SPEC.md B6: the adequacy gate gates the search. Fail loudly and stop the
    # DAG rather than filtering with a bank that cannot represent the events.
    if gate_failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
