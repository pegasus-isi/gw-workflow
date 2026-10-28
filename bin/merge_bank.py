#!/usr/bin/env python3
"""Concatenate the per-q bank rows and apply the SPEC.md B6 adequacy gate.

Fan-in over bank_row. Rows are independent by construction (Usman §3.1 places at
fixed mass ratio), so merging is concatenation -- but the gate that decides
whether the search may proceed can only be applied to the whole bank, which is
why it lives here and not in the row jobs.

Emits the same bank.hdf5 / bank_manifest.json contract as the serial
build_bank.py, so matched_filter and aggregate_report are indifferent to which
path built the bank.
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
    p = argparse.ArgumentParser(description="Merge bank rows, gate on B6 adequacy")
    p.add_argument("--row", action="append", required=True,
                   help="Row JSON from bank_row (repeatable)")
    p.add_argument("--psd", required=True, help="PSD HDF5 from bank_psd")
    p.add_argument("--psd-report", required=True,
                   help="Provenance JSON from bank_psd")
    p.add_argument("--reference", required=True,
                   help="gwtc1_reference.json, for the B6 adequacy gate")
    p.add_argument("--tier-a", nargs="+", required=True,
                   help="Tier A event names the bank must cover, space- or "
                        "comma-separated")
    p.add_argument("--adequacy-tol", type=float, default=0.01,
                   help="SPEC B6: template within this fractional Mc of each event")
    p.add_argument("--output", required=True, help="Output bank HDF5")
    p.add_argument("--manifest", required=True, help="Output manifest JSON")
    args = p.parse_args()

    args.tier_a = [e for raw in args.tier_a
                   for e in str(raw).replace(",", " ").split()]
    if not args.tier_a:
        p.error("--tier-a is empty")

    inputs = list(args.row) + [args.psd, args.psd_report, args.reference]
    missing = [q for q in inputs if not os.path.exists(q)]
    if missing:
        p.error("input file(s) not found: " + ", ".join(missing))
    for path in (args.output, args.manifest):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    rows = []
    for path in args.row:
        with open(path) as fh:
            rows.append(json.load(fh))
    rows.sort(key=lambda r: r["row"])

    # A row that silently disagreed about the band or the minimal match would
    # produce a bank with no single measured minimal match to report.
    for key in ("band_hz", "minimal_match_requested", "approximant",
                "mc_range_requested"):
        vals = {json.dumps(r[key]) for r in rows}
        if len(vals) > 1:
            logger.error("rows disagree on %s: %s -- they were not placed "
                         "against the same bank definition", key, sorted(vals))
            sys.exit(1)

    templates, adjacent = [], []
    for r in rows:
        templates.extend([tuple(t) for t in r["templates"]])
        adjacent.extend(r["adjacent_overlaps"])

    if not templates:
        logger.error("merged bank is empty -- every row placed zero templates")
        sys.exit(1)

    used_fallback = any(r["waveform_fallback_used"] for r in rows)
    truncated_rows = [r["row"] for r in rows if r["truncated"]]
    q_grid = [r["q"] for r in rows]
    approximant = rows[0]["approximant"]
    f_low, f_high = rows[0]["band_hz"]
    mc_min, mc_max = rows[0]["mc_range_requested"]
    minimal_match = rows[0]["minimal_match_requested"]
    df = rows[0]["df_hz"]

    with h5py.File(args.psd, "r") as f:
        freqs = f["freq"][:]
        psd = f["psd"][:]
    with open(args.psd_report) as fh:
        psd_report = json.load(fh)

    mcs = np.array([t[0] for t in templates])
    logger.info("Merged %d row(s): %d templates, Mc in [%.2f, %.2f]",
                len(rows), len(templates), mcs.min(), mcs.max())
    for r in rows:
        logger.info("  row %d q=%.4f: %d templates, %.1fs",
                    r["row"], r["q"], r["n_templates"], r["build_seconds"])

    # --- SPEC.md B6 adequacy gate -----------------------------------------
    with open(args.reference) as fh:
        ref = json.load(fh)
    adequacy, gate_failures = gwlib.bank_adequacy(mcs, ref["events"],
                                                  args.tier_a, args.adequacy_tol)

    manifest = {
        "n_templates": len(templates),
        "chirp_mass_detector_range": [float(mcs.min()), float(mcs.max())],
        "requested_mc_range": [mc_min, mc_max],
        "q_grid": [round(float(q), 4) for q in q_grid],
        "spins": "zero (SPEC.md M3, non-spinning bank)",
        "approximant": approximant,
        "waveform_fallback_used": bool(used_fallback),
        "minimal_match_requested": minimal_match,
        "minimal_match_measured_min": (round(float(np.min(adjacent)), 5)
                                       if adjacent else None),
        "minimal_match_measured_median": (round(float(np.median(adjacent)), 5)
                                          if adjacent else None),
        "band_hz": [f_low, f_high],
        "df_hz": df,
        "psd_sources": psd_report["sources"],
        "psd_average": psd_report["psd_average"],
        "truncated": bool(truncated_rows),
        "truncated_rows": truncated_rows,
        "built_by": "merge_bank (parallel row fan-out)",
        "rows": [{"row": r["row"], "q": round(r["q"], 4),
                  "n_templates": r["n_templates"],
                  "build_seconds": r["build_seconds"]} for r in rows],
        "build_seconds_total_cpu": round(sum(r["build_seconds"] for r in rows), 2),
        "build_seconds_slowest_row": round(max(r["build_seconds"] for r in rows), 2),
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
        d = f.create_dataset("templates", data=np.array(templates, dtype=np.float64))
        d.attrs["columns"] = "chirp_mass_detector,q,mass1,mass2"
        g = f.create_group("psd")
        g.create_dataset("freq", data=freqs)
        g.create_dataset("psd", data=psd)
        f["/"].attrs["approximant"] = approximant
        f["/"].attrs["f_low"] = f_low
        f["/"].attrs["f_high"] = f_high
        f["/"].attrs["minimal_match"] = minimal_match
        f["/"].attrs["waveform_fallback_used"] = bool(used_fallback)

    with open(args.manifest, "w") as fh:
        json.dump(manifest, fh, indent=2)

    for msg in gate_failures:
        logger.error("B6 adequacy failure: %s", msg)
    logger.info("Bank: %s (%d templates), manifest: %s",
                args.output, len(templates), args.manifest)

    # SPEC.md B6: the gate gates the search. Stop the DAG rather than filter
    # with a bank that cannot represent the events.
    if gate_failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
