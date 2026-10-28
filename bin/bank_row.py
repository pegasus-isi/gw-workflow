#!/usr/bin/env python3
"""Place one constant-mass-ratio row of the template bank.

SPEC.md B1-B5. Rows at fixed q are independent -- the chirp-mass walk at one q
never consults another -- so the bank build fans out across q and merge_bank
concatenates. The placement itself lives in gwlib.place_bank_row so the serial
build_bank path and this parallel path cannot diverge.

Each next template is chosen by bisection so the MEASURED overlap with its
neighbour equals the minimal match, which is why the manifest can report a
measured minimal match rather than asserting one (SPEC.md B4, V2d).
"""

import argparse
import json
import logging
import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(0, os.getcwd())
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser(description="Place one q row of the bank")
    p.add_argument("--psd", required=True, help="PSD HDF5 from bank_psd")
    p.add_argument("--q", type=float, required=True, help="Mass ratio for this row")
    p.add_argument("--row", type=int, required=True, help="Row index (for labelling)")
    p.add_argument("--mc-min", type=float, default=5.0)
    p.add_argument("--mc-max", type=float, default=60.0)
    p.add_argument("--minimal-match", type=float, default=0.97)
    p.add_argument("--f-low", type=float, default=20.0)
    p.add_argument("--f-high", type=float, default=512.0)
    p.add_argument("--approximant", default="IMRPhenomD")
    p.add_argument("--allow-fallback-waveform", action="store_true",
                   help="Permit the restricted-PN fallback. SMOKE TESTING ONLY: "
                        "validation still refuses to PASS from a fallback bank.")
    p.add_argument("--max-templates", type=int, default=20000,
                   help="Per-row stop so a collapsed spacing cannot run away")
    p.add_argument("--output", required=True, help="Output row JSON")
    args = p.parse_args()

    if not os.path.exists(args.psd):
        p.error(f"input file not found: {args.psd}")
    if args.mc_min >= args.mc_max:
        p.error("--mc-min must be below --mc-max")
    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    with h5py.File(args.psd, "r") as f:
        freqs = f["freq"][:]
        psd = f["psd"][:]
        df = float(f["/"].attrs["df"])
    mask = gwlib.band_mask(freqs, args.f_low, args.f_high)

    # Probe once: a missing approximant would otherwise place the whole row from
    # the PN toy and only surface as a verdict at the end of the run.
    _, fb = gwlib.fd_waveform(38.6, 33.3, args.f_low, df, args.f_high,
                              args.approximant)
    if fb and not args.allow_fallback_waveform:
        logger.error(
            "approximant %r is unavailable, so this row would be placed from the "
            "restricted-PN fallback and no result could be validated (SPEC.md M2). "
            "If this is SEOBNRv4_ROM, LALSuite >= 7.25 needs its interpolation "
            "grids in $LAL_DATA_PATH -- see https://zenodo.org/records/14999310. "
            "Otherwise install PyCBC, or pass --allow-fallback-waveform for a "
            "smoke-test bank whose results can never be validated.",
            args.approximant)
        sys.exit(1)
    if fb:
        logger.warning("approximant %r unavailable; placing a FALLBACK row "
                       "because --allow-fallback-waveform was given",
                       args.approximant)

    logger.info("row %d: q=%.4f, Mc %.2f-%.2f, minimal match %.4f",
                args.row, args.q, args.mc_min, args.mc_max, args.minimal_match)
    t0 = time.time()
    templates, adjacent, used_fb, truncated = gwlib.place_bank_row(
        args.q, args.mc_min, args.mc_max, args.minimal_match, psd, freqs, df,
        mask, args.f_low, args.f_high, args.approximant, args.max_templates,
        progress=50,
    )
    elapsed = time.time() - t0
    logger.info("row %d: %d templates in %.1fs (%.2f s/template)",
                args.row, len(templates), elapsed,
                elapsed / max(1, len(templates)))
    if truncated:
        logger.warning("row %d hit --max-templates=%d; this row is incomplete "
                       "and merge_bank will report it", args.row, args.max_templates)

    with open(args.output, "w") as fh:
        json.dump({
            "row": args.row, "q": args.q,
            "mc_range_requested": [args.mc_min, args.mc_max],
            "minimal_match_requested": args.minimal_match,
            "approximant": args.approximant,
            "waveform_fallback_used": bool(used_fb or fb),
            "truncated": bool(truncated),
            "n_templates": len(templates),
            "build_seconds": round(elapsed, 2),
            "band_hz": [args.f_low, args.f_high],
            "df_hz": df,
            # columns: chirp_mass_detector, q, mass1, mass2
            "templates": templates,
            "adjacent_overlaps": [round(float(a), 6) for a in adjacent],
        }, fh)
    logger.info("Output: %s", args.output)


if __name__ == "__main__":
    main()
