#!/usr/bin/env python3
"""Compute the bank's placement PSD once, from every conditioned segment.

Usman §3.1: the pipeline places one bank from a noise PSD averaged over time AND
over detectors via the harmonic mean (eq. 7), so a single bank serves every
detector -- which is what makes the exact-match coincidence of §3.5 possible.

Split out from the bank build so the fan-out over mass-ratio rows does not
restage every conditioned segment. The conditioned files are ~30 MB each; with 14
segments and 6 rows, recomputing the mean per row would move ~2.5 GB instead of
the ~2 MB this job emits.
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
    p = argparse.ArgumentParser(description="Harmonic-mean placement PSD (Usman eq. 7)")
    p.add_argument("--psd", action="append", required=True,
                   help="Conditioned HDF5 to draw a PSD from (repeatable)")
    p.add_argument("--output", required=True, help="Output PSD HDF5")
    p.add_argument("--report", required=True, help="Output provenance JSON")
    args = p.parse_args()

    missing = [q for q in args.psd if not os.path.exists(q)]
    if missing:
        p.error("input file(s) not found: " + ", ".join(missing))
    for path in (args.output, args.report):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    raw, sources = [], []
    for path in args.psd:
        with h5py.File(path, "r") as f:
            fr = f["psd/freq"][:]
            ps = f["psd/psd_truncated"][:]
            sources.append({"file": os.path.basename(path),
                            "event": f["meta"].attrs["event"],
                            "detector": f["meta"].attrs["detector"],
                            "psd_bins": int(len(fr))})
        raw.append((fr, ps))

    # Analysis segments may differ in length -- check_quality shortens a window
    # that cannot fit inside a clean science-mode run -- so interpolate onto the
    # finest grid rather than demanding a shared resolution.
    freqs = max((fr for fr, _ in raw), key=len)
    psds = [ps if len(fr) == len(freqs)
            else gwlib.interpolate_psd(fr, ps, freqs) for fr, ps in raw]
    bins = sorted({s["psd_bins"] for s in sources})
    if len(bins) > 1:
        logger.info("PSD bin counts differ across segments (%s); interpolated "
                    "onto the finest grid (%d bins)", bins, len(freqs))

    psd = gwlib.harmonic_mean_psd(psds)
    df = float(freqs[1] - freqs[0])
    logger.info("Harmonic-mean PSD from %d segment(s), df=%.6f Hz, %d bins",
                len(psds), df, len(freqs))

    with h5py.File(args.output, "w") as f:
        f.create_dataset("freq", data=freqs)
        f.create_dataset("psd", data=psd)
        f["/"].attrs["n_segments"] = len(psds)
        f["/"].attrs["df"] = df
        f["/"].attrs["average"] = "harmonic mean, Usman 2016 eq. (7)"

    with open(args.report, "w") as fh:
        json.dump({"n_segments": len(psds), "df_hz": df, "n_bins": int(len(freqs)),
                   "psd_average": "harmonic mean, Usman 2016 eq. (7)",
                   "interpolated_to_finest_grid": len(bins) > 1,
                   "segment_bin_counts": bins, "sources": sources}, fh, indent=2)
    logger.info("Output: %s, %s", args.output, args.report)


if __name__ == "__main__":
    main()
