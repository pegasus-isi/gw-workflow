#!/usr/bin/env python3
"""Matched-filter one conditioned segment against one chunk of the bank.

SPEC.md M2, M3, M7, M11. Implements Usman §2 eq. (1)-(2) and the §3.3
threshold-and-cluster step. Fanned out over (event, detector, bank chunk) so the
16-core pool stays saturated (SPEC.md W2).

SPEC.md M11: templates are generated one at a time inside the loop and never
cached. A frequency-domain template at 2048 Hz over a 256 s segment is ~8 MB as
complex128, so caching even a hundred would eat most of the 1.5 GB job budget.

SPEC.md M5: the first and last 16 s of the segment are filter-corrupted (the
inverse-PSD truncation length) and triggers there are discarded.
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


def load_conditioned(path, dataset="analysis"):
    """Read one segment. `dataset` is "analysis" or "offsource/span_NN".

    The off-source null test (SPEC.md V5a) runs through this same function and
    the same filter_segment() below, so "identical configuration" is guaranteed
    by construction rather than by two implementations agreeing.
    """
    with h5py.File(path, "r") as f:
        if dataset == "analysis":
            strain = f["analysis/strain"][:]
            gps_start = float(f["analysis/strain"].parent.attrs["gps_start"])
        else:
            d = f[dataset]
            strain = d[:]
            gps_start = float(d.attrs["gps_start"])
        fs = int(f["analysis/strain"].parent.attrs["sample_rate"])
        event_gps = float(f["analysis/strain"].parent.attrs["event_gps"])
        psd = f["psd/psd_truncated"][:]
        psd_freq = f["psd/freq"][:]
        meta = dict(f["meta"].attrs)
    return strain, gps_start, fs, event_gps, psd, psd_freq, meta


def filter_segment(strain, fs, psd_full, psd_freq_full, templates, f_low, f_high,
                   threshold, approximant, edge_guard, gps_start, progress_every=50):
    """Filter one segment against a template list; return triggers."""
    n = len(strain)
    dt = 1.0 / fs
    df = 1.0 / (n * dt)
    f = gwlib.freqs(n, dt)
    # The conditioned PSD is stored at the analysis segment's resolution; an
    # off-source span of a different length needs it re-interpolated.
    psd = (psd_full if len(psd_full) == len(f)
           else gwlib.interpolate_psd(psd_freq_full, psd_full, f))
    mask = gwlib.band_mask(f, f_low, f_high)
    data_fd = gwlib.rfft_norm(strain, dt)

    guard = int(edge_guard * fs)
    if n - 2 * guard <= 0:
        raise RuntimeError(
            f"segment of {n / fs:.0f}s leaves no usable samples after "
            f"discarding {edge_guard:.0f}s at each end as filter-corrupted "
            f"(SPEC.md M5); every trigger would be dropped silently")
    triggers, fallback_any, t_start = [], False, time.time()

    for i, (mc, q, m1, m2) in enumerate(templates):
        h, fb = gwlib.fd_waveform(m1, m2, f_low, df, f_high, approximant)
        fallback_any |= fb
        if len(h) < len(f):
            h = np.pad(h, (0, len(f) - len(h)))
        h = h[: len(f)]
        snr = gwlib.matched_filter_snr(data_fd, h, psd, df, n, mask)
        for idx, val in gwlib.cluster_triggers(snr, fs, threshold):
            if idx < guard or idx > n - guard:
                continue           # SPEC.md M5: filter-corrupted edges
            triggers.append({
                "template_index": i, "chirp_mass_detector": mc, "q": q,
                "mass1": m1, "mass2": m2, "snr": float(val),
                "gps_time": gps_start + idx * dt, "sample_index": int(idx),
            })
        if progress_every and (i + 1) % progress_every == 0:
            rate = (i + 1) / max(1e-9, time.time() - t_start)
            logger.info("  %d/%d templates (%.1f/s), %d triggers",
                        i + 1, len(templates), rate, len(triggers))
    return triggers, fallback_any


def main():
    p = argparse.ArgumentParser(description="Matched filter one bank chunk")
    p.add_argument("--conditioned", required=True, help="Conditioned HDF5")
    p.add_argument("--dataset", default="analysis",
                   help='Segment to filter: "analysis" (default) or '
                        '"offsource/span_NN" for the SPEC.md V5a null test')
    p.add_argument("--bank", required=True, help="Bank HDF5 from build_bank")
    p.add_argument("--chunk", type=int, required=True, help="0-based chunk index")
    p.add_argument("--n-chunks", type=int, required=True, help="Total chunks")
    p.add_argument("--snr-threshold", type=float, default=5.5,
                   help="Single-detector SNR threshold (GWTC-1 §III A uses 5.5)")
    p.add_argument("--output", required=True, help="Output triggers JSON")
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)
    if not 0 <= args.chunk < args.n_chunks:
        logger.error("chunk %d out of range for %d chunks", args.chunk, args.n_chunks)
        sys.exit(2)

    strain, gps_start, fs, event_gps, psd, psd_freq, meta = load_conditioned(
        args.conditioned, args.dataset)
    f_low = float(meta.get("f_low", 20.0))
    f_high = float(meta.get("f_high", 512.0))

    with h5py.File(args.bank, "r") as f:
        all_tmpl = f["templates"][:]
        approximant = f["/"].attrs.get("approximant", "IMRPhenomD")
        bank_fallback = bool(f["/"].attrs.get("waveform_fallback_used", False))

    idx = np.array_split(np.arange(len(all_tmpl)), args.n_chunks)[args.chunk]
    chunk = [tuple(all_tmpl[i]) for i in idx]
    logger.info("%s %s [%s]: chunk %d/%d -> %d of %d templates",
                meta.get("event"), meta.get("detector"), args.dataset,
                args.chunk, args.n_chunks, len(chunk), len(all_tmpl))
    logger.info("Segment %d samples @ %d Hz, band %.0f-%.0f Hz, threshold %.2f",
                len(strain), fs, f_low, f_high, args.snr_threshold)

    triggers, fallback = filter_segment(
        strain, fs, psd, psd_freq, chunk, f_low, f_high, args.snr_threshold,
        approximant, gwlib.INV_PSD_TRUNC_S, gps_start,
    )
    # Global template indices, so coincidence can require the SAME template
    # across detectors (Usman §3.5 exact-match coincidence).
    for t in triggers:
        t["template_index"] = int(idx[t["template_index"]])

    out = {
        "event": meta.get("event"), "detector": meta.get("detector"),
        "dataset": args.dataset,
        "is_offsource": args.dataset != "analysis",
        "event_gps": event_gps, "segment_gps_start": gps_start,
        "sample_rate": fs, "band_hz": [f_low, f_high],
        "chunk": args.chunk, "n_chunks": args.n_chunks,
        "n_templates_in_chunk": len(chunk), "n_templates_total": int(len(all_tmpl)),
        "snr_threshold": args.snr_threshold,
        "approximant": approximant,
        "waveform_fallback_used": bool(fallback or bank_fallback),
        "edge_guard_s": gwlib.INV_PSD_TRUNC_S,
        "n_triggers": len(triggers),
        "triggers": sorted(triggers, key=lambda t: -t["snr"])[:2000],
    }
    if len(triggers) > 2000:
        # CLAUDE.md §2 / SPEC.md: never truncate silently.
        out["truncated"] = True
        out["truncation_note"] = (
            f"kept the 2000 loudest of {len(triggers)} triggers in this chunk"
        )
        logger.warning("kept 2000 loudest of %d triggers", len(triggers))

    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    logger.info("Output: %s (%d triggers)", args.output, len(triggers))


if __name__ == "__main__":
    main()
