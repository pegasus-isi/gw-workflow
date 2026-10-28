#!/usr/bin/env python3
"""Condition one detector's strain: decimate, estimate PSD, gate.

SPEC.md D4-D5, M4-M6. Implements the conditioning chain of Usman §3.2-3.3:

  1. Decimate 4096 -> 2048 Hz (SPEC.md D4). The highest-mass event in scope
     merges near 48 Hz and the lowest-mass near 221 Hz, so a 1024 Hz Nyquist has
     ample margin and halves the FFT cost against 4096 Hz.
  2. Estimate the PSD from an off-source stretch that EXCLUDES +/-16 s around the
     event (SPEC.md M4) -- otherwise the signal biases its own PSD.
  3. Truncate the inverse PSD to 16 s in the time domain (Usman §3.2).
  4. Gate short, loud transients with a Tukey window (Usman §3.2).

Also carves out off-source spans for the SPEC.md V5a null test and the V5b
injection test, so those controls reuse this download instead of re-fetching
(CLAUDE.md §4).

Output HDF5 layout:
  analysis/strain, analysis/gps_start      conditioned analysis segment
  psd/freq, psd/psd, psd/psd_truncated     PSD and its inverse-truncated form
  offsource/span_NN                        off-source spans, same length
  meta/*                                   provenance and span bookkeeping
"""

import argparse
import json
import logging
import os
import sys

import h5py
import numpy as np
from scipy import signal

sys.path.insert(0, os.getcwd())      # support-file import (PEGASUS.md)
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

EDGE_GUARD_S = gwlib.INV_PSD_TRUNC_S     # SPEC.md M5: discard 16 s per end


def allocate_spans(clean_runs, analysis_span, event_gps, psd_s, offsource_s,
                   n_offsource, exclude_s, edge_guard=EDGE_GUARD_S):
    """Place the PSD stretch and off-source spans inside CLEAN data.

    SPEC.md M4 keeps the PSD off the event; V5a keeps off-source spans >= 200 s
    from it. Both must additionally lie in science-mode data: a span straddling
    a science-mode gap is NaN-filled, and after zero-filling it looks like a
    giant glitch -- which would poison the PSD and fake a loud off-source
    trigger, defeating the null test it is supposed to support.

    Candidate regions are the clean runs minus the event exclusion zone, shrunk
    by `edge_guard` at each run boundary.
    """
    regions = []
    for a, b in clean_runs:
        a, b = a + edge_guard, b - edge_guard
        if b - a <= 0:
            continue
        # Remove [event - exclude, event + exclude] from this run.
        for lo, hi in ((a, min(b, event_gps - exclude_s)),
                       (max(a, event_gps + exclude_s), b)):
            if hi - lo > 0:
                regions.append((lo, hi))
    regions.sort(key=lambda r: r[1] - r[0], reverse=True)

    if not regions:
        raise RuntimeError(
            f"no clean data at least {exclude_s:.0f}s from the event is available "
            f"for PSD estimation or off-source spans")

    longest = regions[0][1] - regions[0][0]
    if longest < psd_s:
        raise RuntimeError(
            f"longest clean off-source region is {longest:.0f}s, short of the "
            f"{psd_s:.0f}s requested for PSD estimation (--psd-length). Reduce "
            f"--psd-length or accept a noisier PSD.")
    psd_start = regions[0][0]

    spans = []
    for a, b in sorted(regions):
        t = a
        while t + offsource_s <= b and len(spans) < n_offsource:
            spans.append(t)
            t += offsource_s
    return psd_start, spans


def main():
    p = argparse.ArgumentParser(description="Condition strain for matched filtering")
    p.add_argument("--strain", required=True, help="GWOSC HDF5 from fetch_strain")
    p.add_argument("--meta", required=True, help="Provenance JSON from fetch_strain")
    p.add_argument("--quality", required=True, help="Verdict JSON from check_quality")
    p.add_argument("--sample-rate", type=int, default=2048, help="Target rate (SPEC D4)")
    p.add_argument("--f-low", type=float, default=20.0, help="Analysis band low (SPEC D5)")
    p.add_argument("--f-high", type=float, default=512.0, help="Analysis band high (SPEC D5)")
    p.add_argument("--analysis-length", type=float, default=256.0,
                   help="Requested length; the ACTUAL span comes from the "
                        "quality verdict, which places it inside clean data")
    p.add_argument("--psd-length", type=float, default=1024.0)
    p.add_argument("--offsource-length", type=float, default=256.0)
    p.add_argument("--n-offsource", type=int, default=10, help="SPEC V5a wants 10")
    p.add_argument("--exclude", type=float, default=200.0,
                   help="Keep off-source spans this far from the event (SPEC V5a)")
    p.add_argument("--highpass", type=float, default=15.0,
                   help="High-pass cutoff in Hz, below the analysis band "
                        "(SoftwareX §2.2). 0 disables it.")
    p.add_argument("--gate-threshold", type=float, default=50.0,
                   help="Whitened-strain gate threshold. Usman §3.2 fixes the "
                        "window but not this value; SPEC.md §9.3 records it as a "
                        "tuning choice validated against the off-source control.")
    p.add_argument("--output", required=True, help="Output conditioned HDF5")
    p.add_argument("--report", required=True, help="Output conditioning report JSON")
    args = p.parse_args()

    for path in (args.output, args.report):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    with open(args.meta) as fh:
        meta = json.load(fh)
    with open(args.quality) as fh:
        quality = json.load(fh)
    if not quality.get("pass"):
        logger.error("data-quality gate failed for %s %s: %s",
                     meta["event"], meta["detector"], quality.get("failures"))
        sys.exit(1)

    # Single source of truth: check_quality already placed the analysis window
    # inside a clean run, shifting it off-centre if the event sat near a science
    # -mode boundary. Do not recompute it here.
    a_span = quality.get("analysis_span")
    clean = [(r["start"], r["end"]) for r in quality["checks"]["clean_runs"]]
    if not a_span:
        logger.error("quality verdict carries no analysis_span; re-run check_quality")
        sys.exit(1)
    a_start = float(a_span[0])
    analysis_len = float(quality.get("analysis_length_s", a_span[1] - a_span[0]))

    event_gps = float(meta["event_gps"])

    with h5py.File(args.strain, "r") as f:
        raw = f["strain/Strain"][:]
        t0 = float(f["meta/GPSstart"][()])
        file_len = float(f["meta/Duration"][()])
    fs_in = int(round(len(raw) / file_len))
    logger.info("Loaded %d samples at %d Hz from GPS %.1f (%.0f s)",
                len(raw), fs_in, t0, file_len)

    # NaNs are how GWOSC marks CAT1 failures (SoftwareX §3.2). check_quality
    # already vetoed the analysis span; zero any elsewhere so decimation's IIR
    # filter cannot propagate them across the whole file.
    n_nan = int(np.count_nonzero(~np.isfinite(raw)))
    if n_nan:
        logger.warning("zeroing %d non-finite samples outside the analysis span", n_nan)
        raw = np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)

    # --- decimate (SPEC.md D4) ---------------------------------------------
    if fs_in % args.sample_rate:
        raise RuntimeError(f"{fs_in} Hz is not an integer multiple of {args.sample_rate} Hz")
    q = fs_in // args.sample_rate
    fs = args.sample_rate
    if q > 1:
        # scipy's order-8 Chebyshev type I anti-alias filter -- the same
        # decimation family the SoftwareX paper documents for GWOSC's own
        # 16 kHz -> 4 kHz downsampling (§3, footnote 14).
        data = signal.decimate(raw, q, ftype="iir", zero_phase=True)
    else:
        data = raw.astype(np.float64)
    logger.info("Decimated %d -> %d Hz (q=%d), %d samples", fs_in, fs, q, len(data))

    # --- high-pass (SoftwareX §2.2) ----------------------------------------
    # Applied AFTER decimation so the filter is designed at the analysis rate.
    # Without it, sub-band power leaks across the segment FFT and inflates every
    # SNR (see gwlib.highpass for the failure this prevents).
    if args.highpass and args.highpass > 0:
        if args.highpass >= args.f_low:
            raise RuntimeError(
                f"--highpass {args.highpass} must sit below --f-low {args.f_low}")
        data = gwlib.highpass(data, fs, args.highpass)
        logger.info("High-passed at %.1f Hz (below the %.0f Hz analysis band)",
                    args.highpass, args.f_low)

    psd_start, off_starts = allocate_spans(
        clean, a_span, event_gps, args.psd_length, args.offsource_length,
        args.n_offsource, args.exclude,
    )
    logger.info("analysis @ %.1f (%.0fs, event %+.0fs from centre) | psd @ %.1f "
                "| %d off-source spans",
                a_start, analysis_len,
                event_gps - (a_start + analysis_len / 2.0), psd_start,
                len(off_starts))
    if len(off_starts) < args.n_offsource:
        # CLAUDE.md §2: report shortfalls, never silently reduce coverage.
        logger.warning("placed only %d of %d off-source spans (SPEC.md §9.4)",
                       len(off_starts), args.n_offsource)

    def cut(start, length):
        i0 = int(round((start - t0) * fs))
        return data[i0 : i0 + int(round(length * fs))]

    # --- PSD from off-source data (SPEC.md M4) -----------------------------
    psd_data = cut(psd_start, args.psd_length)
    f_psd, psd_raw, n_blocks = gwlib.estimate_psd(psd_data, fs)
    logger.info("PSD from %.0f s off-source: %d blocks of %.0f s",
                args.psd_length, n_blocks, gwlib.PSD_BLOCK_S)

    analysis = cut(a_start, analysis_len)
    n_seg = len(analysis)
    f_seg = gwlib.freqs(n_seg, 1.0 / fs)
    psd_interp = gwlib.interpolate_psd(f_psd, psd_raw, f_seg)
    psd_trunc = gwlib.truncate_inverse_psd(psd_interp, fs, n_seg)

    # Inverse-PSD truncation must preserve the in-band level. Getting it wrong
    # raises nothing -- it silently rescales every SNR by 1/sqrt(ratio) -- so
    # assert it here rather than trusting the implementation.
    trunc_ok, trunc_ratio = gwlib.check_truncation(
        psd_interp, psd_trunc, f_seg, args.f_low, args.f_high)
    logger.info("Inverse-PSD truncation preserved the in-band level to %.4f",
                trunc_ratio)
    if not trunc_ok:
        logger.error("inverse-PSD truncation changed the in-band level by "
                     "%.4gx; every SNR downstream would be wrong by "
                     "%.4gx. Refusing to continue.",
                     trunc_ratio, 1.0 / max(trunc_ratio, 1e-300) ** 0.5)
        sys.exit(1)

    # --- gate (Usman §3.2) -------------------------------------------------
    def condition(seg):
        f_s = gwlib.freqs(len(seg), 1.0 / fs)
        p_s = gwlib.truncate_inverse_psd(
            gwlib.interpolate_psd(f_psd, psd_raw, f_s), fs, len(seg))
        white = gwlib.whiten(seg, p_s, fs)
        gates = gwlib.find_gate_times(white, fs, args.gate_threshold)
        return gwlib.apply_gates(seg, fs, gates), gates

    analysis_gated, analysis_gates = condition(analysis)
    logger.info("Gated %d transient(s) in the analysis segment", len(analysis_gates))

    off_gated, off_gate_counts = [], []
    for s in off_starts:
        seg, gates = condition(cut(s, args.offsource_length))
        off_gated.append(seg)
        off_gate_counts.append(len(gates))

    # --- write -------------------------------------------------------------
    with h5py.File(args.output, "w") as f:
        g = f.create_group("analysis")
        g.create_dataset("strain", data=analysis_gated, compression="gzip")
        g.attrs["gps_start"] = a_start
        g.attrs["sample_rate"] = fs
        g.attrs["event_gps"] = event_gps
        g.attrs["edge_guard_s"] = EDGE_GUARD_S

        g = f.create_group("psd")
        g.create_dataset("freq", data=f_seg)
        g.create_dataset("psd", data=psd_interp)
        g.create_dataset("psd_truncated", data=psd_trunc)
        g.attrs["n_blocks"] = n_blocks
        g.attrs["block_s"] = gwlib.PSD_BLOCK_S
        g.attrs["psd_gps_start"] = psd_start
        g.attrs["psd_length_s"] = args.psd_length

        g = f.create_group("offsource")
        for i, (s, seg) in enumerate(zip(off_starts, off_gated)):
            d = g.create_dataset(f"span_{i:02d}", data=seg, compression="gzip")
            d.attrs["gps_start"] = s
        g.attrs["n_spans"] = len(off_gated)
        g.attrs["span_length_s"] = args.offsource_length

        g = f.create_group("meta")
        for k, v in (("event", meta["event"]), ("detector", meta["detector"])):
            g.attrs[k] = v
        g.attrs["f_low"] = args.f_low
        g.attrs["f_high"] = args.f_high
        g.attrs["sample_rate"] = fs
        g.attrs["source_sample_rate"] = fs_in
        g.attrs["gate_threshold"] = args.gate_threshold

    report = {
        "event": meta["event"], "detector": meta["detector"],
        "event_gps": event_gps,
        "sample_rate_in": fs_in, "sample_rate_out": fs, "decimation_q": q,
        "analysis": {"gps_start": a_start, "length_s": analysis_len,
                     "requested_length_s": args.analysis_length,
                     "shifted_off_centre": bool(quality.get("analysis_shifted")),
                     "truncated": bool(quality.get("analysis_truncated")),
                     "samples": int(n_seg), "gates_applied": len(analysis_gates),
                     "gate_sample_indices": analysis_gates},
        "psd": {"gps_start": psd_start, "length_s": args.psd_length,
                "n_blocks": int(n_blocks), "block_s": gwlib.PSD_BLOCK_S,
                "inverse_truncation_s": gwlib.INV_PSD_TRUNC_S,
                "truncation_inband_ratio": round(trunc_ratio, 5)},
        "offsource": {"requested": args.n_offsource, "placed": len(off_starts),
                      "length_s": args.offsource_length,
                      "exclusion_s": args.exclude,
                      "gps_starts": off_starts,
                      "gates_per_span": off_gate_counts},
        "band_hz": [args.f_low, args.f_high],
        "highpass_hz": args.highpass,
        "gate": {"threshold_whitened": args.gate_threshold,
                 "tukey_alpha": gwlib.GATE_TAPER_ALPHA,
                 "zero_s": gwlib.GATE_ZERO_S,
                 "provenance": "Usman 2016 §3.2 (window); threshold is a tuning "
                               "choice recorded in SPEC.md §9.3"},
        "nonfinite_samples_zeroed": n_nan,
        "science_mode_fraction": quality["checks"].get("science_mode_fraction"),
        "clean_runs": quality["checks"].get("clean_runs"),
    }
    with open(args.report, "w") as fh:
        json.dump(report, fh, indent=2)

    for path in (args.output, args.report):
        if not os.path.exists(path):
            logger.error("expected output missing: %s", path)
            sys.exit(1)
    logger.info("Output: %s, %s", args.output, args.report)


if __name__ == "__main__":
    main()
