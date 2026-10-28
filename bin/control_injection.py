#!/usr/bin/env python3
"""Injection-recovery control (SPEC.md V5b) — a HARD gate.

Inject a waveform of KNOWN parameters into off-source data, search the full bank,
and check the recovered chirp mass and time.

This is the criterion that separates a working pipeline from a lucky one: the
true answer is known by construction, so there is no published value to
accidentally fit toward. If V5a (null test) says noise is quiet and V5b says a
known signal is recovered correctly, then an on-source recovery means something.

Runs on a shorter crop than the on-source segment (segment length does not affect
chirp-mass or time recovery) so the full bank can be searched per injection
inside the job budget. The crop length is reported.
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
from matched_filter import filter_segment

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser(description="Injection recovery control (SPEC.md V5b)")
    p.add_argument("--conditioned", required=True, help="Conditioned HDF5 (off-source source)")
    p.add_argument("--bank", required=True)
    p.add_argument("--span", default="offsource/span_00",
                   help="Off-source dataset to inject into")
    p.add_argument("--chirp-mass", type=float, required=True,
                   help="Injected DETECTOR-frame chirp mass (Msun)")
    p.add_argument("--q", type=float, default=1.2, help="Injected mass ratio")
    p.add_argument("--target-snr", type=float, default=15.0,
                   help="Scale the injection to this optimal SNR")
    p.add_argument("--crop-length", type=float, default=64.0,
                   help="Seconds of off-source data to inject into and filter")
    p.add_argument("--snr-threshold", type=float, default=5.5)
    p.add_argument("--mc-tolerance", type=float, default=0.10,
                   help="SPEC.md V5b: recover Mc within this fraction (V2b tolerance)")
    p.add_argument("--time-tolerance-ms", type=float, default=100.0,
                   help="SPEC.md V5b: recover time within this (V3a tolerance)")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    with h5py.File(args.conditioned, "r") as f:
        span = f[args.span]
        data = span[:]
        span_gps = float(span.attrs["gps_start"])
        fs = int(f["analysis/strain"].parent.attrs["sample_rate"])
        psd_full = f["psd/psd_truncated"][:]
        psd_freq = f["psd/freq"][:]
        meta = dict(f["meta"].attrs)
        f_low = float(meta.get("f_low", 20.0))
        f_high = float(meta.get("f_high", 512.0))

    with h5py.File(args.bank, "r") as f:
        templates = [tuple(t) for t in f["templates"][:]]
        approximant = f["/"].attrs.get("approximant", "IMRPhenomD")
        bank_fallback = bool(f["/"].attrs.get("waveform_fallback_used", False))

    # --- crop, well inside the filter-corrupted edges ---------------------
    n_crop = int(args.crop_length * fs)
    if n_crop > len(data):
        logger.error("crop %.0f s exceeds the %.0f s span", args.crop_length,
                     len(data) / fs)
        sys.exit(2)
    # SPEC.md M5 discards INV_PSD_TRUNC_S at EACH end, so a crop of exactly
    # 2 x that length leaves no usable samples and every trigger is silently
    # dropped -- which looks identical to "the injection was not recovered".
    min_crop = 2 * gwlib.INV_PSD_TRUNC_S
    if args.crop_length <= min_crop:
        logger.error("--crop-length %.0f s leaves no usable samples: %.0f s is "
                     "discarded at each end as filter-corrupted (SPEC.md M5), "
                     "so the crop must exceed %.0f s.",
                     args.crop_length, gwlib.INV_PSD_TRUNC_S, min_crop)
        sys.exit(2)
    off = (len(data) - n_crop) // 2
    seg = data[off : off + n_crop].copy()
    seg_gps = span_gps + off / fs

    n = len(seg)
    dt = 1.0 / fs
    df = 1.0 / (n * dt)
    f = gwlib.freqs(n, dt)
    psd = gwlib.interpolate_psd(psd_freq, psd_full, f)
    mask = gwlib.band_mask(f, f_low, f_high)

    # --- build the injection ---------------------------------------------
    m1, m2 = gwlib.mass1_mass2(args.chirp_mass, args.q)
    h, fb = gwlib.fd_waveform(m1, m2, f_low, df, f_high, approximant)
    if len(h) < len(f):
        h = np.pad(h, (0, len(f) - len(h)))
    h = h[: len(f)]
    sigma = np.sqrt(gwlib.sigma_sq(h, psd, df, mask))
    if sigma <= 0:
        logger.error("injection template has zero norm")
        sys.exit(1)
    h_scaled = h * (args.target_snr / sigma)

    # Place the merger at the crop centre. The time-domain waveform from an
    # f-domain template wraps at the segment edge, so roll it to the centre.
    h_td = np.fft.irfft(h_scaled, n=n) / dt
    h_td = np.roll(h_td, n // 2)
    inj_gps = seg_gps + (n // 2) / fs
    seg_inj = seg + h_td

    logger.info("Injected Mc_det=%.3f q=%.2f (m1=%.2f m2=%.2f) at optimal SNR "
                "%.1f, GPS %.4f", args.chirp_mass, args.q, m1, m2,
                args.target_snr, inj_gps)

    triggers, fb2 = filter_segment(
        seg_inj, fs, psd, f, templates, f_low, f_high, args.snr_threshold,
        approximant, gwlib.INV_PSD_TRUNC_S, seg_gps, progress_every=200,
    )

    out = {
        "control": "injection_recovery",
        "criterion": "V5b",
        "hard_gate": True,
        "event_context": meta.get("event"),
        "detector": meta.get("detector"),
        "span": args.span,
        "crop_length_s": args.crop_length,
        "crop_note": ("shorter than the on-source segment; segment length does "
                      "not affect chirp-mass or time recovery, and this keeps a "
                      "full-bank search inside the job budget"),
        "injected": {"chirp_mass_detector_msun": args.chirp_mass, "q": args.q,
                     "mass1": round(m1, 4), "mass2": round(m2, 4),
                     "optimal_snr": args.target_snr, "gps_time": round(inj_gps, 6)},
        "tolerances": {"chirp_mass_fractional": args.mc_tolerance,
                       "time_ms": args.time_tolerance_ms},
        "waveform_fallback_used": bool(fb or fb2 or bank_fallback),
        "n_triggers": len(triggers),
    }

    if not triggers:
        out["recovered"] = None
        out["pass"] = False
        out["failure_reason"] = ("no trigger above threshold; a known-parameter "
                                 "injection at optimal SNR "
                                 f"{args.target_snr} was not recovered")
    else:
        best = max(triggers, key=lambda t: t["snr"])
        mc_frac = abs(best["chirp_mass_detector"] - args.chirp_mass) / args.chirp_mass
        dt_ms = abs(best["gps_time"] - inj_gps) * 1e3
        mc_ok = mc_frac <= args.mc_tolerance
        t_ok = dt_ms <= args.time_tolerance_ms
        out["recovered"] = {
            "chirp_mass_detector_msun": round(best["chirp_mass_detector"], 5),
            "q": round(best["q"], 4),
            "snr": round(best["snr"], 4),
            "gps_time": round(best["gps_time"], 6),
            "chirp_mass_fractional_error": round(mc_frac, 5),
            "time_error_ms": round(dt_ms, 3),
            "chirp_mass_within_tolerance": bool(mc_ok),
            "time_within_tolerance": bool(t_ok),
        }
        out["pass"] = bool(mc_ok and t_ok)
        if not out["pass"]:
            bad = []
            if not mc_ok:
                bad.append(f"Mc error {mc_frac:.1%} > {args.mc_tolerance:.0%}")
            if not t_ok:
                bad.append(f"time error {dt_ms:.1f} ms > {args.time_tolerance_ms:.0f} ms")
            out["failure_reason"] = "; ".join(bad)

    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    logger.info("V5b Mc=%.2f: %s", args.chirp_mass,
                "PASS" if out.get("pass") else f"FAIL ({out.get('failure_reason')})")


if __name__ == "__main__":
    main()
