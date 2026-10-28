#!/usr/bin/env python3
"""Render a time-frequency spectrogram of whitened strain around one event.

SPEC.md V4. A constant-Q transform of the whitened analysis segment, which is
the representation the PRL 116 061102 Fig. 1 bottom row and GWTC-1 Fig. 10 use.

SPEC.md V4a/V4b are explicitly NON-BLOCKING: a spectrogram is a visualisation,
and GWTC-1's own Fig. 10 shows several events with no visible track. This job
therefore never fails the run on the absence of a chirp.
"""

import argparse
import json
import logging
import os
import sys

import h5py
import numpy as np

# Point matplotlib's cache somewhere writable BEFORE importing it. In the
# container $HOME is not writable, which otherwise emits a warning in every job
# log and slows the import.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".mplcache"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg")          # headless: no DISPLAY on a worker node
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def q_transform(x, fs, f_low, f_high, q=12.0, n_freq=180):
    """Constant-Q transform via frequency-domain Gaussian windows.

    Each row is the analytic-signal energy in a band of width f0/q centred on
    f0, so the time-frequency resolution tracks frequency the way a chirp needs.
    """
    n = len(x)
    xf = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    rows, freqs = [], np.geomspace(f_low, f_high, n_freq)
    for f0 in freqs:
        sigma = f0 / q
        win = np.exp(-0.5 * ((f - f0) / sigma) ** 2)
        # Analytic signal: one-sided spectrum, doubled, complex inverse FFT.
        full = np.zeros(n, dtype=np.complex128)
        full[: len(xf)] = xf * win * 2.0
        rows.append(np.abs(np.fft.ifft(full)) ** 2)
    return freqs, np.array(rows)


def main():
    p = argparse.ArgumentParser(description="Q-transform spectrogram for one event/detector")
    p.add_argument("--conditioned", required=True)
    p.add_argument("--half-width", type=float, default=1.0,
                   help="Seconds either side of the event to plot")
    p.add_argument("--q", type=float, default=12.0)
    p.add_argument("--f-low", type=float, default=20.0)
    p.add_argument("--f-high", type=float, default=512.0)
    p.add_argument("--output", required=True, help="Output PNG")
    p.add_argument("--report", required=True, help="Output JSON summary")
    args = p.parse_args()

    for path in (args.output, args.report):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    with h5py.File(args.conditioned, "r") as f:
        strain = f["analysis/strain"][:]
        attrs = dict(f["analysis/strain"].parent.attrs)
        psd = f["psd/psd_truncated"][:]
        meta = dict(f["meta"].attrs)

    fs = int(attrs["sample_rate"])
    gps_start = float(attrs["gps_start"])
    event_gps = float(attrs["event_gps"])
    event, det = meta.get("event"), meta.get("detector")

    white = gwlib.whiten(strain, psd, fs)

    # Crop around the event, staying inside the filter-corrupted edges.
    c = int((event_gps - gps_start) * fs)
    half = int(args.half_width * fs)
    guard = int(gwlib.INV_PSD_TRUNC_S * fs)
    lo = max(guard, c - half)
    hi = min(len(white) - guard, c + half)
    if hi - lo < fs // 4:
        logger.error("crop around the event is too short (%d samples)", hi - lo)
        sys.exit(1)
    seg = white[lo:hi]
    t = (np.arange(lo, hi) / fs) + gps_start - event_gps

    freqs, power = q_transform(seg, fs, args.f_low, args.f_high, args.q)
    # Normalise each row so a loud band does not swamp the colour scale.
    norm = power / np.median(power, axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    m = ax.pcolormesh(t, freqs, norm, shading="auto", cmap="viridis",
                      vmin=0, vmax=max(8.0, float(np.percentile(norm, 99.9))))
    ax.set_yscale("log")
    ax.set_ylim(args.f_low, args.f_high)
    ax.set_xlabel(f"Time from GPS {event_gps:.3f} (s)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title(f"{event} — {det} — whitened Q-transform (Q={args.q:g})")
    fig.colorbar(m, ax=ax, label="Normalised energy")
    fig.tight_layout()
    fig.savefig(args.output, dpi=130)
    plt.close(fig)

    peak = np.unravel_index(int(np.argmax(norm)), norm.shape)
    report = {
        "event": event, "detector": det, "png": os.path.basename(args.output),
        "q": args.q, "band_hz": [args.f_low, args.f_high],
        "half_width_s": args.half_width,
        "peak_normalised_energy": round(float(norm[peak]), 3),
        "peak_frequency_hz": round(float(freqs[peak[0]]), 2),
        "peak_time_from_event_s": round(float(t[peak[1]]), 4),
        "blocking": False,
        "note": ("SPEC.md V4a/V4b are non-blocking; absence of a visible chirp "
                 "is not a failure. GWTC-1 Fig. 10 shows several events with no "
                 "visible track."),
    }
    with open(args.report, "w") as fh:
        json.dump(report, fh, indent=2)
    logger.info("%s %s: peak %.1f at %.1f Hz -> %s", event, det,
                report["peak_normalised_energy"], report["peak_frequency_hz"],
                args.output)


if __name__ == "__main__":
    main()
