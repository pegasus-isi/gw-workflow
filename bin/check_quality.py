#!/usr/bin/env python3
"""Verify data quality for the analysis span of one event/detector.

SPEC.md D6-D7, V5c. Reads the data-quality and injection bitmasks directly from
the GWOSC HDF5 file rather than a second API call -- the SoftwareX paper §3.2
documents both masks as 1 Hz bitmasks carried inside the same file.

Rather than pass/fail a fixed event-centred window, this reports the contiguous
CLEAN RUNS -- seconds passing DATA, CBC_CAT1 and CBC_CAT2 -- and identifies the
run containing the event. condition_strain then places the analysis segment
inside that run.

That distinction is load-bearing. GW151226's H1 detector comes online 2025 s into
its 4096 s file: only 53% is science mode, in a single clean run of 2171 s, with
the event just 123 s from its start. A fixed +/-128 s centred window misses by
five seconds even though 2048 s of clean data follows the event. Real pipelines
analyse contiguous science segments, not fixed windows around a time.

Checks:
  * a clean run contains the event and is long enough to be usable
  * no CBC hardware injection overlaps the usable window
  * no NaNs in the strain over that window (CAT1 failures are NaN by
    construction, SoftwareX §3.2)

SPEC.md D7 is explicit that passing all flags does NOT prove the span is free of
transient artifacts, so the verdict records the bits rather than declaring the
data clean.
"""

import argparse
import json
import logging
import os
import sys

import h5py
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# SoftwareX Table 3 / Table 4: bit names as published.
REQUIRED_DQ = ("DATA", "CBC_CAT1", "CBC_CAT2")
REQUIRED_INJ = ("NO_CBC_HW_INJ",)


def _names(grp, key):
    raw = grp[key][:]
    return [n.decode() if isinstance(n, (bytes, bytearray)) else str(n) for n in raw]


def clean_runs(f, gps_start, required=REQUIRED_DQ):
    """Contiguous [start, end) GPS intervals where every required bit is set."""
    grp = f["quality/simple"]
    mask = grp["DQmask"][:]
    names = _names(grp, "DQShortnames")
    good = np.ones(len(mask), dtype=bool)
    for name in required:
        if name not in names:
            raise RuntimeError(f"{name} absent from DQShortnames {names}")
        good &= ((mask >> names.index(name)) & 1).astype(bool)
    runs, i = [], 0
    while i < len(good):
        if good[i]:
            j = i
            while j < len(good) and good[j]:
                j += 1
            runs.append((gps_start + i, gps_start + j))
            i = j
        else:
            i += 1
    return runs, float(good.mean())


def place_window(run_start, run_end, event_gps, want, edge_guard, margin):
    """Largest usable window <= `want` inside [run_start, run_end).

    Prefers an event-centred window; shifts it when the event sits near a run
    boundary. Returns (start, length) or None when the run cannot host a usable
    window at all.
    """
    avail = run_end - run_start
    need_min = 2 * (edge_guard + margin)
    if avail < need_min:
        return None
    length = min(want, avail)
    start = event_gps - length / 2.0
    start = max(run_start, min(start, run_end - length))   # shift, don't shrink
    # The event must clear the filter-corrupted edges by `margin`.
    if not (start + edge_guard + margin <= event_gps <= start + length - edge_guard - margin):
        return None
    return start, length


def bit_report(f, group, mask_key, names_key, gps_start, lo, hi, required):
    """Fraction of seconds in [lo, hi) with each named bit set."""
    grp = f[group]
    mask = grp[mask_key][:]
    names = _names(grp, names_key)
    i0 = max(0, int(np.floor(lo - gps_start)))
    i1 = min(len(mask), int(np.ceil(hi - gps_start)))
    if i1 <= i0:
        raise RuntimeError(f"analysis span outside file coverage ({group})")
    window = mask[i0:i1]
    out, failures = {}, []
    for name in names:
        bit = names.index(name)
        frac = float(np.mean((window >> bit) & 1))
        out[name] = frac
        if name in required and frac < 1.0:
            failures.append(f"{name} set for only {frac*100:.3f}% of the span")
    return out, failures, (i0, i1)


def main():
    p = argparse.ArgumentParser(description="Data-quality gate for one event/detector")
    p.add_argument("--strain", required=True, help="GWOSC HDF5 file")
    p.add_argument("--meta", required=True, help="Provenance JSON from fetch_strain")
    p.add_argument("--analysis-length", type=float, default=256.0,
                   help="Requested analysis span length in seconds (default: 256)")
    p.add_argument("--edge-guard", type=float, default=16.0,
                   help="Filter-corrupted samples discarded per end (SPEC.md M5)")
    p.add_argument("--event-margin", type=float, default=4.0,
                   help="Extra seconds the event must clear the edge guard by")
    p.add_argument("--output", required=True, help="Output verdict JSON")
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    with open(args.meta) as fh:
        meta = json.load(fh)
    event_gps = float(meta["event_gps"])

    verdict = {
        "event": meta["event"], "detector": meta["detector"],
        "event_gps": event_gps,
        "requested_analysis_length_s": args.analysis_length,
        "checks": {}, "failures": [], "pass": False,
        "caveat": ("SPEC.md D7: passing every data-quality flag does not imply "
                   "the span is free of transient artifacts."),
    }

    with h5py.File(args.strain, "r") as f:
        gps_start = float(f["meta/GPSstart"][()])
        duration = float(f["meta/Duration"][()])
        fs = int(round(len(f["strain/Strain"]) / duration))
        verdict["file_gps_start"] = gps_start
        verdict["file_sample_rate"] = fs
        det = f["meta/Detector"][()]
        verdict["channel"] = det.decode() if isinstance(det, bytes) else str(det)

        runs, frac = clean_runs(f, gps_start)
        verdict["checks"]["science_mode_fraction"] = round(frac, 5)
        verdict["checks"]["n_clean_runs"] = len(runs)
        verdict["checks"]["clean_runs"] = [
            {"start": a, "end": b, "length_s": b - a} for a, b in runs]

        host = next((r for r in runs if r[0] <= event_gps < r[1]), None)
        if host is None:
            verdict["failures"].append(
                f"no clean run contains the event GPS {event_gps:.1f}; the "
                f"detector was not in science mode at that time")
            placed = None
        else:
            verdict["checks"]["event_clean_run"] = {
                "start": host[0], "end": host[1], "length_s": host[1] - host[0],
                "event_offset_from_start_s": round(event_gps - host[0], 3),
                "event_offset_from_end_s": round(host[1] - event_gps, 3),
            }
            placed = place_window(host[0], host[1], event_gps,
                                  args.analysis_length, args.edge_guard,
                                  args.event_margin)
            if placed is None:
                verdict["failures"].append(
                    f"clean run [{host[0]:.0f}, {host[1]:.0f}) ({host[1]-host[0]:.0f}s) "
                    f"cannot host a usable window: the event needs "
                    f"{args.edge_guard + args.event_margin:.0f}s of clearance at "
                    f"each end after the filter-corrupted edges")

        if placed is not None:
            a_start, a_len = placed
            lo, hi = a_start, a_start + a_len
            verdict["analysis_span"] = [lo, hi]
            verdict["analysis_length_s"] = a_len
            verdict["analysis_shifted"] = abs((lo + a_len / 2.0) - event_gps) > 1e-6
            verdict["analysis_truncated"] = a_len < args.analysis_length - 1e-6
            if verdict["analysis_truncated"]:
                logger.warning("analysis span shortened to %.0fs (requested %.0fs) "
                               "to fit the clean run", a_len, args.analysis_length)
            if verdict["analysis_shifted"]:
                logger.info("analysis span shifted off-centre: event sits %.0fs "
                            "into a %.0fs window", event_gps - lo, a_len)

            try:
                inj, inj_fail, _ = bit_report(f, "quality/injections", "Injmask",
                                              "InjShortnames", gps_start, lo, hi,
                                              REQUIRED_INJ)
                verdict["checks"]["injection_bits"] = inj
                verdict["failures"] += inj_fail
            except KeyError:
                # SoftwareX §3.2: Virgo performed no hardware injections in O2.
                verdict["checks"]["injection_bits"] = "absent from file"

            s0 = int(round((lo - gps_start) * fs))
            s1 = int(round((hi - gps_start) * fs))
            seg = f["strain/Strain"][max(0, s0):s1]
            n_nan = int(np.count_nonzero(~np.isfinite(seg)))
            verdict["checks"]["nan_samples_in_span"] = n_nan
            verdict["checks"]["span_samples"] = int(len(seg))
            if n_nan:
                verdict["failures"].append(
                    f"{n_nan} non-finite samples inside a span the DQ bits call "
                    f"clean -- the mask and the data disagree")

    verdict["pass"] = not verdict["failures"]
    with open(args.output, "w") as fh:
        json.dump(verdict, fh, indent=2)

    for msg in verdict["failures"]:
        logger.error("DQ failure: %s", msg)
    logger.info("Verdict: %s -> %s", args.output, "PASS" if verdict["pass"] else "FAIL")
    # A DQ failure vetoes this detector but must not abort the DAG: the verdict
    # is recorded and downstream jobs read it (SPEC.md §7.1).
    sys.exit(0)


if __name__ == "__main__":
    main()
