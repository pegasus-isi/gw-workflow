#!/usr/bin/env python3
"""Combine per-detector triggers into coincident candidates for one event.

SPEC.md M8-M9. Usman §3.5: a candidate must appear in both detectors within
15 ms (10 ms inter-site light travel + 5 ms timing allowance) AND in the SAME
template -- exact-match coincidence, which is only possible because one bank is
shared across detectors (§3.1).

Usman §3.6: a final clustering step keeps only the loudest coincidence in a
10 s window, so one signal yields one statistically independent candidate rather
than a correlated cluster across the bank.

The ranking statistic is the plain quadrature sum of single-detector SNR
(Usman eq. 16). SPEC.md N15/M12: the GWTC-1 ranking statistics are NOT
reimplemented, and no reweighting is applied, so this value is directly
comparable to the network SNR published in GWTC-1 Table I.
"""

import argparse
import itertools
import json
import logging
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.getcwd())
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

FINAL_CLUSTER_S = 10.0     # Usman §3.6


def coincide(by_det, window_s, detectors):
    """Exact-template, time-coincident pairs across detector pairs."""
    cands = []
    for d1, d2 in itertools.combinations(detectors, 2):
        a = defaultdict(list)
        for t in by_det[d1]:
            a[t["template_index"]].append(t)
        for t2 in by_det[d2]:
            for t1 in a.get(t2["template_index"], ()):
                dt = t2["gps_time"] - t1["gps_time"]
                if abs(dt) <= window_s:
                    cands.append({
                        "detectors": [d1, d2],
                        "template_index": t1["template_index"],
                        "chirp_mass_detector": t1["chirp_mass_detector"],
                        "q": t1["q"], "mass1": t1["mass1"], "mass2": t1["mass2"],
                        "snr": {d1: t1["snr"], d2: t2["snr"]},
                        "gps_time": {d1: t1["gps_time"], d2: t2["gps_time"]},
                        "dt_seconds": dt,
                        # Usman eq. (16)
                        "network_snr": (t1["snr"] ** 2 + t2["snr"] ** 2) ** 0.5,
                    })
    return cands


def final_cluster(cands, window_s):
    """Usman §3.6: keep the loudest candidate per `window_s` window."""
    kept = []
    for c in sorted(cands, key=lambda c: -c["network_snr"]):
        t = min(c["gps_time"].values())
        if all(abs(t - min(k["gps_time"].values())) > window_s for k in kept):
            kept.append(c)
    return kept


def main():
    p = argparse.ArgumentParser(description="Coincidence for one event")
    p.add_argument("--triggers", action="append", required=True,
                   help="Per-chunk trigger JSON (repeatable)")
    p.add_argument("--event", required=True)
    p.add_argument("--coincidence-window-ms", type=float, default=gwlib.COINC_WINDOW_MS,
                   help="Usman §3.5: 15 ms = 10 ms light travel + 5 ms timing")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    by_det, meta, fallback, truncated = defaultdict(list), {}, False, False
    for path in args.triggers:
        with open(path) as fh:
            blob = json.load(fh)
        det = blob["detector"]
        by_det[det].extend(blob["triggers"])
        fallback |= bool(blob.get("waveform_fallback_used"))
        truncated |= bool(blob.get("truncated"))
        meta.setdefault(det, {
            "event_gps": blob["event_gps"],
            "segment_gps_start": blob["segment_gps_start"],
            "snr_threshold": blob["snr_threshold"],
            "band_hz": blob["band_hz"],
            "n_templates_total": blob["n_templates_total"],
        })

    detectors = sorted(by_det)
    for d_ in detectors:
        logger.info("%s: %d single-detector triggers", d_, len(by_det[d_]))

    window_s = args.coincidence_window_ms / 1000.0
    result = {
        "event": args.event,
        "detectors": detectors,
        "n_single_triggers": {d_: len(by_det[d_]) for d_ in detectors},
        "coincidence_window_ms": args.coincidence_window_ms,
        "final_cluster_window_s": FINAL_CLUSTER_S,
        "template_match": "exact (Usman §3.5)",
        "ranking_statistic": "quadrature sum of single-detector SNR (Usman eq. 16)",
        "waveform_fallback_used": fallback,
        "input_chunks_truncated": truncated,
        "meta": meta,
        "no_significance": ("SPEC.md N1: false-alarm rate and significance are "
                            "NOT computable from a single segment and are not "
                            "reported."),
    }

    if len(detectors) < 2:
        result["candidates"] = []
        result["n_candidates"] = 0
        result["note"] = (
            f"only {len(detectors)} detector(s) produced triggers; a two-detector "
            f"coincidence is impossible. For GW170818 this reproduces PyCBC's "
            f"published non-detection (SPEC.md §1.2)."
        )
        logger.warning("%s", result["note"])
    else:
        cands = final_cluster(coincide(by_det, window_s, detectors), FINAL_CLUSTER_S)
        cands.sort(key=lambda c: -c["network_snr"])
        result["candidates"] = cands[:100]
        result["n_candidates"] = len(cands)
        if cands:
            logger.info("Loudest: network SNR %.2f, Mc_det %.3f, dt %.2f ms",
                        cands[0]["network_snr"], cands[0]["chirp_mass_detector"],
                        cands[0]["dt_seconds"] * 1e3)
        else:
            logger.warning("no coincident candidates for %s", args.event)

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)
    logger.info("Output: %s (%d candidates)", args.output, result["n_candidates"])


if __name__ == "__main__":
    main()
