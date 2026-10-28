#!/usr/bin/env python3
"""Off-source null test (SPEC.md V5a) — a HARD gate.

Compares the loudest coincident candidate found in each off-source span against
the network SNR recovered at the true event time. If noise alone produces
something as loud, the on-source "recovery" carries no weight, whatever it
matched in GWTC-1.

The off-source spans were filtered by matched_filter.py with --dataset, i.e. the
identical code path and configuration as the on-source segment, and then passed
through the same coincidence.py. This job only compares.

SPEC.md §6.6: a control failure is a HARD fail regardless of how well the event
comparisons went.
"""

import argparse
import json
import logging
import os
import sys

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def loudest(path):
    with open(path) as fh:
        blob = json.load(fh)
    cands = blob.get("candidates") or []
    return (max(c["network_snr"] for c in cands) if cands else 0.0), blob


def was_searched(blob):
    """A span with zero coincident candidates is only a quiet result if the
    single-detector filters actually produced triggers. Zero triggers in every
    detector means the span was never really searched (empty input, wrong
    dataset, a silent upstream failure), and counting it as quiet would let
    the null test pass vacuously."""
    if blob.get("n_candidates", 0) > 0:
        return True
    singles = blob.get("n_single_triggers") or {}
    return any(v > 0 for v in singles.values())


def main():
    p = argparse.ArgumentParser(description="Off-source null test (SPEC.md V5a)")
    p.add_argument("--offsource", action="append", required=True,
                   help="Coincidence JSON for one off-source span (repeatable)")
    p.add_argument("--onsource", required=True,
                   help="Coincidence JSON for the true event time")
    p.add_argument("--event", required=True, help="Event the spans came from")
    p.add_argument("--min-spans", type=int, default=10,
                   help="Minimum off-source spans for a valid null test "
                        "(SPEC.md V5a requires 10)")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    on_snr, on_blob = loudest(args.onsource)
    spans = []
    unsearched = []
    for path in args.offsource:
        snr, blob = loudest(path)
        singles = blob.get("n_single_triggers") or {}
        spans.append({"file": os.path.basename(path),
                      "loudest_network_snr": round(snr, 4),
                      "n_candidates": blob.get("n_candidates", 0),
                      "n_single_triggers": singles})
        if not was_searched(blob):
            unsearched.append(os.path.basename(path))
        logger.info("%s: loudest off-source network SNR %.2f "
                    "(single-detector triggers: %s)",
                    spans[-1]["file"], snr, singles or "none recorded")

    worst = max((s["loudest_network_snr"] for s in spans), default=0.0)
    ok = (worst < on_snr and on_snr > 0.0
          and len(spans) >= args.min_spans and not unsearched)

    out = {
        "control": "offsource_null",
        "criterion": "V5a",
        "event": args.event,
        "n_spans": len(spans),
        "min_spans": args.min_spans,
        "n_unsearched_spans": len(unsearched),
        "onsource_network_snr": round(on_snr, 4),
        "loudest_offsource_network_snr": round(worst, 4),
        "spans": spans,
        "pass": bool(ok),
        "hard_gate": True,
        "separation": round(on_snr - worst, 4),
        "note": ("SPEC.md §6.6: a control failure is a hard fail. Off-source "
                 "spans were filtered by the same matched_filter.py code path "
                 "with --dataset, then the same coincidence.py."),
    }
    if on_snr <= 0.0:
        out["failure_reason"] = ("no on-source candidate to compare against; the "
                                "null test is vacuous")
    elif len(spans) < args.min_spans:
        out["failure_reason"] = (
            f"only {len(spans)} off-source span(s) supplied, fewer than the "
            f"{args.min_spans} SPEC.md V5a requires: the null test is "
            f"under-sampled and cannot pass")
    elif unsearched:
        out["failure_reason"] = (
            f"{len(unsearched)} span(s) show zero coincident candidates AND "
            f"zero single-detector triggers ({', '.join(unsearched)}): those "
            f"spans were never actually searched, so their quiet result is "
            f"vacuous and the null test cannot pass")
    elif not ok:
        out["failure_reason"] = (
            f"an off-source span reached network SNR {worst:.2f} >= the on-source "
            f"{on_snr:.2f}: noise alone is as loud as the claimed recovery")

    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    logger.info("V5a null test: %s (on-source %.2f vs loudest off-source %.2f)",
                "PASS" if ok else "FAIL", on_snr, worst)


if __name__ == "__main__":
    main()
