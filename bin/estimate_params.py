#!/usr/bin/env python3
"""Point estimates from the loudest coincident candidate for one event.

SPEC.md M10, N1-N8. Reports exactly three things plus provenance:
  * network SNR
  * DETECTOR-FRAME chirp mass
  * merger time per detector, and the inter-detector delay

Deliberately absent, and each for a stated reason:
  * false-alarm rate / significance   -- SPEC.md N1, not computable here
  * credible intervals                -- SPEC.md N3, this is a point estimate
  * component masses m1, m2 separately -- SPEC.md N4, needs full PE
  * spins, final mass/spin, distance   -- SPEC.md N5-N7

SPEC.md M10: chirp mass is DETECTOR frame. Converting to source frame needs a
redshift, which needs a distance, which needs parameter estimation this workflow
does not do -- so a source-frame value appears only as a clearly-labelled
cross-check using the PUBLISHED redshift.
"""

import argparse
import json
import logging
import os
import sys

sys.path.insert(0, os.getcwd())
import gwlib

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser(description="Point estimates for one event")
    p.add_argument("--coincidence", required=True, help="JSON from coincidence.py")
    p.add_argument("--reference", required=True, help="gwtc1_reference.json")
    p.add_argument("--event", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    with open(args.coincidence) as fh:
        coinc = json.load(fh)
    with open(args.reference) as fh:
        ref = json.load(fh)
    ref_ev = {e["name"]: e for e in ref["events"]}.get(args.event)

    out = {
        "event": args.event,
        "detectors": coinc.get("detectors", []),
        "waveform_fallback_used": bool(coinc.get("waveform_fallback_used")),
        "frame": "detector",
        "recovered": None,
        "not_reported": {
            "false_alarm_rate": "SPEC.md N1 -- not computable from one segment",
            "significance_sigma": "SPEC.md N1",
            "credible_intervals": "SPEC.md N3 -- point estimate, not a posterior",
            "component_masses": "SPEC.md N4 -- requires full parameter estimation",
            "spins_distance_sky": "SPEC.md N5-N7",
        },
    }

    cands = coinc.get("candidates") or []
    if not cands:
        out["status"] = "no_candidate"
        out["note"] = coinc.get("note", "no coincident candidate above threshold")
        with open(args.output, "w") as fh:
            json.dump(out, fh, indent=2)
        logger.warning("%s: no candidate -> %s", args.event, args.output)
        return

    best = cands[0]
    dets = best["detectors"]
    times = best["gps_time"]
    # Earliest arrival is the reference; the published UTC in GWTC-1 Table I is
    # quoted to 0.1 s, so this is compared at that precision (SPEC.md V3a).
    t_ref = min(times.values())
    first = min(times, key=lambda k: times[k])

    rec = {
        "network_snr": round(best["network_snr"], 4),
        "single_detector_snr": {k: round(v, 4) for k, v in best["snr"].items()},
        "chirp_mass_detector_msun": round(best["chirp_mass_detector"], 5),
        "mass_ratio_q": round(best["q"], 4),
        "template_index": best["template_index"],
        "template_masses_msun": [round(best["mass1"], 4), round(best["mass2"], 4)],
        "merger_gps_time": {k: round(v, 6) for k, v in times.items()},
        "merger_gps_reference": round(t_ref, 6),
        "first_arrival_detector": first,
        "inter_detector_delay_ms": round(best["dt_seconds"] * 1e3, 4),
        "template_total_mass_detector_msun": round(best["mass1"] + best["mass2"], 4),
        "template_f_isco_hz": round(gwlib.f_isco(best["mass1"] + best["mass2"]), 2),
    }

    # Source-frame cross-check -- published redshift, NOT measured (SPEC.md N8).
    if ref_ev:
        z = ref_ev["source_parameters"]["redshift"]["median"]
        rec["chirp_mass_source_msun_using_published_z"] = round(
            best["chirp_mass_detector"] / (1.0 + z), 5)
        rec["published_redshift_used"] = z
        rec["source_frame_caveat"] = (
            "Derived with the PUBLISHED redshift, not a measured one "
            "(SPEC.md N8). Not an independent measurement.")

    out["recovered"] = rec
    out["status"] = "ok"
    out["n_candidates"] = coinc.get("n_candidates", 0)

    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    logger.info("%s: network SNR %.2f, Mc_det %.3f Msun, dt %.2f ms -> %s",
                args.event, rec["network_snr"], rec["chirp_mass_detector_msun"],
                rec["inter_detector_delay_ms"], args.output)


if __name__ == "__main__":
    main()
