#!/usr/bin/env python3
"""Validate one event's recovered values against GWTC-1.

Implements SPEC.md §6.1-6.3. Every criterion emits recovered value, reference
value, tolerance and signed margin -- never a bare boolean.

Three things this file exists to get right:

1. FRAME CONVERSION (SPEC.md §6.2). A matched filter recovers a DETECTOR-frame
   chirp mass; GWTC-1 Table III publishes SOURCE frame. (1+z) runs 1.07-1.49
   across this event set, so conflating them is up to a 49% error at GW170729 --
   it would fail every event while looking like a physics bug. The published
   interval is propagated INTO the detector frame here, never the other way.

2. ASYMMETRIC SNR TOLERANCE (SPEC.md §6.1). GWTC-1 Table I's caption states its
   network SNR is "that of the template ranked highest by that search, which is
   not necessarily the template with the highest SNR" -- ranked under a
   reweighted statistic this workflow does not implement. Our maximum-SNR value
   can legitimately come out HIGHER, so exceeding the reference is not a failure.

3. TIMING WITHOUT LEAP SECONDS. GWTC-1 quotes UTC to 0.1 s. Rather than
   implementing GPS-UTC leap-second arithmetic, the recovered GPS merger time is
   compared against the event GPS the GWOSC event API returned -- the same
   instant the published UTC names, so the comparison is exact and no conversion
   can go wrong.
"""

import argparse
import json
import logging
import os
import sys

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

TIER_TOL = {
    # SPEC.md §6.1-6.3
    "A": {"snr_low": -0.15, "snr_high": 0.25, "mc_frac": 0.10, "time_ms": 100.0},
    "B": {"snr_low": -0.30, "snr_high": 0.40, "mc_frac": 0.20, "time_ms": 200.0},
}
SNR_FLOOR = 8.0          # SPEC.md V1a
SINGLE_SNR_FLOOR = 5.5   # SPEC.md V1c
LIGHT_TRAVEL_MS = 10.0   # SPEC.md V3b
# SPEC.md §6.1: GW170729 is exempt from the SNR band -- chi_eff = 0.37 against a
# non-spinning bank (M3) plus an ISCO of 35 Hz against f_low = 20 Hz (D5) make
# the deficit expected and physical.
SNR_BAND_EXEMPT = {"GW170729"}


def crit(name, ok, recovered, reference, tolerance, required=True, note=None):
    rec = {"criterion": name, "pass": bool(ok), "required": required,
           "recovered": recovered, "reference": reference, "tolerance": tolerance}
    if note:
        rec["note"] = note
    return rec


def main():
    p = argparse.ArgumentParser(description="Validate one event against GWTC-1")
    p.add_argument("--estimates", required=True, help="JSON from estimate_params")
    p.add_argument("--coincidence", required=True, help="JSON from coincidence")
    p.add_argument("--bank-manifest", required=True, help="JSON from build_bank")
    p.add_argument("--reference", required=True, help="gwtc1_reference.json")
    p.add_argument("--event", required=True)
    p.add_argument("--tier", choices=["A", "B"], required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    d = os.path.dirname(args.output)
    if d:
        os.makedirs(d, exist_ok=True)

    tol = TIER_TOL[args.tier]
    with open(args.estimates) as fh:
        est = json.load(fh)
    with open(args.coincidence) as fh:
        coinc = json.load(fh)
    with open(args.bank_manifest) as fh:
        manifest = json.load(fh)
    with open(args.reference) as fh:
        ref = json.load(fh)

    ev = {e["name"]: e for e in ref["events"]}.get(args.event)
    if ev is None:
        logger.error("%s not present in %s", args.event, args.reference)
        sys.exit(1)
    sp = ev["source_parameters"]
    pub_snr = (ev["searches"]["pycbc"] or {}).get("network_snr")

    out = {
        "event": args.event, "tier": args.tier,
        "tolerances": tol,
        "reference_source": "GWTC-1 (Abbott 2019, PRX 9, 031040), Tables I and III",
        "criteria": [], "failures": [], "verdict": None,
        "waveform_fallback_used": bool(est.get("waveform_fallback_used")),
    }

    rec = est.get("recovered")
    if rec is None:
        out["criteria"].append(crit(
            "V1a recovered as a coincidence", False, None,
            f"network SNR >= {SNR_FLOOR}", f">= {SNR_FLOOR}",
            required=(args.tier == "A"),
            note=coinc.get("note") or est.get("note")))
        out["failures"].append("V1a: no coincident candidate")
        out["verdict"] = "FAIL" if args.tier == "A" else "EXPECTED_MISS"
        if args.event == "GW170818":
            out["verdict"] = "EXPECTED_MISS"
            out["note"] = ("GW170818 was found by GstLAL alone (Hanford SNR 4.1); "
                           "PyCBC did not report it. A two-detector search "
                           "missing it reproduces the published result "
                           "(SPEC.md §1.2).")
        with open(args.output, "w") as fh:
            json.dump(out, fh, indent=2)
        logger.warning("%s -> %s", args.event, out["verdict"])
        return

    # ---------------- V1: network SNR ----------------------------------
    snr = rec["network_snr"]
    out["criteria"].append(crit(
        "V1a coincidence with network SNR >= 8", snr >= SNR_FLOOR,
        snr, f">= {SNR_FLOOR}", f">= {SNR_FLOOR}", required=True))
    if snr < SNR_FLOOR:
        out["failures"].append(f"V1a: network SNR {snr:.2f} < {SNR_FLOOR}")

    if pub_snr:
        lo = pub_snr * (1 + tol["snr_low"])
        hi = pub_snr * (1 + tol["snr_high"])
        ok = lo <= snr <= hi
        exempt = args.event in SNR_BAND_EXEMPT
        out["criteria"].append(crit(
            "V1b network SNR within band", ok or exempt, snr, pub_snr,
            f"[{lo:.2f}, {hi:.2f}] ({tol['snr_low']:+.0%}/{tol['snr_high']:+.0%})",
            required=not exempt,
            note=("exempt: chi_eff=0.37 against a non-spinning bank and ISCO "
                  "35 Hz vs f_low 20 Hz make the deficit physical "
                  "(SPEC.md §6.1)") if exempt else
                 (f"margin {(snr - pub_snr) / pub_snr:+.1%}; exceeding the "
                  f"reference is not a failure -- GWTC-1's value is the "
                  f"highest-RANKED template, not the highest-SNR one")))
        if not ok and not exempt:
            out["failures"].append(
                f"V1b: network SNR {snr:.2f} outside [{lo:.2f}, {hi:.2f}]")

    single_ok = any(v >= SINGLE_SNR_FLOOR for v in rec["single_detector_snr"].values())
    out["criteria"].append(crit(
        "V1c single-detector SNR >= 5.5 somewhere", single_ok,
        rec["single_detector_snr"], f">= {SINGLE_SNR_FLOOR}",
        f">= {SINGLE_SNR_FLOOR}", required=(args.tier == "A")))
    if not single_ok and args.tier == "A":
        out["failures"].append("V1c: no detector reached SNR 5.5")

    # ---------------- V2: chirp mass, DETECTOR frame -------------------
    mc = sp["chirp_mass_source_msun"]
    z = sp["redshift"]
    mc_det_med = mc["median"] * (1.0 + z["median"])
    # Conservative interval propagation into the detector frame.
    mc_det_hi = (mc["median"] + mc["plus"]) * (1.0 + z["median"] + z["plus"])
    mc_det_lo = (mc["median"] - mc["minus"]) * (1.0 + z["median"] - z["minus"])
    got = rec["chirp_mass_detector_msun"]

    inside = mc_det_lo <= got <= mc_det_hi
    frac = abs(got - mc_det_med) / mc_det_med
    within = frac <= tol["mc_frac"]
    out["criteria"].append(crit(
        "V2a/V2b chirp mass (detector frame)", inside or within, got,
        {"detector_frame_median": round(mc_det_med, 4),
         "detector_frame_90pct": [round(mc_det_lo, 4), round(mc_det_hi, 4)],
         "source_frame_published": mc, "redshift_published": z},
        f"inside the propagated interval, or within {tol['mc_frac']:.0%} of the median",
        required=True,
        note=(f"inside 90% interval: {inside}; fractional offset from median "
              f"{frac:+.2%}. Frame conversion applied (SPEC.md §6.2): "
              f"(1+z) = {1 + z['median']:.3f}")))
    if not (inside or within):
        out["failures"].append(
            f"V2: Mc_det {got:.3f} outside [{mc_det_lo:.3f}, {mc_det_hi:.3f}] "
            f"and {frac:.1%} from the median {mc_det_med:.3f}")

    lo_edge, hi_edge = manifest["chirp_mass_detector_range"]
    at_edge = abs(got - lo_edge) / lo_edge < 1e-3 or abs(got - hi_edge) / hi_edge < 1e-3
    out["criteria"].append(crit(
        "V2c recovered value not at a bank edge", not at_edge, got,
        {"bank_mc_range": [lo_edge, hi_edge]}, "not within 0.1% of either edge",
        required=True,
        note="catches a rail-pinned fit masquerading as a measurement"))
    if at_edge:
        out["failures"].append(f"V2c: Mc_det {got:.3f} sits on a bank edge")

    spacing = next((e.get("local_spacing_fractional")
                    for e in manifest["adequacy_gate"]["events"]
                    if e["event"] == args.event), None)
    sp_ok = spacing is not None and spacing < tol["mc_frac"]
    out["criteria"].append(crit(
        "V2d bank spacing finer than the claimed tolerance", sp_ok, spacing,
        f"< {tol['mc_frac']}", f"< {tol['mc_frac']}", required=True,
        note="prevents claiming an agreement the bank's own resolution cannot support"))
    if not sp_ok:
        out["failures"].append(
            f"V2d: local bank spacing {spacing} not finer than {tol['mc_frac']}")

    # ---------------- V3: timing --------------------------------------
    det_meta = coinc.get("meta") or {}
    event_gps = next((m["event_gps"] for m in det_meta.values() if "event_gps" in m), None)
    if event_gps is not None:
        off_ms = abs(rec["merger_gps_reference"] - event_gps) * 1e3
        ok = off_ms <= tol["time_ms"]
        out["criteria"].append(crit(
            "V3a merger time vs published event time", ok, round(off_ms, 2),
            {"event_gps_from_gwosc_api": event_gps,
             "published_utc": f"{ev['date_utc']} {ev['time_utc']}"},
            f"<= {tol['time_ms']} ms", required=True,
            note="compared in GPS against the API value for the published UTC, "
                 "so no leap-second conversion is involved"))
        if not ok:
            out["failures"].append(
                f"V3a: merger time {off_ms:.1f} ms from the published event time")

    dt_ms = abs(rec["inter_detector_delay_ms"])
    ok = dt_ms <= LIGHT_TRAVEL_MS
    out["criteria"].append(crit(
        "V3b |dt| <= inter-site light travel", ok, round(dt_ms, 3),
        f"<= {LIGHT_TRAVEL_MS} ms", f"<= {LIGHT_TRAVEL_MS} ms", required=True,
        note="internal consistency check; needs no published value"))
    if not ok:
        out["failures"].append(f"V3b: |dt| {dt_ms:.2f} ms exceeds light travel")

    if args.event == "GW150914":
        # PRL 116 061102: GW150914 arrived at L1 first, 6.9 +0.5/-0.4 ms before H1.
        # The only per-event inter-detector delay any of the four papers publishes.
        first = rec["first_arrival_detector"]
        ok = first == "L1" and 6.9 - 0.4 - 1.0 <= dt_ms <= 6.9 + 0.5 + 1.0
        out["criteria"].append(crit(
            "V3c GW150914 arrival order and delay", ok,
            {"first_arrival": first, "delay_ms": round(dt_ms, 3)},
            {"expected_first": "L1", "delay_ms": "6.9 +0.5 -0.4"},
            "L1 first, delay within the published interval (+/-1 ms allowance)",
            required=True,
            note="PRL 116 061102 -- the only published per-event delay"))
        if not ok:
            out["failures"].append(
                f"V3c: got {first} first with dt {dt_ms:.2f} ms; expected L1 "
                f"first at 6.9 +0.5 -0.4 ms")

    # ---------------- verdict ------------------------------------------
    required_failed = [c for c in out["criteria"] if c["required"] and not c["pass"]]
    if out["waveform_fallback_used"]:
        # SPEC.md M2 / CLAUDE.md §2: a toy waveform must never yield a PASS.
        out["verdict"] = "INVALID_FALLBACK_WAVEFORM"
        out["note"] = ("The restricted-PN fallback waveform was used, not the "
                       "requested approximant. Adequate for a smoke test only; "
                       "no validation verdict can be drawn (SPEC.md M2).")
    elif not required_failed:
        out["verdict"] = "PASS"
    elif args.tier == "B":
        out["verdict"] = "PASS_WITH_NOTES"
    else:
        out["verdict"] = "FAIL"

    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    for msg in out["failures"]:
        logger.error("%s %s", args.event, msg)
    logger.info("%s (tier %s) -> %s", args.event, args.tier, out["verdict"])


if __name__ == "__main__":
    main()
