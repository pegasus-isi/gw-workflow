#!/usr/bin/env python3
"""Aggregate every per-event validation and control into one report.

SPEC.md §5.1, §6.6. Emits report.md (for humans) and report.json (for machines).

SPEC.md §5.3 requires the report ITSELF -- not just the spec -- to state that
this workflow produces no statement of statistical significance, and that
recovering a signal at a known time is not a detection in the sense the
published papers mean. That text is written into the Markdown below.

Verdict rules (SPEC.md §6.6):
  PASS             all V5 controls pass and every Tier A event passes its
                   required criteria
  PASS_WITH_NOTES  the above except for named, explained deviations, or any
                   Tier B outcome
  FAIL             any V5 control fails, or any Tier A event fails
                   V1a / V2c / V2d / V3b
  INCOMPLETE       controls were skipped (SPEC.md §7.3) -- never PASS
"""

import argparse
import json
import logging
import os
import sys
from datetime import datetime, timezone

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

HARD_FAIL_CRITERIA = ("V1a", "V2c", "V2d", "V3b")   # SPEC.md §6.6

NON_EXPECTATIONS = """\
## What this workflow does NOT establish

- **No statistical significance.** No false-alarm rate, no inverse FAR, no
  sigma value appears anywhere in this report. Those come from time-shifting an
  entire observing run (`T_b = T^2/delta`; roughly 15 days of coincident data to
  reach 5 sigma, and 608,000 equivalent years of background for GW150914). A
  4096 s segment cannot produce a background estimate, so the quantity is not
  merely unreported here -- it is not computable from this input.
- **Recovering a signal at a known time is not a detection** in the sense the
  published papers mean. This is a targeted reanalysis at times taken from the
  catalogue; it cannot discover anything.
- **No posteriors.** Every number below is a point estimate from the
  best-matching template, not a credible interval.
- **GWTC-1 supersedes the individual detection papers.** Its GW150914 values
  differ from the 2016 PRL (for instance d_L = 440 +150 -170 against
  410 +160 -180 Mpc) because of improved PSD estimation, spline calibration
  envelopes, and O2 noise subtraction.
"""


def load_all(paths):
    out = []
    for path in paths or []:
        try:
            with open(path) as fh:
                out.append((os.path.basename(path), json.load(fh)))
        except Exception as exc:
            logger.error("could not read %s: %s", path, exc)
            out.append((os.path.basename(path), {"_read_error": str(exc)}))
    return out


def main():
    p = argparse.ArgumentParser(description="Aggregate the workflow's results")
    p.add_argument("--validation", action="append", help="Per-event validation JSON")
    p.add_argument("--control", action="append", help="Control JSON (V5a / V5b)")
    p.add_argument("--bank-manifest", required=True)
    p.add_argument("--conditioning", action="append", help="Conditioning report JSON")
    p.add_argument("--spectrogram", action="append", help="Spectrogram report JSON")
    p.add_argument("--controls-skipped", action="store_true",
                   help="Set when --skip-controls was used; forces INCOMPLETE")
    p.add_argument("--min-injections", type=int, default=5,
                   help="Injection controls required for a valid run "
                        "(SPEC.md V5b requires 5, at distinct chirp masses)")
    p.add_argument("--expected-event", action="append", default=[],
                   help="Event that MUST have a validation result (repeatable; "
                        "the generator passes the full event list)")
    p.add_argument("--output-md", required=True)
    p.add_argument("--output-json", required=True)
    args = p.parse_args()

    for path in (args.output_md, args.output_json):
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    with open(args.bank_manifest) as fh:
        bank = json.load(fh)
    validations = load_all(args.validation)
    controls = load_all(args.control)
    conditioning = load_all(args.conditioning)
    spectrograms = load_all(args.spectrogram)

    # The fallback waveform invalidates controls just as it invalidates the
    # bank and the event estimates: an injection recovered with the wrong
    # waveform proves nothing.
    fallback = (bool(bank.get("waveform_fallback_used"))
                or any(v.get("waveform_fallback_used") for _, v in validations)
                or any(c.get("waveform_fallback_used") for _, c in controls))

    # ---- verdict -------------------------------------------------------
    control_failures, tier_a_hard, tier_a_soft, tier_b_notes = [], [], [], []
    for name, c in controls:
        if not c.get("pass"):
            control_failures.append(
                f"{c.get('criterion', '?')} ({name}): {c.get('failure_reason', 'failed')}")
            continue
        # Audit a claimed pass rather than trusting it. A control that says
        # pass:true while its own numbers contradict that is a failure, and
        # one that predates the audit fields earns a note, never a clean PASS.
        if c.get("control") == "offsource_null":
            if "n_unsearched_spans" in c:
                if c["n_unsearched_spans"] > 0:
                    control_failures.append(
                        f"V5a ({name}): claims pass but records "
                        f"{c['n_unsearched_spans']} unsearched span(s)")
                elif c.get("n_spans", 0) < c.get("min_spans", 10):
                    control_failures.append(
                        f"V5a ({name}): claims pass but records only "
                        f"{c.get('n_spans')} of {c.get('min_spans', 10)} spans")
            else:
                tier_b_notes.append(
                    f"offsource control {name} predates the span-audit fields, "
                    f"so span search adequacy is not machine-verified here")
        elif c.get("control") == "injection_recovery":
            if c.get("n_triggers") == 0:
                control_failures.append(
                    f"V5b ({name}): claims pass but recorded zero triggers, "
                    f"so nothing was actually recovered")

    # The bank's own adequacy gate (SPEC.md B6) is a required criterion, not
    # a decoration: a verdict computed over an inadequate bank is not a PASS.
    # The gate must be an explicit boolean. An absent, empty, or null gate is
    # missing evidence (handled in the census below), never a pass.
    adequacy_pass = (bank.get("adequacy_gate") or {}).get("pass")
    if adequacy_pass is False:
        control_failures.append(
            "B6 (bank_manifest): template bank adequacy gate failed: "
            + "; ".join((bank.get("adequacy_gate") or {}).get("failures")
                        or ["unspecified"]))
    broken_validations = []
    for name, v in validations:
        ev = v.get("event", name)
        # A validation that could not be read, or that carries no criteria,
        # is not a passing event -- it is absent evidence.
        if "_read_error" in v or not v.get("criteria"):
            broken_validations.append(
                f"{name}: unreadable or carries no criteria")
            continue
        # A missing tier must not demote hard failures to notes: default to
        # the strict tier rather than the lenient one.
        tier = v.get("tier") or "A"
        for c in v.get("criteria", []):
            if c.get("pass"):
                continue
            code = c["criterion"].split()[0]
            entry = f"{ev} {c['criterion']}"
            if tier == "A" and code in HARD_FAIL_CRITERIA:
                tier_a_hard.append(entry)
            elif tier == "A" and c.get("required"):
                tier_a_soft.append(entry)
            else:
                tier_b_notes.append(entry)

    # ---- evidence census: PASS must be earned, not defaulted ------------
    # SPEC.md §6.6 makes the V5 controls hard gates. A run whose control (or
    # validation) results simply never reached this job must therefore be
    # INCOMPLETE, never PASS: absent evidence is not passing evidence.
    missing_evidence = []
    if not args.controls_skipped:
        n_null = sum(1 for _, c in controls
                     if c.get("control") == "offsource_null")
        inj_masses = {(c.get("injected") or {}).get("chirp_mass_detector_msun")
                      for _, c in controls
                      if c.get("control") == "injection_recovery"}
        inj_masses.discard(None)
        if n_null < 1:
            missing_evidence.append(
                "no offsource_null control result supplied (SPEC.md V5a)")
        if len(inj_masses) < args.min_injections:
            missing_evidence.append(
                f"injection_recovery controls cover only {len(inj_masses)} "
                f"distinct chirp mass(es) of the {args.min_injections} "
                f"SPEC.md V5b requires")
    if not validations:
        missing_evidence.append("no per-event validation results supplied")
    missing_evidence.extend(
        f"validation not usable: {b}" for b in broken_validations)
    if adequacy_pass is not True and adequacy_pass is not False:
        missing_evidence.append(
            "bank manifest carries no usable adequacy gate result: pass must "
            "be an explicit boolean (SPEC.md B6)")
    # An empty expectation list would make the coverage census vacuous, so
    # its absence is itself missing evidence. The generator always passes
    # the full event list; a standalone caller must state one too.
    if not args.expected_event:
        missing_evidence.append(
            "no --expected-event list supplied, so event coverage cannot "
            "be verified")
    validated_events = {v.get("event") for _, v in validations
                        if "_read_error" not in v and v.get("criteria")}
    for ev in args.expected_event:
        if ev not in validated_events:
            missing_evidence.append(
                f"expected event {ev} has no usable validation result")

    if args.controls_skipped:
        verdict = "INCOMPLETE"
    elif fallback:
        verdict = "INVALID_FALLBACK_WAVEFORM"
    elif control_failures or tier_a_hard:
        verdict = "FAIL"
    elif missing_evidence:
        verdict = "INCOMPLETE"
    elif tier_a_soft or tier_b_notes:
        verdict = "PASS_WITH_NOTES"
    else:
        verdict = "PASS"

    summary = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "verdict": verdict,
        "n_events_validated": len(validations),
        "n_controls": len(controls),
        "controls_skipped": bool(args.controls_skipped),
        "missing_evidence": missing_evidence,
        "waveform_fallback_used": fallback,
        "control_failures": control_failures,
        "tier_a_hard_failures": tier_a_hard,
        "tier_a_other_failures": tier_a_soft,
        "tier_b_notes": tier_b_notes,
        "bank": {k: bank.get(k) for k in
                 ("n_templates", "chirp_mass_detector_range", "approximant",
                  "minimal_match_requested", "minimal_match_measured_min",
                  "waveform_fallback_used", "truncated")},
        # Only an explicit boolean is a gate result; anything else is
        # recorded as None so no consumer can mistake it for a pass.
        "bank_adequacy_pass": (adequacy_pass
                               if isinstance(adequacy_pass, bool) else None),
        "events": [
            {"event": v.get("event"), "tier": v.get("tier"),
             "verdict": v.get("verdict"), "failures": v.get("failures", []),
             "recovered": {
                 c["criterion"].split()[0]: c.get("recovered")
                 for c in v.get("criteria", [])}}
            for _, v in validations],
        "controls": [{"file": n, **{k: c.get(k) for k in
                                    ("control", "criterion", "pass",
                                     "failure_reason")}} for n, c in controls],
        "not_established": [
            "false alarm rate / significance (SPEC.md N1)",
            "credible intervals (SPEC.md N3)",
            "component masses, spins, distance, sky position (SPEC.md N4-N7)",
        ],
    }
    with open(args.output_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- markdown ------------------------------------------------------
    L = []
    w = L.append
    w("# gw-workflow — results")
    w("")
    w(f"**Verdict: {verdict}**")
    w("")
    w(f"Generated {summary['generated_utc']} · {len(validations)} event(s) "
      f"validated · {len(controls)} control(s)")
    w("")
    if verdict == "INVALID_FALLBACK_WAVEFORM":
        w("> ⚠️ **The restricted-PN fallback waveform was used, not the requested "
          "approximant.** These numbers are a smoke test of the plumbing and "
          "carry no physics claim. Install PyCBC in the container and re-run "
          "(SPEC.md M2).")
        w("")
    if args.controls_skipped:
        w("> ⚠️ **Controls were skipped**, so the run is INCOMPLETE by "
          "construction and can never be PASS (SPEC.md §7.3).")
        w("")
    if missing_evidence and not args.controls_skipped:
        w("> ⚠️ **Required evidence is missing**, so this run cannot be "
          "PASS (SPEC.md §6.6):")
        for m in missing_evidence:
            w(f"> - {m}")
        w("")

    w("## Controls (hard gates)")
    w("")
    if not controls:
        w("_No control results supplied._")
    else:
        w("| Criterion | Control | Result | Detail |")
        w("|---|---|---|---|")
        for name, c in controls:
            res = "PASS" if c.get("pass") else "**FAIL**"
            detail = c.get("failure_reason") or (
                f"on-source {c.get('onsource_network_snr')} vs loudest off-source "
                f"{c.get('loudest_offsource_network_snr')}"
                if c.get("control") == "offsource_null" else
                f"Mc err {(c.get('recovered') or {}).get('chirp_mass_fractional_error')}, "
                f"dt {(c.get('recovered') or {}).get('time_error_ms')} ms")
            w(f"| {c.get('criterion','?')} | {c.get('control','?')} | {res} | {detail} |")
    w("")
    w("A control failure is a hard fail regardless of how well the event "
      "comparisons went: a pipeline that cannot pass its own null and injection "
      "tests has not validated anything (SPEC.md §6.6).")
    w("")

    w("## Template bank")
    w("")
    b = summary["bank"]
    w(f"- Templates: **{b.get('n_templates')}**, detector-frame chirp mass "
      f"{b.get('chirp_mass_detector_range')}")
    w(f"- Approximant: `{b.get('approximant')}`; minimal match requested "
      f"{b.get('minimal_match_requested')}, measured minimum "
      f"{b.get('minimal_match_measured_min')}")
    gate = summary["bank_adequacy_pass"]
    w(f"- B6 adequacy gate: "
      f"**{'PASS' if gate is True else 'FAIL' if gate is False else 'NOT RECORDED (no usable gate result)'}**")
    if b.get("truncated"):
        w("- ⚠️ Bank hit `--max-templates`; coverage is incomplete.")
    w("")

    w("## Events")
    w("")
    w("| Event | Tier | Verdict | Network SNR | Mc det (Msun) | dt (ms) | Failures |")
    w("|---|---|---|---|---|---:|---|")
    for _, v in validations:
        got = {c["criterion"].split()[0]: c for c in v.get("criteria", [])}
        snr = (got.get("V1a") or {}).get("recovered")
        mc = (got.get("V2a/V2b") or {}).get("recovered")
        dt = (got.get("V3b") or {}).get("recovered")
        fails = "; ".join(v.get("failures") or []) or "—"
        w(f"| {v.get('event')} | {v.get('tier')} | {v.get('verdict')} | "
          f"{snr if snr is not None else '—'} | {mc if mc is not None else '—'} | "
          f"{dt if dt is not None else '—'} | {fails} |")
    w("")

    if conditioning:
        w("## Conditioning")
        w("")
        w("| Event | Det | Rate in→out | PSD blocks | Gates | Off-source spans |")
        w("|---|---|---|---:|---:|---:|")
        for _, c in conditioning:
            w(f"| {c.get('event')} | {c.get('detector')} | "
              f"{c.get('sample_rate_in')}→{c.get('sample_rate_out')} | "
              f"{(c.get('psd') or {}).get('n_blocks')} | "
              f"{(c.get('analysis') or {}).get('gates_applied')} | "
              f"{(c.get('offsource') or {}).get('placed')}"
              f"/{(c.get('offsource') or {}).get('requested')} |")
        w("")

    if spectrograms:
        w("## Spectrograms")
        w("")
        w("Non-blocking (SPEC.md V4). Peak normalised energy and its location:")
        w("")
        w("| Event | Det | Peak | Freq (Hz) | t−t_event (s) | File |")
        w("|---|---|---:|---:|---:|---|")
        for _, s in spectrograms:
            w(f"| {s.get('event')} | {s.get('detector')} | "
              f"{s.get('peak_normalised_energy')} | {s.get('peak_frequency_hz')} | "
              f"{s.get('peak_time_from_event_s')} | `{s.get('png')}` |")
        w("")

    w(NON_EXPECTATIONS)
    with open(args.output_md, "w") as fh:
        fh.write("\n".join(L) + "\n")

    logger.info("Verdict: %s", verdict)
    for msg in control_failures + tier_a_hard:
        logger.error("%s", msg)
    logger.info("Output: %s, %s", args.output_md, args.output_json)


if __name__ == "__main__":
    main()
