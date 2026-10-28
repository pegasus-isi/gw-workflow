# gw-workflow — results

**Verdict: PASS**

Generated 2026-08-23T01:01:11Z · 1 event(s) validated · 6 control(s)

## Controls (hard gates)

| Criterion | Control | Result | Detail |
|---|---|---|---|
| V5a | offsource_null | PASS | on-source 23.6136 vs loudest off-source 0.0 |
| V5b | injection_recovery | PASS | Mc err 0.00245, dt 0.488 ms |
| V5b | injection_recovery | PASS | Mc err 0.01765, dt 0.488 ms |
| V5b | injection_recovery | PASS | Mc err 0.00197, dt 0.0 ms |
| V5b | injection_recovery | PASS | Mc err 0.00357, dt 0.0 ms |
| V5b | injection_recovery | PASS | Mc err 0.02512, dt 0.0 ms |

A control failure is a hard fail regardless of how well the event comparisons went: a pipeline that cannot pass its own null and injection tests has not validated anything (SPEC.md §6.6).

## Template bank

- Templates: **49**, detector-frame chirp mass [25.0, 38.95468823410252]
- Approximant: `IMRPhenomD`; minimal match requested 0.97, measured minimum 0.96993
- B6 adequacy gate: **PASS**

## Events

| Event | Tier | Verdict | Network SNR | Mc det (Msun) | dt (ms) | Failures |
|---|---|---|---|---|---:|---|
| GW150914 | A | PASS | 23.6136 | 31.92637 | 6.836 | — |

## Conditioning

| Event | Det | Rate in→out | PSD blocks | Gates | Off-source spans |
|---|---|---|---:|---:|---:|
| GW150914 | H1 | 4096→2048 | 127 | 2 | 5/5 |
| GW150914 | L1 | 4096→2048 | 127 | 2 | 5/5 |

## Spectrograms

Non-blocking (SPEC.md V4). Peak normalised energy and its location:

| Event | Det | Peak | Freq (Hz) | t−t_event (s) | File |
|---|---|---:|---:|---:|---|
| GW150914 | H1 | 86.139 | 53.19 | -0.0122 | `GW150914_H1_spectrogram.png` |
| GW150914 | L1 | 33.369 | 93.27 | 0.0088 | `GW150914_L1_spectrogram.png` |

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

