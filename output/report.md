# gw-workflow — results

**Verdict: PASS**

Generated 2026-08-23T02:29:43Z · 7 event(s) validated · 6 control(s)

## Controls (hard gates)

| Criterion | Control | Result | Detail |
|---|---|---|---|
| V5a | offsource_null | PASS | on-source 23.6942 vs loudest off-source 0.0 |
| V5b | injection_recovery | PASS | Mc err 0.00227, dt 0.0 ms |
| V5b | injection_recovery | PASS | Mc err 0.00256, dt 0.0 ms |
| V5b | injection_recovery | PASS | Mc err 0.0021, dt 0.0 ms |
| V5b | injection_recovery | PASS | Mc err 0.05037, dt 0.977 ms |
| V5b | injection_recovery | PASS | Mc err 0.01764, dt 0.488 ms |

A control failure is a hard fail regardless of how well the event comparisons went: a pipeline that cannot pass its own null and injection tests has not validated anything (SPEC.md §6.6).

## Template bank

- Templates: **3858**, detector-frame chirp mass [5.0, 59.94094156401742]
- Approximant: `IMRPhenomD`; minimal match requested 0.97, measured minimum 0.96805
- B6 adequacy gate: **PASS**

## Events

| Event | Tier | Verdict | Network SNR | Mc det (Msun) | dt (ms) | Failures |
|---|---|---|---|---|---:|---|
| GW150914 | A | PASS | 23.6942 | 31.57544 | 7.324 | — |
| GW151226 | A | PASS | 12.2649 | 9.47313 | 1.074 | — |
| GW170104 | A | PASS | 13.5263 | 26.0222 | 2.93 | — |
| GW170608 | A | PASS | 19.1436 | 8.15802 | 6.836 | — |
| GW170809 | A | PASS | 12.4562 | 29.44581 | 8.789 | — |
| GW170814 | A | PASS | 16.9908 | 26.2131 | 7.812 | — |
| GW170823 | A | PASS | 11.5105 | 38.92656 | 0.977 | — |

## Conditioning

| Event | Det | Rate in→out | PSD blocks | Gates | Off-source spans |
|---|---|---|---:|---:|---:|
| GW150914 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW150914 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW151226 | H1 | 4096→2048 | 127 | 2 | 7/10 |
| GW151226 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170104 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170104 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170608 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170608 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170809 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170809 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170814 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170814 | L1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170823 | H1 | 4096→2048 | 127 | 2 | 10/10 |
| GW170823 | L1 | 4096→2048 | 127 | 2 | 10/10 |

## Spectrograms

Non-blocking (SPEC.md V4). Peak normalised energy and its location:

| Event | Det | Peak | Freq (Hz) | t−t_event (s) | File |
|---|---|---:|---:|---:|---|
| GW150914 | H1 | 86.139 | 53.19 | -0.0122 | `GW150914_H1_spectrogram.png` |
| GW150914 | L1 | 33.369 | 93.27 | 0.0088 | `GW150914_L1_spectrogram.png` |
| GW151226 | H1 | 17.162 | 55.16 | -0.2641 | `GW151226_H1_spectrogram.png` |
| GW151226 | L1 | 15.895 | 163.54 | 0.0298 | `GW151226_L1_spectrogram.png` |
| GW170104 | H1 | 24.065 | 102.11 | -0.0122 | `GW170104_H1_spectrogram.png` |
| GW170104 | L1 | 27.786 | 100.28 | -0.0083 | `GW170104_L1_spectrogram.png` |
| GW170608 | H1 | 1094.312 | 31.46 | 0.9976 | `GW170608_H1_spectrogram.png` |
| GW170608 | L1 | 18.263 | 122.39 | -0.0391 | `GW170608_L1_spectrogram.png` |
| GW170809 | H1 | 24.458 | 63.76 | -0.0576 | `GW170809_H1_spectrogram.png` |
| GW170809 | L1 | 17.99 | 50.38 | -0.105 | `GW170809_L1_spectrogram.png` |
| GW170814 | H1 | 21.559 | 55.16 | -0.0132 | `GW170814_H1_spectrogram.png` |
| GW170814 | L1 | 43.287 | 89.95 | 0.0098 | `GW170814_L1_spectrogram.png` |
| GW170823 | H1 | 13.453 | 502.81 | -0.6162 | `GW170823_H1_spectrogram.png` |
| GW170823 | L1 | 28.906 | 76.42 | -0.0029 | `GW170823_L1_spectrogram.png` |

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

