# Validation record — GW150914, real GWOSC data

Output of an actual run of this workflow against real GWOSC strain, retained as
evidence. Reproduce with the commands in the parent `README.md`.

| | |
|---|---|
| Date | 2026-08-23 |
| Event | GW150914 (Tier A) |
| Data | GWOSC `GWTC-1-confident/GW150914/v3`, 4096 s @ 4096 Hz, H1 + L1 |
| Approximant | `IMRPhenomD` (real library model, `waveform_fallback_used: false`) |
| Bank | 49 templates, ℳ_det 25–39, measured minimal match 0.96993 |
| **Verdict** | **PASS** — all 9 criteria, both hard gates |

## Recovered vs published

| Quantity | Recovered | GWTC-1 / PRL | Agreement |
|---|---:|---:|---|
| Network SNR | 23.614 | 23.6 (Table I, PyCBC) | **+0.1%** |
| ℳ detector-frame (M☉) | 31.926 | 31.174 (Table III × (1+z)) | +2.4%, inside the 90% interval |
| Merger time | — | published event GPS | within **21 ms** |
| H1–L1 delay | 6.836 ms, L1 first | 6.9 +0.5 −0.4 ms, L1 first (PRL) | **0.064 ms** |
| Per-detector SNR | H1 19.53, L1 13.28 | — | |

## Hard gates

- **V5a off-source null**: all 5 real off-source spans produced **zero**
  coincident candidates; on-source 23.61. Clean separation.
- **V5b injection recovery**: 5 injections across ℳ_det 26–38, all recovered —
  worst chirp-mass error 2.5%, worst timing error 0.5 ms.

## Caveats

- The bank was restricted to ℳ_det 25–40 for run time. A full ℳ_det 5–60 bank is
  needed to reproduce the whole event set; the B6 adequacy gate confirms this
  bank covers GW150914 specifically (nearest template 0.22% away).
- One event, not the full Tier A set.
- **No statistical significance is claimed** — see the parent README and
  `SPEC.md` N1. `report.md` restates this.
