# Reference summary — transient extraction from open LIGO/Virgo strain data

Summary of the four papers in this directory, focused on (1) how a transient
gravitational-wave signal is extracted from open strain data, and (2) which
published values are quoted per event.

Machine-readable companion: [`gwtc1_reference.json`](gwtc1_reference.json).

| File | Paper | Role |
|---|---|---|
| `1-s2.0-S2352711021000030-main.pdf` | Abbott et al. 2021, *SoftwareX* **13**, 100658 — GWOSC open data for O1/O2 | what you download |
| `Usman_2016_Class._Quantum_Grav._33_215004.pdf` | Usman et al. 2016, *CQG* **33**, 215004 — the PyCBC search pipeline | the algorithm |
| `PhysRevX.9.031040.pdf` | Abbott et al. 2019, *PRX* **9**, 031040 — GWTC-1 | the numbers |
| `PhysRevLett.116.061102.pdf` | Abbott et al. 2016, *PRL* **116**, 061102 — GW150914 | the founding observation |

**Coverage gap:** these four cover O1 and O2 only. For anything after GW170823,
use GWTC-2 / GWTC-2.1 / GWTC-3.

---

## Part 1 — How a transient is extracted

### Layer 0: what the open data actually is (SoftwareX)

One data product: calibrated dimensionless strain `h(t)`, derived from the
interferometer's differential arm length readout.

- **Sample rates.** LIGO DAQ at **16 384 Hz**, Virgo at **20 000 Hz**, digitally
  converted to the LIGO rate before analysis. GWOSC serves the native 16 kHz plus
  a **4096 Hz** decimation (`scipy.signal.decimate`, order-8 Chebyshev-I IIR
  anti-alias filter). Because of the filter roll-off, **4 kHz data is only valid
  to about 1700 Hz** — not the 2048 Hz Nyquist frequency.
- **Usable band.** `h(t)` is calibrated only over **10–5000 Hz** (LIGO) and
  **10–8000 Hz** (Virgo). Outside those ranges it is not a faithful strain
  measurement. The data is **high-passed at 8 Hz** to avoid dynamic-range and
  floating-point-precision problems downstream.
- **Packaging.** Bulk files of **4096 s**; per-event snippets of **32 s and
  4096 s**; formats `hdf5`, `gwf`, and two-column ASCII (GPS time, strain).
- **Provenance channels** (Table 1):

  | Run | Detector | Channel |
  |---|---|---|
  | O1 | H1 / L1 | `H1:DCS-CALIB_STRAIN_C02` / `L1:DCS-CALIB_STRAIN_C02` |
  | O2 | H1 / L1 | `H1:DCH-CLEAN_STRAIN_C02` / `L1:DCH-CLEAN_STRAIN_C02` |
  | O2 | V1 | `V1:Hrec_hoft_V1O2Repro2A_16384Hz` |

  `CLEAN` marks witness-channel noise subtraction. Downsampling and NaN
  substitution are the **only** modifications GWOSC makes to the archive data.
- **Data quality and injections** arrive as **1 Hz bitmasks**: `DATA`,
  `CBC_CAT1/2/3`, `BURST_CAT1/2/3`, plus five injection bits (CBC, burst,
  detchar, CW, stochastic). Categories cascade — failing one fails all higher
  ones. Seconds failing CAT1 are **replaced by NaN in the strain array**.
- **Dead times.** CAT1: 1.683 % (H1) / 1.039 % (L1) in O1; 0.001 / 0.003 /
  0.053 % (H1/L1/V1) in O2. CBC CAT2: 0.890 / 0.007 % (O1), 0.157 / 0.090 % (O2).
  BURST CAT2: 0.624 / 0.021 % (O1), 0.212 / 0.151 % (O2) — not used for Virgo.
  CAT3 was unused in O1/O2 but is retained in the file format.
- **Explicit caveat from the paper:** applying every data-quality flag does *not*
  mean the surviving data is free of transient artifacts. Virgo performed no
  hardware injections in O2, so all its injection bits read 1.

### Layer 1: conditioning (PyCBC §3.1–3.2, GWTC-1 App. A)

**PSD estimation.** PyCBC decouples the PSD from the filter block — the main
departure from the older FindChirp scheme. The inverse PSD is truncated to
**16 s** in the time domain to bound the filter's impulse response, so it is
natural to estimate on 16 s blocks: a **2048 s** input stretch gives **127
averages** (against ~15 overlapping blocks before), sharply cutting PSD variance.
The estimate is then interpolated to the data resolution, inverted, inverse-FFT'd,
truncated to 16 s, and FFT'd back for use in the filter.

Template-bank *placement* uses a different PSD: the **harmonic mean** of 2048 s
median estimates, averaged over time *and* over detectors, so a single bank
serves the whole run and every detector.

**Gating** is the core glitch-removal step, and the paper's motivation for it is
worth internalising: a delta-function glitch convolves with both the template
*and* the inverse PSD, so one glitch produces SNR ~1000 triggers at its own time
plus a skirt of SNR ~10 triggers spread over **±8 s**. The chi-squared
reweighting suppresses but does not remove these. So PyCBC finds the transient —
whiten in the frequency domain, inverse-FFT, threshold `|whitened strain|`,
cluster the peaks — and **zeroes about 1 s of data with a Tukey window, shape
parameter α = 0.5**. Applied to short (~1 s) transients only; longer ones are
removed by data-quality vetoes before the pipeline runs.

GWTC-1 App. A adds the practical picture. Four glitch classes matter: scattered
light, 60–200 Hz nonstationarity, short-duration transients, and blips.
Short-duration DAC-overflow transients occur at **~1/hour** in each LIGO detector
and are gated by the searches — **the LLO glitch during GW170817 is of this
class**. cWB does not gate, which is *why cWB rejected GW170817*. For parameter
estimation GW170817 used full **BayesWave glitch subtraction** instead, which is
more robust than gating. Blip transients (~1/hour, sub-second, band-limited)
remain the limiting noise source for high-mass CBC searches in O1/O2 — they do
not couple into any witness sensor, so only signal-morphology tests mitigate them.

**Noise subtraction.** O2 LIGO data has witness-channel linear-coupling noise
subtracted post-facto: **+~20 % BNS range at LHO** averaged over O2 (traced to a
point absorber on an input test mass coupling laser-table jitter into the
interferometer), negligible at LLO where only narrow line features were targeted.
The SoftwareX paper phrases the same result as up to **+30 % sensitive volume**.

### Layer 2: matched filtering (PyCBC §3.3)

`rho^2(t)` is built from two orthogonal template phases (`h_cos`, `h_sin`), which
maximises over the binary's sky position and orientation analytically. Concrete
numbers: `dt = 1/4096 s`, FFT blocks `T_B = 256 s`, so `N = 2^20`. The bank is
placed for at most **3 %** SNR loss from its discreteness. Triggers are local
maxima above threshold; clustering divides the 256 s SNR series into **1 s
windows**, takes the max in each (parallelisable), then keeps a window's maximum
only if it beats both neighbours.

Bank and threshold configurations differ by paper:

| | GW150914 PRL (2016) | GWTC-1 PyCBC | GWTC-1 GstLAL |
|---|---|---|---|
| SNR threshold | — | 5.5 | 4 (LIGO), 3.5 (Virgo) |
| Component / total mass | 1–99 Msun, M < 100 Msun | M_tot 2–500 Msun | M_tot 2–400 Msun |
| Mass ratio | — | down to 1/98 | 1/98 – 1 |
| Max spin | 0.99 | 0.05 if m < 2 Msun, else 0.998 | 0.05 if m < 2 Msun, else 0.999 |
| Templates | ~250 000 | — | — |
| Waveforms | EOB aligned-spin | SEOBNRv4 ROM > 4 Msun, TaylorF2 below | TaylorF2 if chirp mass < 1.73, else SEOBNRv4 ROM |
| Analyses Virgo | n/a | no | yes (August 2017) |

GWTC-1's PyCBC bank also requires **waveform duration ≥ 0.15 s** — this is what
sets the high-mass boundary, and it exists to suppress false alarms from short
glitches.

### Layer 3: signal consistency (PyCBC §3.4)

Chi-squared splits the template into `p` frequency bins of equal expected power
and compares measured to expected SNR per bin. `p` varies with template
parameters (more bins for longer, lower-mass templates). With
`chi2_r = chi2/(2p-2)`, the SNR is reweighted as

    rho_hat = rho / {[1 + (chi2_r)^3] / 2}^(1/6)     for chi2_r > 1
    rho_hat = rho                                     otherwise

PyCBC's efficiency trick: evaluate chi-squared only at clustered triggers via an
optimised frequency integral rather than `p` extra FFTs, with an iterative
multiply for the exponential term. The crossover is `N_p ~ (5/7) p log N ~ 100`
points for `N = 2^20`, and real LIGO trigger rates sit below that.

### Layer 4: coincidence (PyCBC §3.5)

A **15 ms** window = 10 ms LHO–LLO light travel + 5 ms timing slop, *and*
**exact-match** on template (identical masses and spins) — only possible because
a single bank is used everywhere. This replaced the older metric-based test and
measurably improved sensitivity. One consequence: triggers can no longer be
clustered across the bank before coincidence, or quiet signals are lost when
noise shifts the best-fit template between detectors.

Ranking is `rho_hat_c = sqrt(rho_hat_1^2 + rho_hat_2^2)`. GWTC-1 uses more
elaborate per-pipeline statistics — a single-detector rank plus a network
likelihood-ratio approximation for PyCBC; a full log-likelihood-ratio for GstLAL
that folds in SNR, signal consistency, time and phase differences, which
detectors participated, per-detector sensitivity, and per-detector trigger rate.

### Layer 5: significance (PyCBC §3.6, GWTC-1 §III)

Time shifts with an offset larger than the coincidence window give background
time `T_b = T^2 / delta`. With **delta = 0.1 s**, about **15 days** of coincident
data reaches FAR ~1 in 2e5 yr, i.e. **5 sigma** — the paper states outright that
this is why 16 days was the GW150914 analysis length. The p-value is

    p_b(rho_hat_c; A) = 1 - exp[ -n_C * T * (1 + n_b^A(rho_hat_c)) / T_b ]

with a **trials factor `n_C` = number of search classes** (3, split by template
length). A final coincident clustering over a **10 s** window keeps candidates
statistically independent. Background triggers below a statistic threshold are
decimated by a stored factor to save disk.

GWTC-1 tightens this: analysis periods each containing at least **5.2 days** of
coincident LIGO data, and shifts performed **only within a period** (detector
noise changed too much across O1 to O2), which raises the floor on the smallest
measurable FAR. Event selection then runs in two stages — FAR < **1 per 30 days**
in at least one matched-filter search to be listed at all, then **p_astro > 0.5**
for the "GW" designation, computed from a four-category (terrestrial / BNS /
NSBH / BBH) Poisson mixture fit.

### Layer 6: the unmodelled path (cWB)

Coherent excess power in multiresolution time-frequency maps, up to **1 kHz** and
durations up to a few seconds. Detection statistic
`eta_c = sqrt(2 Ec / (1 + En/Ec))`, where `Ec` is coherent and `En` residual
energy. Veto on the network correlation coefficient
`cc = Ec / (Ec + En) < 0.7`, plus a "chirping" cut requiring frequency to
increase with time. GW150914's PRL classes were C1 (known noise morphology), C3
(rising frequency), C2 (the remainder).

### Layer 7: parameter estimation (GWTC-1 App. B)

`LALInference` (MCMC and nested sampling), with `RapidPE` for expensive models.
PSD from **BayesWave** — cubic spline for broadband structure plus a sum of
Lorentzians for line features, median taken over the PSD posterior.
Frequency-dependent **spline calibration envelopes** are marginalised over (O1
had used uniform-in-frequency calibration priors).

- Likelihood low-frequency cutoff **20 Hz**, with exceptions: **GW170817 at
  23 Hz, GW170818 at 16 Hz, GW170608 Hanford at 30 Hz**.
- Spin quantities quoted at reference frequency **20 Hz**.
- Cosmology: flat LCDM, **H0 = 67.9 km/s/Mpc, Omega_m = 0.306**. Switching to
  Planck 2018 changes redshift by < 0.3 % and source-frame total mass by < 0.03 %.
- BBH posteriors are the **average of IMRPhenomPv2 and SEOBNRv3**; GW170817 uses
  IMRPhenomPv2NRT / SEOBNRv4NRT / TaylorF2 / SEOBNRv4T / TEOBResumS.
- O2 final calibration uncertainties: **3.8 % amplitude / 2.1 deg phase (LLO)**,
  **2.6 % / 2.4 deg (LHO)**, **5.1 % / 2.3 deg (Virgo, plus 20 us timestamp)**.
  LIGO timing uncertainty < 1 us.

---

## Part 2 — Published values quoted per event

**GWTC-1 Table III is the authoritative set and explicitly supersedes every
earlier per-event paper.** Values are medians with symmetric 90 % credible
intervals, including statistical error, waveform-model systematics, and
calibration uncertainty. Masses are source-frame.

| Event | m1/Msun | m2/Msun | chirp M/Msun | chi_eff | M_f/Msun | a_f | E_rad/Msun c^2 | l_peak (erg/s) | d_L/Mpc | z | sky area/deg^2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GW150914 | 35.6 +4.7 -3.1 | 30.6 +3.0 -4.4 | 28.6 +1.7 -1.5 | -0.01 +0.12 -0.13 | 63.1 +3.4 -3.0 | 0.69 +0.05 -0.04 | 3.1 +0.4 -0.4 | 3.6 +0.4 -0.4 e56 | 440 +150 -170 | 0.09 +0.03 -0.03 | 182 |
| GW151012 | 23.2 +14.9 -5.5 | 13.6 +4.1 -4.8 | 15.2 +2.1 -1.2 | 0.05 +0.31 -0.20 | 35.6 +10.8 -3.8 | 0.67 +0.13 -0.11 | 1.6 +0.6 -0.5 | 3.2 +0.8 -1.7 e56 | 1080 +550 -490 | 0.21 +0.09 -0.09 | 1523 |
| GW151226 | 13.7 +8.8 -3.2 | 7.7 +2.2 -2.5 | 8.9 +0.3 -0.3 | 0.18 +0.20 -0.12 | 20.5 +6.4 -1.5 | 0.74 +0.07 -0.05 | 1.0 +0.1 -0.2 | 3.4 +0.7 -1.7 e56 | 450 +180 -190 | 0.09 +0.04 -0.04 | 1033 |
| GW170104 | 30.8 +7.3 -5.6 | 20.0 +4.9 -4.6 | 21.4 +2.2 -1.8 | -0.04 +0.17 -0.21 | 48.9 +5.1 -4.0 | 0.66 +0.08 -0.11 | 2.2 +0.5 -0.5 | 3.3 +0.6 -1.0 e56 | 990 +440 -430 | 0.20 +0.08 -0.08 | 921 |
| GW170608 | 11.0 +5.5 -1.7 | 7.6 +1.4 -2.2 | 7.9 +0.2 -0.2 | 0.03 +0.19 -0.07 | 17.8 +3.4 -0.7 | 0.69 +0.04 -0.04 | 0.9 +0.0 -0.1 | 3.5 +0.4 -1.3 e56 | 320 +120 -110 | 0.07 +0.02 -0.02 | 392 |
| GW170729 | 50.2 +16.2 -10.2 | 34.0 +9.1 -10.1 | 35.4 +6.5 -4.8 | 0.37 +0.21 -0.25 | 79.5 +14.7 -10.2 | 0.81 +0.07 -0.13 | 4.8 +1.7 -1.7 | 4.2 +0.9 -1.5 e56 | 2840 +1400 -1360 | 0.49 +0.19 -0.21 | 1041 |
| GW170809 | 35.0 +8.3 -5.9 | 23.8 +5.1 -5.2 | 24.9 +2.1 -1.7 | 0.08 +0.17 -0.17 | 56.3 +5.2 -3.8 | 0.70 +0.08 -0.09 | 2.7 +0.6 -0.6 | 3.5 +0.6 -0.9 e56 | 1030 +320 -390 | 0.20 +0.05 -0.07 | 308 |
| GW170814 | 30.6 +5.6 -3.0 | 25.2 +2.8 -4.0 | 24.1 +1.4 -1.1 | 0.07 +0.12 -0.12 | 53.2 +3.2 -2.4 | 0.72 +0.07 -0.05 | 2.7 +0.4 -0.3 | 3.7 +0.4 -0.5 e56 | 600 +150 -220 | 0.12 +0.03 -0.04 | 87 |
| GW170817 | 1.46 +0.12 -0.10 | 1.27 +0.09 -0.09 | 1.186 +0.001 -0.001 | 0.00 +0.02 -0.01 | <= 2.8 | <= 0.89 | >= 0.04 | >= 0.1 e56 | 40 +7 -15 | 0.01 | 16 |
| GW170818 | 35.4 +7.5 -4.7 | 26.7 +4.3 -5.2 | 26.5 +2.1 -1.7 | -0.09 +0.18 -0.21 | 59.4 +4.9 -3.8 | 0.67 +0.07 -0.08 | 2.7 +0.5 -0.5 | 3.4 +0.5 -0.7 e56 | 1060 +420 -380 | 0.21 +0.07 -0.07 | 39 |
| GW170823 | 39.5 +11.2 -6.7 | 29.0 +6.7 -7.8 | 29.2 +4.6 -3.6 | 0.09 +0.22 -0.26 | 65.4 +10.1 -7.4 | 0.72 +0.09 -0.12 | 3.3 +1.0 -0.9 | 3.6 +0.7 -1.1 e56 | 1940 +970 -900 | 0.35 +0.15 -0.15 | 1666 |

GW170817's remnant entries are deliberately **conservative bounds, not medians**:
`E_rad` and `l_peak` are lower limits at 1 % credibility, `M_f` an upper limit at
99 %, and `a_f` follows from a phenomenological universal relation. Source-frame
masses for GW170817 use the NGC 4993 redshift; for the BBHs they use the LCDM
conversion from luminosity distance. Sky area is the 90 % credible region.

### Detection statistics per event (GWTC-1 Table I)

FAR is per-pipeline; a blank means that pipeline did not report the event.

| Event | UTC | FAR PyCBC (1/yr) | FAR GstLAL | FAR cWB | SNR PyCBC / GstLAL / cWB |
|---|---|---|---|---|---|
| GW150914 | 09:50:45.4 | < 1.53e-5 | < 1.00e-7 | < 1.63e-4 | 23.6 / 24.4 / 25.2 |
| GW151012 | 09:54:43.4 | 0.17 | 7.92e-3 | — | 9.5 / 10.0 / — |
| GW151226 | 03:38:53.6 | < 1.69e-5 | < 1.00e-7 | 0.02 | 13.1 / 13.1 / 11.9 |
| GW170104 | 10:11:58.6 | < 1.37e-5 | < 1.00e-7 | 2.91e-4 | 13.0 / 13.0 / 13.0 |
| GW170608 | 02:01:16.5 | < 3.09e-4 | < 1.00e-7 | 1.44e-4 | 15.4 / 14.9 / 14.1 |
| GW170729 | 18:56:29.3 | 1.36 | 0.18 | 0.02 | 9.8 / 10.8 / 10.2 |
| GW170809 | 08:28:21.8 | 1.45e-4 | < 1.00e-7 | — | 12.2 / 12.4 / — |
| GW170814 | 10:30:43.5 | < 1.25e-5 | < 1.00e-7 | < 2.08e-4 | 16.3 / 15.9 / 17.2 |
| GW170817 | 12:41:04.4 | < 1.25e-5 | < 1.00e-7 | — | 30.9 / 33.0 / — |
| GW170818 | 02:25:09.1 | — | 4.20e-5 | — | — / 11.3 / — |
| GW170823 | 13:13:58.5 | < 3.29e-5 | < 1.00e-7 | 2.14e-3 | 11.1 / 11.5 / 10.8 |

Three cases matter if you are reproducing anything:

- **GW170818** was found by GstLAL alone, as a triple coincidence (SNR: H 4.1,
  L 9.7, V 4.2). PyCBC saw only a "chirp-like" single-detector Livingston trigger.
- **GW170809** and **GW170817** were dropped by cWB — the former failed a
  signal-independent veto in the final noise-subtracted analysis, the latter
  because of the un-gated glitch.
- **GW170608**'s GstLAL FAR is assigned from its own local background alone,
  because its analysis period is unique.

Probability of astrophysical origin (GWTC-1 Table IV, `astrophysical` column):
every confident event is 1.0 or near it in GstLAL, but **GW170729 is the weakest
of the eleven** — 0.98 in GstLAL and only **0.52 in PyCBC**. The full
four-category breakdown per pipeline is in the JSON, along with all 14 marginal
triggers from Table II (FAR, network SNR, detector-frame chirp mass, and the
detector-characterisation note per trigger).

### Where the earlier per-event values differ

The 2016 PRL quotes for GW150914: **m1 = 36 +5 -4, m2 = 29 +4 -4,
M_f = 62 +4 -4 Msun, a_f = 0.67 +0.05 -0.07, E_rad = 3.0 +0.5 -0.5 Msun c^2,
l_peak = 3.6 +0.5 -0.4 e56 erg/s (= 200 +30 -20 Msun c^2/s),
d_L = 410 +160 -180 Mpc, z = 0.09 +0.03 -0.04**, sky area ~600 deg^2, combined
network SNR 24, `rho_hat_c` = 23.6, `eta_c` = 20.0, arrival 6.9 +0.5 -0.4 ms
after L1, peak strain 1.0e-21, frequency sweep 35 to 250 Hz, primary spin
< 0.7 (90 % credible). Significance: FAR < 1 per 203 000 yr (**5.1 sigma**) from
the CBC search, < 1 per 22 500 yr (**4.6 sigma**) from the burst search in class
C3, and 1 per 8400 yr (4.4 sigma) in the unconstrained C2+C3 class. The second
most significant event of that analysis had FAR 1 per 2.3 yr, p = 0.02.

GWTC-1 §V F itemises the revisions and their causes:

| Event | Total-mass median shift | Note |
|---|---|---|
| GW150914 | +~1.0 Msun | improved PSD estimation, spline calibration envelopes |
| GW170104 | +0.3 Msun | recalibration and noise subtraction |
| GW151012 | +0.2 Msun | |
| GW151226 | -0.3 Msun | |
| GW170608 | -0.2 Msun | |

GW170104 benefited most from noise subtraction: PE matched-filter SNR rose from
**13.3 +0.2 -0.3 to 14.0 +0.2 -0.3**, and `chi_eff` moved from
**-0.12 +0.21 -0.30 to -0.04 +0.17 -0.21** (the low-value tail was lost, upper
limit roughly unchanged). For GW151226, the fully precessing model gives
`chi_eff = 0.15 +0.25 -0.11` and the effective-precession model
`0.20 +0.18 -0.08`, with at least one spin magnitude > 0.28 at 99 % probability,
against 0.2 in the O1 analysis. GW170817's tidal deformability moved from the
published `300 +420 -230` to `330 +438 -251` (IMRPhenomPv2NRT, low-spin prior,
90 % highest-posterior-density interval); the other four models give 305, 394,
349 and 405.

### Population values (GWTC-1)

- 10 BBH + 1 BNS confident, 0 NSBH, 14 marginal triggers.
- Roughly **one detection per 15 days** of searched data.
- BBH total masses span **18.6 +3.2 -0.7** to **84.4 +15.8 -11.1 Msun**;
  distances **320 +120 -110** to **2840 +1400 -1360 Mpc**.
- Merger rates, 90 % CI: **BNS 110–3840**, **BBH 9.7–101 Gpc^-3 yr^-1**;
  **NSBH < 610 Gpc^-3 yr^-1** (90 % upper limit).
- Primary-mass power-law index **alpha = 1.3 +1.4 -1.7**.

### Run-level context

| | O1 | O2 |
|---|---|---|
| Dates | 2015-09-12 to 2016-01-19 | 2016-11-30 to 2017-08-25 (Virgo from 2017-08-01) |
| Duty factor | 64.6 % H1, 57.4 % L1 | 65.3 % H1, 61.8 % L1, 85.1 % V1 |
| Network duty factor | 42.8 % | 46.1 % (LIGO), 35 % (LIGO+Virgo) |
| Detections | 3 | 8 |

O2 yielded about **118 days** of data suitable for coincident LIGO analysis, of
which about **15 days** were coincident with Virgo.
