# SPEC — Independent reproduction of GWTC-1 binary black hole detections

**Project:** `gw-workflow`
**Goal:** Recover the binary black hole detections of GWTC-1 from publicly
released GWOSC strain data using a Pegasus workflow built from the published
methods, and validate the recovered signal-to-noise ratio, chirp mass, and
merger time against the values published in GWTC-1 for each event.

Reference material and published values: [`references/`](references/), with
per-event numbers in machine-readable form at
[`references/gwtc1_reference.json`](references/gwtc1_reference.json) and the
method summary at [`references/SUMMARY.md`](references/SUMMARY.md).

Sourcing rules for all implementation work are in [`CLAUDE.md`](CLAUDE.md) and
are binding.

---

## 1. Scope

### 1.1 What this workflow does

For each selected GWTC-1 event and each participating detector:

1. **Fetch** the public 4096 s calibrated strain segment and its data-quality
   and injection bitmasks from GWOSC.
2. **Verify** data quality — no NaNs in the analysis span, CBC category flags
   satisfied, no hardware injection overlapping the analysis span.
3. **Condition** the strain — decimate, high-pass, estimate the noise power
   spectral density, whiten, and discard filter-corrupted edges.
4. **Matched filter** the conditioned data against a template bank, producing a
   signal-to-noise ratio time series per template, then threshold and cluster to
   triggers.
5. **Combine detectors** into coincident candidates and a network SNR.
6. **Estimate** the detector-frame chirp mass and the merger time from the
   best-matching coincident template.
7. **Render** a time-frequency spectrogram of the whitened data around the event.
8. **Validate** the recovered SNR, chirp mass, and merger time against the
   GWTC-1 published values, per event, with explicit pass/fail criteria.
9. **Report** an aggregate table across all events.

### 1.2 Event set

Detector-frame chirp mass is `M_source * (1 + z)`; ISCO frequency is estimated as
`4400 Hz / (M_total_detector / Msun)`. Both are computed from GWTC-1 Table III
values, not measured.

| Tier | Event | ℳ_det (Msun) | M_tot,det (Msun) | f_ISCO (Hz) | PyCBC network SNR | chi_eff |
|---|---|---|---|---|---|---|
| A | GW150914 | 31.2 | 72.2 | 61 | 23.6 | -0.01 |
| A | GW151226 | 9.7 | 23.3 | 189 | 13.1 | 0.18 |
| A | GW170104 | 25.7 | 61.0 | 72 | 13.0 | -0.04 |
| A | GW170608 | 8.5 | 19.9 | 221 | 15.4 | 0.03 |
| A | GW170809 | 29.9 | 70.6 | 62 | 12.2 | 0.08 |
| A | GW170814 | 27.0 | 62.5 | 70 | 16.3 | 0.07 |
| A | GW170823 | 39.4 | 92.5 | 48 | 11.1 | 0.09 |
| B | GW151012 | 18.4 | 44.5 | 99 | 9.5 | 0.05 |
| B | GW170729 | 52.8 | 125.5 | **35** | 9.8 | **0.37** |
| B | GW170818 | 32.1 | 75.1 | 59 | **not found by PyCBC** | -0.09 |

**Tier A (7 events)** are the primary targets and are subject to the full
validation criteria in §6.

**Tier B (3 events)** are run but held to relaxed criteria, for reasons that are
properties of the events rather than of this workflow:

- **GW151012** has the lowest published SNR of the catalog (9.5 PyCBC) and a
  PyCBC astrophysical probability of only 0.96 against a terrestrial 0.04.
- **GW170729** is the worst case on three counts at once: the largest published
  `chi_eff` (0.37), which a non-spinning bank cannot match; the highest redshift
  (0.49); and an ISCO frequency of **35 Hz**, only 15 Hz above the analysis
  low-frequency cutoff, so very little of the signal is in band. Its PyCBC
  astrophysical probability is 0.52 — the catalog's weakest.
- **GW170818** was found by **GstLAL alone**, as a triple coincidence with
  Hanford SNR 4.1, Livingston 9.7, Virgo 4.2. PyCBC did not report it; a
  two-detector Hanford-Livingston coincidence search of the kind specified here
  is not expected to recover it. It is included as an honest negative: recovering
  it would be a bonus, failing to recover it reproduces the published result.

**GW170817 is excluded.** It is a binary neutron star, requires tidal terms in
the waveform model, a different bank, and BayesWave glitch subtraction of the
Livingston transient. It is out of scope (§4).

### 1.3 Detectors

H1 and L1 for every event. Virgo (V1) is **optional and per-event**: Virgo joined
O2 on 2017-08-01 with a BNS range near 25 Mpc against 80-100 Mpc for the LIGO
instruments, so its contribution to network SNR is small. Detector availability
must be **queried from the GWOSC event API at runtime, not hardcoded** — the
papers in `references/` do not tabulate per-event detector participation, and
GWTC-1 Table I notes that its own quoted network SNR sometimes comes from only
two detectors even when three were operating.

---

## 2. Hardware envelope

The execution pool is fixed and small. Every sizing decision below derives from it.

| Resource | Value |
|---|---|
| Nodes | 4 |
| Cores per node | 4 |
| Total cores | **16** |
| Memory per node | 8 GB |
| **Memory per core (fully packed)** | **2 GB** |
| GPU | **none** |

Derived budgets:

| Budget | Target | Rationale |
|---|---|---|
| Per-filter-job memory request | **1.5 GB**, hard ceiling 2 GB | 4 concurrent jobs per node within 8 GB, leaving headroom for the OS and page cache |
| Per-filter-job cores | **1** | maximum job-level parallelism; no threaded FFT contention between the 4 jobs on a node |
| Total matched-filter cost | **<= 40 core-hours** | about 2.5 h wall clock on 16 cores |
| Total wall clock, full event set | **<= 6 h** | overnight-runnable, re-runnable within a working day |
| Peak scratch | **<= 50 GB** | |
| Staged-out products | **<= 2 GB** | |
| Concurrent GWOSC downloads | **<= 4** | GWOSC is a shared public service; rate-limit deliberately |

---

## 3. Constraints

### 3.1 Data constraints

| # | Constraint | Source / rationale |
|---|---|---|
| D1 | Input is GWOSC public strain only. No proprietary frames, no LIGO-internal channels, no pre-computed posterior samples. | project rule, `CLAUDE.md` |
| D2 | Use the **4096 Hz** GWOSC files, not 16384 Hz. | The 4 kHz data is valid to about 1700 Hz (anti-alias roll-off); the analysis band tops out at 512 Hz, so the limit does not bite, and bandwidth and I/O drop 4x. |
| D3 | Use the **4096 s** segment files, not the 32 s snippets. | PSD estimation uses 16 s blocks with an inverse-PSD truncation of 16 s; a 32 s file cannot support this and leaves no off-source stretch. |
| D4 | Decimate to **2048 Hz** before filtering. | The highest-mass event merges near 48 Hz (GW170823) and the lowest-mass near 221 Hz (GW170608); a 1024 Hz Nyquist has ample margin, and halves FFT cost against 4096 Hz. |
| D5 | Analysis band **20-512 Hz**. | 20 Hz matches the GWTC-1 parameter-estimation default (GWTC-1 App. B). GWOSC data is high-passed at 8 Hz and calibrated from 10 Hz, so 20 Hz is safely inside the trustworthy band. |
| D6 | Reject any analysis span containing NaNs, failing CBC_CAT1 or CBC_CAT2, or overlapping a set CBC hardware-injection bit. | GWOSC bitmasks; CAT1 failures are NaN by construction. |
| D7 | Passing all data-quality flags must **not** be treated as proof the span is artifact-free. | Stated explicitly in the GWOSC paper. |
| D8 | Event GPS times must be resolved from the **GWOSC event API**, not from `references/gwtc1_reference.json`. | The papers publish UTC to 0.1 s only; the JSON deliberately omits GPS times rather than asserting derived ones. |

### 3.2 Method constraints

| # | Constraint | Source / rationale |
|---|---|---|
| M1 | Pipeline logic — conditioning, PSD estimation, whitening, the filtering loop, clustering, coincidence, chirp-mass estimation — is implemented from the equations and stated parameters in `references/`, with the paper section or equation cited in a code comment. | `CLAUDE.md` |
| M2 | Waveform generation uses an established library approximant (`SEOBNRv4_ROM`, matching the GWTC-1 PyCBC bank above 4 Msun; `IMRPhenomD` acceptable as a fallback). Re-deriving a waveform model is out of scope. | §4 |
| M3 | **Non-spinning bank**: component spins set to zero. | Fits the compute envelope. Cost: an SNR deficit for events with non-zero `chi_eff`, largest for GW170729 (0.37). This is a named, accepted limitation, not a defect. |
| M4 | PSD from a median-averaged Welch estimate over **16 s** blocks, taken from an off-source stretch that **excludes +/-16 s around the event**, with the inverse PSD truncated to **16 s** in the time domain. | PyCBC paper §3.2-3.3: the 16 s truncation bounds the filter impulse response, so 16 s estimation blocks follow naturally. Excluding the event prevents the signal biasing its own PSD. |
| M5 | Discard **16 s** at each end of every filtered segment as filter-corrupted. | Consequence of M4. |
| M6 | Gate short high-amplitude transients before filtering: zero about 1 s with a **Tukey window, alpha = 0.5**, identified by thresholding the whitened strain magnitude. | PyCBC paper §3.2. A single glitch otherwise produces SNR ~1000 triggers at its own time plus an SNR ~10 skirt over **+/-8 s** from the template and inverse-PSD impulse responses. |
| M7 | Cluster the SNR time series in **1 s** windows, keeping a window maximum only if it exceeds both neighbours. | PyCBC paper §3.3. |
| M8 | Coincidence requires the **same template** in both detectors and arrival times within **15 ms** (10 ms light travel + 5 ms timing allowance). | PyCBC paper §3.5. |
| M9 | Network SNR is the quadrature sum of per-detector SNR over the coincident detectors. | PyCBC paper §3.5, eq. 16. |
| M10 | Chirp mass is reported in the **detector frame**, and compared against GWTC-1 source-frame values only after conversion. See §6.2 — this is the single easiest thing to get wrong in the whole workflow. | |
| M11 | Templates are generated **one at a time inside the filter loop** and not cached across the bank. | A frequency-domain template at 2048 Hz over 512 s is about 8.4 MB as complex128; caching 100 would be 840 MB, over half the per-job memory budget. |
| M12 | No chi-squared reweighting of SNR is required for the primary result. If implemented, the unweighted SNR must also be reported. | GWTC-1 compares against a *reweighted* ranking statistic, but its Table I quotes plain network SNR. Mixing the two would corrupt the §6.1 comparison. |

### 3.3 Template bank constraints

| # | Constraint | Value |
|---|---|---|
| B1 | Detector-frame chirp mass coverage | **5-60 Msun** (spans the required 8.45-52.75 with margin) |
| B2 | Mass ratio coverage | 1:1 to 1:6 |
| B3 | Spins | zero (M3) |
| B4 | Minimal match | **>= 0.97** |
| B5 | Target size | **3,000-6,000 templates**; the realized count must be reported |
| B6 | **Bank adequacy gate** | The bank must contain a template within **1%** in chirp mass of every Tier-A event's expected detector-frame value from §1.2. This is checkable *before* any filtering runs and must gate the search. |
| B7 | Bank is generated **once** per workflow run and shared by all events and detectors | Enables the exact-template coincidence of M8 (PyCBC paper §3.5) |

### 3.4 Workflow and platform constraints

| # | Constraint |
|---|---|
| W1 | Pegasus workflow; all five catalogs generated by `workflow_generator.py`. |
| W2 | Fan-out granularity is `(event, detector, bank_chunk)` for filtering, so the 16-core pool stays saturated. Target roughly 150-250 filter jobs. |
| W3 | Every job requests `cores=1` and `memory <= 1.5 GB` except the bank build and the final report. |
| W4 | Fetch jobs are rate-limited by a DAGMan category with `maxjobs = 4`. |
| W5 | `stage_out=True` only on final products: per-event validation verdicts, the aggregate report, spectrograms, the bank manifest, and the trigger summaries. All conditioned strain, whitened series, and per-chunk trigger files stay in scratch with `stage_out=False`. `register_replica=False` throughout. |
| W6 | Files pass between jobs as explicit `File` objects and explicit path arguments. No directory scanning. |
| W7 | One container image shared by all jobs, defined in `Apptainer/`, based on Ubuntu 24.04 or a python:3.12-era image to match the glibc 2.39 host. |
| W8 | The workflow must run to completion with a single event and a reduced bank as a smoke test, driven by a CLI flag, without editing code. |
| W9 | No GPU, no MPI, no multi-node single-job parallelism. |

---

## 4. Non-constraints — explicitly out of scope

Listing these matters as much as the constraints: several are things a reader
would reasonably expect a "reproduce LIGO" workflow to do, and it cannot.

| # | Out of scope | Why |
|---|---|---|
| N1 | **False alarm rate, inverse FAR, and significance in sigma** | Not reproducible from event snippets, at all. The published FARs come from time-shifting **the full run**: `T_b = T^2 / delta` with `delta = 0.1 s`, requiring about 15 days of coincident data to reach 5 sigma, against 608,000 years of equivalent background for GW150914. A 4096 s segment cannot produce a background estimate. **Do not report a FAR or a sigma value.** |
| N2 | **p_astro** | Requires the four-category Poisson mixture fit over the whole run's trigger distribution (GWTC-1 §VII). |
| N3 | **Bayesian parameter estimation and posteriors** | No `LALInference`, `bilby`, or `RapidPE`. This workflow produces **point estimates from the best-matching template**, not credible intervals. Any comparison against a GWTC-1 90% interval is a containment check, never an interval-vs-interval comparison. |
| N4 | **Component masses m1 and m2 individually** | A matched filter constrains chirp mass tightly and mass ratio weakly; splitting the two needs full PE. Report chirp mass only. |
| N5 | **Spin, `chi_eff`, `chi_p`, precession** | Excluded by the non-spinning bank (M3). |
| N6 | **Final mass, final spin, radiated energy, peak luminosity** | Derived in GWTC-1 from numerical-relativity fits applied to PE posteriors. |
| N7 | **Sky localization and luminosity distance** | Needs coherent multi-detector PE and the antenna response; two detectors localize to hundreds of square degrees at best. |
| N8 | **Source-frame masses as a primary result** | Requires a redshift, which requires a distance, which requires PE (N7). Source frame appears only as a derived cross-check using the **published** redshift, clearly labelled as such. |
| N9 | **Higher-order modes, eccentricity, tidal effects** | |
| N10 | **A blind search for new events** | The workflow is a targeted reanalysis at known times. It cannot discover anything. |
| N11 | **GW170817 and any BNS** | §1.2. |
| N12 | **O3 and later events** | The papers in `references/` cover O1 and O2 only. GWTC-2, 2.1 and 3 are not in scope. |
| N13 | **Low-latency or streaming operation** | Offline, archival, batch. |
| N14 | **GPU acceleration** | No GPU in the pool. |
| N15 | **Reproducing GWTC-1's ranking statistics** | Neither the PyCBC network statistic nor the GstLAL log-likelihood ratio is reimplemented. Consequence: our trigger *ranking* differs from the published one, which feeds directly into the §6.1 tolerance. |
| N16 | **Detector characterization from auxiliary channels** | GWOSC publishes distilled bitmasks, not the hundreds of thousands of auxiliary channels behind them. |

---

## 5. Expected outcomes

### 5.1 Artifacts

| Artifact | Staged out | Content |
|---|---|---|
| `output/bank/manifest.json` | yes | template count, parameter coverage, minimal match, B6 adequacy report |
| `output/<event>/triggers.json` | yes | coincident candidates: per-detector SNR, network SNR, template parameters, times |
| `output/<event>/estimates.json` | yes | best-candidate detector-frame chirp mass, per-detector merger GPS time, inter-detector delay, network SNR |
| `output/<event>/<det>_spectrogram.png` | yes | Q-transform of whitened strain, about +/-1 s around the event |
| `output/<event>/validation.json` | yes | per-criterion pass/fail against GWTC-1 with recovered value, reference value, tolerance, margin |
| `output/report.md` + `report.json` | yes | aggregate table across events, overall verdict, named failures |
| `output/controls/*.json` | yes | off-source null test and injection-recovery results (§6.5) |
| conditioned strain, PSDs, whitened series, per-chunk triggers | no | scratch only (W5) |

### 5.2 Quantitative expectations

- **Tier A: 7 of 7** events recovered as H1-L1 coincidences with network SNR >= 8.
- **Recovered network SNR** within tolerance of GWTC-1 Table I (§6.1). Recovered
  values **may legitimately exceed** the published ones — see §6.1.
- **Detector-frame chirp mass** consistent with GWTC-1 Table III converted to the
  detector frame (§6.2), for all Tier A events.
- **Merger time** within tolerance of the published UTC (§6.3), and
  `|dt_H1L1| <= 10 ms` for every recovered event.
- **Visible chirp** in the spectrogram for events with recovered network SNR > 13:
  expected for GW150914, GW170608, GW170814, GW151226, GW170104. **Not** expected
  for GW170809, GW170823, or any Tier B event — absence there is not a failure.
- **Tier B:** GW151012 and GW170729 recovered at reduced SNR, chirp mass possibly
  outside the tight tolerance. **GW170818 is expected to fail** two-detector
  coincidence, reproducing PyCBC's published non-detection.

### 5.3 Explicit non-expectations

The final report must state, in the report itself and not only here, that this
workflow produces **no statement of statistical significance** (N1), and that
recovering a signal at a known time is not evidence of detection in the sense the
published papers mean. It must also record that GWTC-1's own values supersede the
individual detection papers, and that its GW150914 numbers differ from the 2016
PRL (for instance `d_L = 440 +150 -170` against `410 +160 -180 Mpc`) because of
improved PSD estimation, spline calibration envelopes, and O2 noise subtraction.

---

## 6. Validation criteria

All reference values come from `references/gwtc1_reference.json`. Every criterion
emits a structured record with recovered value, reference value, tolerance, and
signed margin — never a bare boolean.

### 6.1 Network SNR

**Reference:** GWTC-1 Table I, `searches.pycbc.network_snr`, the H1-L1 quadrature
sum.

**Two caveats that set the tolerance, both quoted from the Table I caption:**

> "The network SNR for the two matched-filter searches is that of the template
> ranked highest by that search, which is not necessarily the template with the
> highest SNR."

> "the network SNR is the quadrature sum of the detectors coincident in the
> highest-ranked trigger; in some cases, only two detectors contribute, even if
> all three are operating nominally at the time of that event."

So the published number is the SNR of the **highest-ranked** template under a
reweighted statistic this workflow does not implement (N15, M12), not the maximum
over the bank. **Our maximum-SNR value can legitimately come out higher.** The
criterion is therefore two-sided but asymmetric, and exceeding the reference is
not a failure.

| Criterion | Tier A | Tier B |
|---|---|---|
| V1a Recovered as an H1-L1 coincidence, network SNR >= 8 | required | recorded, not required |
| V1b Recovered network SNR within **-15% / +25%** of reference | required | -30% / +40% |
| V1c Per-detector SNR >= 5.5 in at least one detector | required | recorded |

GW170729 is exempt from V1b even at Tier B tolerance: `chi_eff = 0.37` against a
non-spinning bank (M3), combined with an ISCO of 35 Hz against `f_low = 20 Hz`
(D5), makes the SNR deficit expected and physical. Record the value, flag the
reason, do not fail the run.

### 6.2 Chirp mass — frame conversion is mandatory

**This is the criterion most likely to be silently wrong.** A matched filter
recovers a **detector-frame** chirp mass. GWTC-1 Table III publishes
**source-frame**. They differ by `(1 + z)`, which ranges from 1.07 to 1.49 across
this event set — up to a **49% discrepancy** if the frames are conflated. A naive
comparison would fail every event while looking like a physics bug.

**Procedure:**

1. Take source-frame `chirp_mass_source_msun` and `redshift` from
   `references/gwtc1_reference.json`.
2. Convert to the detector frame: `Mc_det = Mc_source * (1 + z)`.
3. Propagate the 90% interval bounds through the same conversion, using the
   redshift interval bounds, to obtain a detector-frame reference interval.
4. Compare the recovered best-template chirp mass against that interval.

| Criterion | Tier A | Tier B |
|---|---|---|
| V2a Recovered `Mc_det` inside the propagated detector-frame 90% interval | required | recorded |
| V2b Or, if V2a fails, within **10%** of the propagated median | required | within 20% |
| V2c Recovered value not at a bank edge (5 or 60 Msun) | required | required |
| V2d Bank chirp-mass spacing near the recovered value is finer than the V2b tolerance | required | required |

V2c catches a rail-pinned fit masquerading as a measurement. V2d prevents
claiming a 10% agreement that the bank's own resolution cannot support.

The report may additionally quote a source-frame value obtained by dividing by
`(1 + z_published)`, but must label it as **using the published redshift, not
measured** (N8).

### 6.3 Merger time

**Reference:** GWTC-1 Table I `time_utc`, published to 0.1 s, plus `date_utc`.

| Criterion | Tier A | Tier B |
|---|---|---|
| V3a Recovered coalescence time within **+/-100 ms** of the published UTC | required | +/-200 ms |
| V3b `\|dt_H1L1\| <= 10 ms` (inter-site light travel) | required | required |
| V3c For GW150914 only: `dt` consistent with **L1 first, 6.9 +0.5 -0.4 ms** (2016 PRL) | required | n/a |

V3a's tolerance is set by the published precision (0.1 s), not by our timing
resolution, which is much finer. V3b is a strong internal consistency check that
needs no published value. V3c is the only per-event inter-detector delay any of
the four papers publishes, so it is the only place a tight timing check is
possible — use it.

### 6.4 Spectrogram

Qualitative and explicitly non-blocking.

| Criterion | Applies to |
|---|---|
| V4a A monotonically rising time-frequency track is visually identifiable | events with recovered network SNR > 13 |
| V4b Track terminates below the event's `f_ISCO` from §1.2, within a factor of about 2 | same |
| V4c Rendered without error, correct axes, units, and event label | all events |

V4a and V4b **must not gate the run.** A spectrogram is a visualization, and
GWTC-1's own Fig. 10 shows several events with no visible track.

### 6.5 Controls — independent of GWTC-1

These validate the pipeline without reference to the published answers, and are
the only defense against a pipeline that "recovers" a signal because it was told
where to look.

| Criterion | Test |
|---|---|
| V5a **Off-source null test.** Filter 10 spans of equal length, each >= 200 s from any catalog event, with the identical configuration. No coincident candidate may exceed the network SNR recovered at the true event time. | required |
| V5b **Injection recovery.** Inject a known-parameter non-spinning waveform into off-source data at 5 chirp masses spanning the bank. Recover chirp mass within V2b tolerance and time within V3a tolerance. | required |
| V5c **Provenance record.** Every event records the GWOSC channel name, calibration version, GPS span, sample rate, and the exact DQ bit values for the analyzed span. | required |
| V5d **Determinism.** Two runs with identical inputs produce identical chirp mass and merger time, and network SNR agreeing to 1e-6. | required |

V5b is the criterion that distinguishes a working pipeline from a lucky one: it
tests recovery where the true answer is known by construction and no published
value exists to be accidentally fitted toward.

### 6.6 Overall verdict

| Verdict | Condition |
|---|---|
| **PASS** | All V5 controls pass, and all Tier A events pass V1a, V1b, V1c, V2a-or-V2b, V2c, V2d, V3a, V3b, and V3c where applicable |
| **PASS WITH NOTES** | The above holds except for named, explained Tier A deviations, or any Tier B outcome |
| **FAIL** | Any V5 control fails, or any Tier A event fails V1a, V2c, V2d, or V3b |

A V5 failure is a **hard** fail regardless of how well the event comparisons went:
a pipeline that cannot pass its own null and injection tests has not validated
anything, however closely its numbers match the published table.

---

## 7. Workflow structure

### 7.1 Steps

| Step | Instances | Cores | Memory | Notes |
|---|---|---|---|---|
| `fetch_strain` | event x detector | 1 | 1 GB | GWOSC download; DAGMan category `download`, maxjobs 4 (W4) |
| `check_quality` | event x detector | 1 | 1 GB | D6, D7; may veto a detector, emitting a recorded verdict |
| `condition_strain` | event x detector | 1 | 1.5 GB | D4, D5, M4, M5, M6 |
| `bank_psd` | 1 | 1 | 1.5 GB | B7 + Usman eq. 7; fan-in over `condition_strain` |
| `bank_row` | N_q (6) | 1 | 1.5 GB | B1-B5; rows are independent, so they fan out |
| `merge_bank` | 1 | 1 | 2 GB | concatenation; gates the search on B6 |
| `matched_filter` | event x detector x chunk | 1 | 1.5 GB | M2, M3, M7, M11 |
| `coincidence` | event | 1 | 1 GB | M8, M9 |
| `estimate_params` | event | 1 | 1 GB | M10, detector-frame chirp mass and merger time |
| `spectrogram` | event x detector | 1 | 1.5 GB | V4 |
| `validate_event` | event | 1 | 1 GB | §6.1-6.3, consumes `references/gwtc1_reference.json` |
| `control_offsource` | 10 | 1 | 1.5 GB | V5a |
| `control_injection` | 5 | 1 | 1.5 GB | V5b |
| `aggregate_report` | 1 | 1 | 1 GB | fan-in over all events and controls |

Approximate job count for 10 events, 2 detectors, 8 bank chunks: 20 + 20 + 20 +
(1 + 6 + 1) + 160 + 10 + 10 + 20 + 10 + 10 + 5 + 1 = **294 jobs**, satisfying W2.
The generator emits 294 for Tier A (7 events, 14 detector-segments) and 375 for
all tiers.

### 7.2 Dependency shape

```
                     ┌──> bank_row x N_q ──┐
bank_psd ────────────┴─────────────────────┴──> merge_bank ─┬──> matched_filter ──> coincidence ──> estimate_params ──> validate_event ──┐
   ▲                                                        │        (per chunk)       (per event)                                       │
fetch_strain ──> check_quality ──> condition_strain ─────────┤                                                                           ├──> aggregate_report
   (per event, detector)             └──> spectrogram        │                                                                           │
                                                            └──> control_offsource, control_injection ─────────────────────────────────────┘
```

`merge_bank` is a genuine fan-in point for the filter stage and must complete
before any filtering starts, because exact-template coincidence (M8) requires
every detector to have been filtered against the identical bank (B7).

The bank build is split across three job types rather than the one job §7.1
originally specified. Placement at fixed mass ratio never consults another row
(Usman §3.1), so the N_q rows are independent and fan out; measured on real
GW150914 PSDs the default bank is 3,453 templates at **28 min 56 s** serial on
one core, with rows costing 4.6-5.3 min each, so the fan-out replaces 29 min of
serial work in front of 192 filter jobs with about 5.5 min of wall clock. The
placement PSD is its own job so each row stages one ~4 MB spectrum instead of
every conditioned segment (~30 MB x 14) once per row. Both the serial
`build_bank.py` and the parallel path call the same placement code in
`gwlib.py`, and were verified to produce bit-identical template sets, PSDs and
adequacy-gate records on real data.

### 7.3 Parameterization

`workflow_generator.py` accepts at minimum:

- `--events` — event names; default is the Tier A list
- `--detectors` — default `H1 L1`; `V1` opt-in
- `--tier` — `A`, `B`, or `all`
- `--bank-chunks` — filter fan-out width, default 8
- `--mc-range`, `--minimal-match`, `--approximant` — bank parameters (B1, B4, M2)
- `--sample-rate`, `--f-low`, `--f-high` — conditioning parameters (D4, D5)
- `--smoke-test` — one event, coarse bank, satisfying W8
- `--skip-controls`, `--skip-spectrograms` — for iteration only; either flag must
  force the overall verdict (§6.6) to `INCOMPLETE`, never `PASS`

---

## 8. Deliverables

Per the `pegasus-scaffold` skill, the project is complete when all of these exist:

- [ ] `workflow_generator.py` — all five catalogs plus the DAG
- [ ] `bin/` — one wrapper per step in §7.1
- [ ] `Apptainer/GW_Container.def` — Ubuntu 24.04 base (W7)
- [ ] `run_manual.sh` — every step end to end for one event, outside Pegasus
- [ ] `README.md` — usage, event selection, how to read the validation report
- [ ] `SPEC.md`, `CLAUDE.md`, `references/` — this document and its inputs

---

## 9. Open questions

Recorded rather than guessed, per `CLAUDE.md`:

1. **Bank size against the compute budget.** §2 allows 40 core-hours and B5
   targets 3,000-6,000 templates. The per-template filter cost must be **measured**
   on one node before the bank is finalized. If measurement shows headroom, prefer
   spending it on a higher minimal match over a wider mass range — coverage is
   already adequate per B6.
2. **Virgo inclusion.** Whether V1 improves or degrades network SNR at this
   sensitivity is an empirical question. Run H1-L1 first; add V1 as a documented
   variant and report both.
3. **Gating threshold.** The PyCBC paper specifies the Tukey window and its shape
   parameter but not a universal amplitude threshold on the whitened strain. This
   must be tuned against the off-source control (V5a) and the chosen value
   recorded, not silently defaulted.
4. **Off-source span selection.** V5a requires spans >= 200 s from any catalog
   event. Within a single 4096 s file around a known event, verify enough clean,
   CAT-passing span exists; if not, additional GWOSC segments are needed and the
   download budget (W4) grows.

---

## 10. Resolved during implementation

Recorded here rather than left implicit; each was verified against the live
service or a known answer, not assumed.

| Question | Resolution |
|---|---|
| **Waveform approximant (M2).** SEOBNRv4_ROM is what GWTC-1's PyCBC bank used. | **IMRPhenomD is the default.** LALSuite >= 7.25 resolves SEOBNRv4_ROM's interpolation grids (`SEOBNRv4ROM_v3.0.hdf5`) from `$LAL_DATA_PATH`, and those ship in a separate multi-GB data release, not the pip wheel. IMRPhenomD is fully analytic, needs no data files, is a complete inspiral-merger-ringdown aligned-spin model, and is permitted by M2. SEOBNRv4_ROM stays selectable via `--approximant` for anyone with the data; `build_bank` probes the approximant once and fails fast rather than silently degrading every template to the PN fallback. Smoke testing on a host without PyCBC requires the explicit `--allow-fallback-waveform` opt-in, and such a bank can still never yield a `PASS`. |
| **GWOSC event API endpoint (D8).** | The catalog segment is **required**: `/eventapi/json/{catalog}/{event}/`. A bare `/eventapi/json/{event}/` returns 404. Default catalog `GWTC-1-confident`, with fallback to enumerating the catalogs the root advertises. |
| **File alignment.** | GWOSC 4096 s event files are **event-centred**, not GPS-aligned: GW150914 sits 2047 s into its file with 2049 s after. `condition_strain` reads the actual `GPSstart` so either layout works, and `fetch_strain` now reports the margins. |
| **§9.3 gate threshold.** | 50 (whitened strain) gates 2 transients per 256 s analysis segment on real GW150914 data and leaves the off-source null test empty. Retained as the default; still a tuning choice, not a published value. |
| **§9.1 bank size vs compute budget.** | Measured, not estimated: 49 templates over ℳ_det 25–40 at minimal match 0.97 build in **26 s**, and filtering a 256 s segment against them takes **5 s per detector** (IMRPhenomD, one core). The 40 core-hour budget is not the binding constraint at this scale; a full 5–60 bank should be re-measured before committing. |

### Bugs found by the verification steps

None of these raised an exception; each silently produced plausible wrong
numbers, which is why the invariants in `run_manual.sh` step 0 exist.

| Bug | Detection | Consequence if unfixed |
|---|---|---|
| Matched-filter SNR built from two real `irfft` calls | zero-noise self-filter must peak at `sqrt(<h\|h>)`; measured 2.000000 | every SNR exactly 2x too large |
| Inverse-PSD truncation applied to `1/S` instead of the amplitude `1/sqrt(S)` | truncation must preserve the in-band PSD level; measured 0.0015 | in-band PSD ~700x low, so every SNR wildly inflated |
| No high-pass before filtering | synthetic segment gave network SNR ~5.7e6 | sub-band power leaks across the segment FFT and swamps the band |
| Crop shorter than 2x the edge guard | V5b reported "injection not recovered" | every trigger discarded, indistinguishable from a real miss |
| Analysis window required to be **event-centred** | first Pegasus run: `condition_strain` exited 1 on GW151226/H1, holding the job and stalling 274 downstream jobs | any event whose detector came online shortly before it is unanalysable |

### The event-centred-window flaw (found by the first real Pegasus run)

Worth separating out, because the data was fine and the code was wrong.

GW151226's H1 detector comes online **2025 s into its 4096 s file**: only 53% is
science mode, in a single clean run of 2171 s, with the event just **123 s from
that run's start** — and 2048 s of clean data after it. A fixed +/-128 s
*event-centred* window missed by five seconds and the DQ gate correctly rejected
it (D6), which then hard-failed the job.

Real pipelines analyse **contiguous science segments**, not fixed windows around
a time. So:

* `check_quality` now computes the contiguous clean runs (DATA & CBC_CAT1 &
  CBC_CAT2), finds the run containing the event, and **places** the analysis
  window inside it — preferring centred, shifting when the event sits near a
  boundary, shortening only if the run cannot host the full length. It reports
  `analysis_span`, `analysis_shifted`, `analysis_truncated` and the run list.
* `condition_strain` takes that span as given rather than recomputing it, and
  restricts the PSD stretch and every off-source span to clean runs. A span
  straddling a science-mode gap is NaN-filled; zero-filled it looks like a giant
  glitch, which would poison the PSD and fake a loud off-source trigger —
  defeating the very null test it feeds.
* `bank_psd` interpolates PSDs onto the finest common grid instead of requiring
  equal lengths, since analysis segments may now legitimately differ. (The
  serial `build_bank` does the same.)

GW151226/H1 needed a shift of only **5.4 s** and kept the full 256 s. 13 of the
14 Tier A detector-segments were unaffected.
