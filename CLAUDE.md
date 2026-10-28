# CLAUDE.md — gw-workflow

Project rules for this workflow. These are binding and take precedence over
general habits or defaults.

**What this project is:** an independent reproduction of the GWTC-1 binary black
hole detections from public GWOSC strain data. See [`SPEC.md`](SPEC.md) for
scope, constraints, and validation criteria.

---

## 1. Sourcing rule — the point of this project

The scientific value of this workflow is that it is derived from primary sources.
A pipeline assembled from tutorial code proves nothing about whether the published
method can be independently reimplemented, which is the entire question being
asked. Copying a tutorial does not just violate a preference here — it destroys
the result.

### 1.1 Permitted sources

1. **GWOSC public data** — strain files, data-quality and injection bitmasks,
   segment lists, and the event API at `https://gw-openscience.org`.
2. **The four papers in [`references/`](references/)**, and the two derived
   documents in that directory (`SUMMARY.md`, `gwtc1_reference.json`) which were
   produced from them:
   - Abbott et al. 2016, *PRL* **116**, 061102 — GW150914
   - Usman et al. 2016, *CQG* **33**, 215004 — the PyCBC search pipeline
   - Abbott et al. 2019, *PRX* **9**, 031040 — GWTC-1
   - Abbott et al. 2021, *SoftwareX* **13**, 100658 — GWOSC open data
3. **Official API documentation** for libraries actually used: PyCBC, LALSuite /
   `lalsimulation`, GWpy, `numpy`, `scipy`, `h5py`, `matplotlib`, Pegasus.
   Reference documentation — function signatures, parameter meanings, return
   types — is permitted.
4. **Pegasus documentation** via the `pegasus-docs` MCP server, and the Pegasus
   skills in this environment.

### 1.2 Forbidden sources

Do not read, copy, adapt, paraphrase, or reconstruct from memory:

- **GWOSC tutorial notebooks**, quickstart guides, or worked examples.
- **GW Open Data Workshop** material — lectures, notebooks, exercises, solutions.
- **PyCBC tutorials and example notebooks**, including anything under an
  `examples/` or `tutorials/` path in the PyCBC distribution or its docs site.
- **GWpy example galleries.**
- Third-party blog posts, Medium articles, Stack Overflow answers, course
  materials, or GitHub repositories that implement a gravitational-wave search or
  a GW150914 reanalysis.
- Any "reproduce GW150914" walkthrough, in any language, from any source.

The line is: **libraries are used as libraries, through their documented API.
Pipelines are written from the papers.**

Calling `pycbc.filter.matched_filter` is permitted — it is a documented library
function. Reproducing the structure of the notebook that demonstrates it is not.
Generating a waveform with `SEOBNRv4_ROM` is permitted; re-deriving the waveform
model is out of scope (`SPEC.md` §4, N-list).

### 1.3 What this requires in the code

- **Cite the source in a comment** for every non-obvious numeric constant and
  every algorithmic choice, naming the paper and section or equation. For example:
  `# Tukey window alpha = 0.5, 1 s zeroed: Usman 2016 §3.2`.
- A constant with no citation and no derivation is a defect. If a value is a
  tuning choice rather than a published one, say so and say what it was tuned
  against.
- Prefer transcribing an equation from a paper over calling an opaque helper that
  happens to do the same thing, **where the equation is the thing being
  reproduced** — conditioning, PSD handling, clustering, coincidence, chirp-mass
  estimation. Do not reimplement FFTs, waveform models, or interpolation.

### 1.4 When a needed detail is not in the permitted sources

Some implementation details are genuinely absent from the papers. `SPEC.md` §9
already records four. When another surfaces:

1. Say plainly that the papers do not specify it.
2. Propose a value with the reasoning that leads to it.
3. State how it will be validated — usually against the off-source or injection
   controls (`SPEC.md` §6.5).
4. Record it in `SPEC.md` §9.

**Do not** resolve such a gap by reaching for a tutorial. An unspecified,
honestly-flagged, control-validated choice is a better outcome than a
tutorial-derived one, even if the tutorial's value is correct.

---

## 2. Result integrity

- **Never report a false alarm rate, an inverse FAR, or a significance in sigma.**
  These are not computable from event snippets — `SPEC.md` §4, N1 explains why in
  detail. If a number resembling one appears in output, it is a bug.
- **Never present a point estimate as a credible interval.** This workflow
  produces best-template point estimates, not posteriors (N3).
- **Chirp mass is detector-frame** unless explicitly converted, and every
  comparison against GWTC-1 must convert (`SPEC.md` §6.2). Source-frame values
  quoted anywhere must be labelled as using the *published* redshift, not a
  measured one.
- **Report failures as failures.** If a Tier A event does not validate, the report
  says so with the recovered and reference values side by side. Do not widen a
  tolerance to turn a failure into a pass; if a tolerance is genuinely wrong,
  change it in `SPEC.md` with the reasoning, as a visible decision.
- **The controls are not optional.** A run with `--skip-controls` reports
  `INCOMPLETE`, never `PASS` (`SPEC.md` §7.3).
- Do not tune any threshold, tolerance, or bank parameter to make a specific
  published value come out right. The off-source and injection controls exist to
  make that kind of tuning detectable.

---

## 3. Hardware budget

The pool is **4 nodes x 4 cores x 8 GB, no GPU**. Every job:

- requests `cores=1` and `memory <= 1.5 GB` (bank build and report may take 2 GB)
- assumes **no GPU** — never emit CUDA, `--nv`, or GPU profiles
- generates templates one at a time and does not cache the bank in memory
  (`SPEC.md` M11)

Before enlarging the bank, the event set, or the sample rate, **measure** the
per-template filter cost on one node and check it against the 40 core-hour budget
(`SPEC.md` §2, §9.1). Do not scale up on an estimate.

---

## 4. GWOSC etiquette

- Concurrent downloads capped at **4** via a DAGMan category (`SPEC.md` W4).
- Cache fetched segments in scratch; never re-download within a run.
- Set an explicit timeout on every HTTP request and check the status code.
- Resolve event GPS times from the **GWOSC event API**, not from
  `references/gwtc1_reference.json`, which deliberately omits them
  (`SPEC.md` D8).
- `README.md` must carry the GWOSC data acknowledgement, as the SoftwareX paper
  requests of anyone publishing on these data.

---

## 5. Workflow conventions

- Everything lives under `~/work/workflows/gw-workflow/`. Nothing in `/tmp`.
- Files pass between jobs as shared `File` objects with explicit path arguments.
  No `glob()`, `os.listdir()`, or directory scanning between jobs.
- `stage_out=True` only on final products; `register_replica=False` throughout
  (`SPEC.md` W5).
- Support files — not the step wrappers themselves — go in the Replica Catalog and
  are located with `os.getcwd()`, never a `__file__`-relative path. This includes
  `references/gwtc1_reference.json`, which the validation job consumes as an input
  file.
- Container base must be Ubuntu 24.04 or a python:3.12-era image, matching the
  glibc 2.39 host.
- Every wrapper: `argparse`, structured logging of inputs and the command run,
  `os.makedirs(..., exist_ok=True)` before writing, and propagation of the
  underlying exit code.

---

## 6. Deliverables

The project is not done until all of these exist on disk and `ls -R` confirms it:

`workflow_generator.py`, `bin/<one wrapper per step>`,
`Apptainer/GW_Container.def`, `run_manual.sh`, `README.md`.

`run_manual.sh` and `README.md` are the two that get dropped when the generator
took a while to write. They are not optional.
