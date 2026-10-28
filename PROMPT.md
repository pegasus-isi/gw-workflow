# Starting prompts: gravitational-wave detection reanalysis

Prepared 2026-08-19. Take these into PegasusAI Studio to build the SPEC.
Structured as a spec-first prompt pair, the same shape as the EHT capstone's
opening prompts, so no code is requested until the specification and the
source rules exist as reviewable text.

## Before prompting

1. Create `references/` in the workflow project and download the five PDFs
   below into it.
2. Have the pool's shape at hand (the prompt states it): 4 nodes, 4 cores
   and 8 GB each, no GPU.
3. Nothing else. GWOSC needs no credentials, and the strain files are a few
   megabytes each (32 s at 16 kHz per detector).

## Papers to put in `references/`

| Role | Paper | DOI |
|---|---|---|
| **Primary methods reference** | A guide to LIGO-Virgo detector noise and extraction of transient gravitational-wave signals, CQG 37, 055002 (2020) | `10.1088/1361-6382/ab685e` |
| The landmark result | Observation of Gravitational Waves from a Binary Black Hole Merger, PRL 116, 061102 (2016) | `10.1103/PhysRevLett.116.061102` |
| **Validation targets** | GWTC-1: A Gravitational-Wave Transient Catalog of Compact Binary Mergers, PRX 9, 031040 (2019) | `10.1103/PhysRevX.9.031040` |
| The data release | Open data from the first and second observing runs of Advanced LIGO and Advanced Virgo, SoftwareX 13, 100658 (2021) | `10.1016/j.softx.2021.100658` |
| Search methodology | The PyCBC search for gravitational waves from compact binary coalescence, CQG 33, 215004 (2016) | `10.1088/0264-9381/33/21/215004` |

All five DOIs were verified against Crossref on 2026-08-17.

The CQG 2020 guide is the anchor, not the discovery paper: the collaboration
wrote it to explain how transient signals are extracted from the open strain
data, so it plays the role EHT Paper IV played for the capstone. GWTC-1 is
what makes this a reference comparison rather than an executed-only run,
because it publishes the per-event values to check against.

## Prompt 1, orientation

```
can you read the papers in references and provide me a summary? focus on
how a transient signal is extracted from open strain data and which
published values are quoted per event.
```

## Prompt 2, specification and source rules

```
I want to reproduce LIGO's binary black hole detections using only the
publicly released GWOSC strain data and open tools, not anyone's tutorial
pipeline. Help me build a SPEC.md for a Pegasus workflow that fetches
strain data for a set of GWTC-1 events, conditions and whitens it per
detector, runs a matched filter against a template bank, combines
detectors into a network SNR, estimates chirp mass and merger time, and
renders the time-frequency spectrogram. Let define constraints,
non-constraints, expected outcomes, and validation criteria, where the
validation criteria compare our recovered SNR, chirp mass and merger time
against the values published in GWTC-1 for each event. Also note that the
pool has 4 nodes of 4 cores and 8GB with no GPU, so the template bank and
event count need to fit that. Please create a CLAUDE.md with an explicit
rule to use only GWOSC data, the published papers, and official library
documentation, and not to copy the GWOSC or PyCBC tutorial notebooks.
```

## Why the prompt is worded this way

**The clean-room clause.** Well-known "find GW150914" tutorial notebooks
hand over the entire pipeline. If the agent leans on one, the run measures
copying rather than construction and stops being independent evidence.
Allowing each library's own API documentation while excluding the worked
notebooks is the line that keeps the claim strong without blocking the
build. This mirrors the EHT capstone's rule against the prior
reproduction's artifacts.

**The hardware sentence.** The natural failure mode in this domain is an
agent proposing a full template-bank search that would never finish on four
nodes. Stating the envelope up front turns that into a scoping decision in
the specification, where it is cheap, instead of a runtime discovery.

**No code in either prompt.** The point of the pattern is that review
happens on plain text. Read the SPEC and the source rules, redirect them,
and only then ask for an implementation.

## Validation targets for the first event

Pulled live from the GWOSC event API on 2026-08-17 for GW150914, as an
example of what the workflow's comparison job checks against:

| Quantity | Published value |
|---|---|
| GPS merger time | 1126259462.4 |
| Network matched-filter SNR | 25.2 |
| Chirp mass (source frame) | 28.6 solar masses |
| Component masses | 35.6 and 30.6 solar masses |
| Luminosity distance | 440 Mpc |

The API also lists 24 strain files per event across detectors and formats
(`gwf`, `hdf5`, `txt` at 4 and 16 kHz), so the fetch job has a documented
choice to make.

## Evidence to capture during the build

If this becomes a case study, these are the records the paper needs, so
capture them as you go rather than reconstructing later.

- Every prompt, in order, including the ones that went nowhere.
- Which lifecycle stage each interaction served: generate, plan, submit,
  monitor, diagnose, or fix.
- Every diagnose-fix-resubmit round and what triggered it.
- `pegasus-statistics` output for the run of record, plus job counts and
  wall time.
- Container image digests and the input checksums.
- The comparison table of recovered against published values.
