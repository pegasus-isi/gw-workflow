#!/bin/bash
#
# Manual test for the gw-workflow: run every pipeline step locally, in order,
# before attempting a Pegasus submission. Catches argument mismatches, missing
# dependencies and broken file handoffs in seconds instead of a plan+submit
# cycle.
#
# Usage:
#   ./run_manual.sh                 # fetch one real event from GWOSC
#   ./run_manual.sh --offline       # synthesize a GWOSC-shaped file, no network
#   ./run_manual.sh --event GW151226
#   ./run_manual.sh --keep          # keep test_output/ afterwards
#
# --offline exists because the host may have no network, and because the
# synthetic file lets the whole chain be exercised deterministically. Offline
# results carry no physics meaning; the run says so.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${SCRIPT_DIR}/bin"
OUT="${SCRIPT_DIR}/test_output"
REF="${SCRIPT_DIR}/references/gwtc1_reference.json"

EVENT="GW150914"
OFFLINE=false
KEEP=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --offline) OFFLINE=true; shift ;;
        --event)   EVENT="$2"; shift 2 ;;
        --keep)    KEEP=true; shift ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

GREEN='\033[0;32m'; RED='\033[0;31m'; BLUE='\033[0;34m'; YEL='\033[1;33m'; NC='\033[0m'
step() { echo -e "\n${BLUE}=== $* ===${NC}"; }
ok()   { echo -e "${GREEN}  OK${NC} $*"; }
warn() { echo -e "${YEL}  WARN${NC} $*"; }
die()  { echo -e "${RED}  FAIL${NC} $*"; exit 1; }
have() { [[ -s "$1" ]] || die "expected output missing or empty: $1"; }

rm -rf "$OUT"; mkdir -p "$OUT"
cd "$OUT"
# gwlib.py and matched_filter.py are support files: Pegasus stages them into the
# job's working directory, so mirror that here rather than importing from bin/.
cp "$BIN/gwlib.py" "$BIN/matched_filter.py" .
export PYTHONPATH="$OUT:${PYTHONPATH:-}"

# ---------------------------------------------------------------- step 0
step "Step 0: gwlib self-test (signal-processing core)"
python3 - <<'PYEOF'
import sys, numpy as np, gwlib
fail = []
fs, T = 2048, 64
n, dt, df = fs*T, 1.0/fs, 1.0/T
f = gwlib.freqs(n, dt)
psd = np.full(len(f), 1e-46); psd[0] = 1e30
mask = gwlib.band_mask(f, 20.0, 512.0)
h, fb = gwlib.fd_waveform(30.0, 30.0, 20.0, df, 512.0)
h = np.pad(h, (0, len(f)-len(h)))[:len(f)]

# The normalization check that caught a factor-of-2 bug: filtering a template
# against itself in zero noise must peak at exactly sqrt(<h|h>).
opt = np.sqrt(gwlib.sigma_sq(h, psd, df, mask))
got = gwlib.matched_filter_snr(h, h, psd, df, n, mask).max()
r = got/opt
print(f"  zero-noise self-filter ratio = {r:.6f}")
if abs(r-1) > 1e-6: fail.append(f"SNR normalization off by {r:.4f}x")

o = gwlib.overlap(h, h, psd, df, mask)
print(f"  overlap(h,h)                 = {o:.8f}")
if abs(o-1) > 1e-9: fail.append(f"overlap(h,h) = {o}")

ts = np.zeros(n); ts[50000] = 20.0
tr = gwlib.cluster_triggers(ts, fs, 5.5)
print(f"  clustering one spike         -> {len(tr)} trigger(s)")
if len(tr) != 1: fail.append(f"clustering gave {len(tr)} triggers for one spike")

# Inverse-PSD truncation must preserve the in-band level (ratio ~1).
psd_c = 1e-46*(1.0 + (30.0/np.clip(f,15,None))**4 + (f/300.0)**2)
pt = gwlib.truncate_inverse_psd(psd_c, fs, n)
ok, ratio = gwlib.check_truncation(psd_c, pt, f, 20.0, 512.0)
print(f"  inverse-PSD truncation ratio = {ratio:.4f}")
if not ok: fail.append(f"inverse-PSD truncation changed the in-band level {ratio:.4g}x")

mc = gwlib.chirp_mass(*gwlib.mass1_mass2(28.6, 1.163))
print(f"  mass roundtrip Mc            = {mc:.4f} (want 28.6)")
if abs(mc-28.6) > 1e-3: fail.append(f"mass roundtrip gave {mc}")

print(f"  f_isco(72.2 Msun)            = {gwlib.f_isco(72.2):.1f} Hz (GWTC-1: ~61)")
if fb: print("  NOTE: PyCBC absent, restricted-PN fallback in use")
if fail:
    print("\n".join("  FAIL: "+m for m in fail)); sys.exit(1)
PYEOF
ok "core verified"

# ---------------------------------------------------------------- step 1
step "Step 1: fetch_strain"
STRAIN="${EVENT}_H1_strain.hdf5"; META="${EVENT}_H1_fetch.json"
if $OFFLINE; then
    warn "offline: synthesizing GWOSC-shaped files with a known injection"
    warn "         (synthetic noise -- verifies the plumbing, not the physics)"
    python3 - "$EVENT" <<'INNER_PY'
"""Synthesize two GWOSC-shaped 4096 s files with one coherent injection.

The injection is what makes this a functional test: with no signal in the data
the on-source path finds nothing, the null test is vacuous, and the PASS branch
of validate_event is never exercised. Parameters are GW150914's published
detector-frame values so the validation logic meets the numbers it will really
see. The run is still forced to INVALID_FALLBACK_WAVEFORM when PyCBC is absent.
"""
import json, sys, numpy as np, h5py
import gwlib

event = sys.argv[1]
fs, dur = 4096, 4096
gps = 1126259462.4                        # GW150914, GWOSC event API
t0 = float((int(gps) // dur) * dur)       # GWOSC files are GPS-aligned
n = fs * dur

MC_DET = 31.174        # GWTC-1: Mc_source 28.6 x (1+z = 1.09)
Q_INJ = 1.16           # from the published m1/m2
TARGET_SNR = 17.0      # per detector -> network ~24 vs Table I's 23.6
# PRL 116 061102: GW150914 reached L1 first, 6.9 ms before H1.
ARRIVAL = {"H1": gps, "L1": gps - 0.0069}
SEED = {"H1": 11, "L1": 29}
f_low, f_high = 20.0, 512.0

m1, m2 = gwlib.mass1_mass2(MC_DET, Q_INJ)
print("injecting Mc_det=%.3f q=%.2f (m1=%.2f m2=%.2f)" % (MC_DET, Q_INJ, m1, m2))

for det in ("H1", "L1"):
    rng = np.random.default_rng(SEED[det])
    f = np.fft.rfftfreq(n, 1.0 / fs); f[0] = f[1]
    # Flatten below 15 Hz: real GWOSC strain is high-passed at 8 Hz
    # (SoftwareX 2.2), so it carries no runaway sub-band tail.
    fk = np.clip(f, 15.0, None)
    asd = 1e-22 * ((20.0 / fk) ** 4 + 1.0 + (fk / 300.0) ** 2)
    sf = (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f))) * asd
    strain = np.fft.irfft(sf, n=n)

    # Measure the realized PSD, then scale the injection against it so the
    # recovered SNR lands where intended.
    fp, pxx, _ = gwlib.estimate_psd(strain[: 256 * fs], fs)
    seg_s = 64
    ns = seg_s * fs
    dfs = 1.0 / seg_s
    fseg = gwlib.freqs(ns, 1.0 / fs)
    psd_seg = gwlib.interpolate_psd(fp, pxx, fseg)
    mask = gwlib.band_mask(fseg, f_low, f_high)
    h, fb = gwlib.fd_waveform(m1, m2, f_low, dfs, f_high)
    h = np.pad(h, (0, len(fseg) - len(h)))[: len(fseg)]
    sigma = np.sqrt(gwlib.sigma_sq(h, psd_seg, dfs, mask))
    h = h * (TARGET_SNR / sigma)
    h_td = np.roll(np.fft.irfft(h, n=ns) * fs, ns // 2)   # merger at centre

    off = int(round((ARRIVAL[det] - t0) * fs)) - ns // 2
    strain[off : off + ns] += h_td
    print("  %s: injected at GPS %.4f, fallback=%s" % (det, ARRIVAL[det], fb))

    out = "%s_%s_strain.hdf5" % (event, det)
    with h5py.File(out, "w") as hf:
        hf.create_dataset("strain/Strain", data=strain, compression="gzip")
        hf["meta/GPSstart"] = t0
        hf["meta/Duration"] = float(dur)
        hf["meta/Detector"] = np.bytes_(det)
        q = hf.create_group("quality/simple")
        q.create_dataset("DQmask", data=np.full(dur, 0xFF, dtype=np.uint16))
        q.create_dataset("DQShortnames", data=np.array(
            [b"DATA", b"CBC_CAT1", b"CBC_CAT2", b"CBC_CAT3",
             b"BURST_CAT1", b"BURST_CAT2", b"BURST_CAT3"]))
        i = hf.create_group("quality/injections")
        i.create_dataset("Injmask", data=np.full(dur, 0x1F, dtype=np.uint16))
        i.create_dataset("InjShortnames", data=np.array(
            [b"NO_CBC_HW_INJ", b"NO_BURST_HW_INJ", b"NO_DETCHAR_HW_INJ",
             b"NO_CW_HW_INJ", b"NO_STOCH_HW_INJ"]))
    json.dump({"event": event, "detector": det, "event_gps": gps,
               "url": "SYNTHETIC-OFFLINE", "gps_start": t0, "duration": dur,
               "sample_rate": fs, "bytes": 0,
               "injected": {"chirp_mass_detector": MC_DET, "q": Q_INJ,
                            "optimal_snr": TARGET_SNR,
                            "arrival_gps": ARRIVAL[det]},
               "gps_source": "synthetic (run_manual.sh --offline)"},
              open("%s_%s_fetch.json" % (event, det), "w"), indent=2)
    print("  wrote %s" % out)
INNER_PY
else
    for det in H1 L1; do
        python3 "$BIN/fetch_strain.py" --event "$EVENT" --detector "$det" \
            --output "${EVENT}_${det}_strain.hdf5" \
            --meta-output "${EVENT}_${det}_fetch.json" \
            || die "fetch_strain failed for $det (no network? try --offline)"
    done
fi
have "$STRAIN"; have "$META"; ok "both detectors ready"

# ---------------------------------------------------------------- step 2
step "Step 2: check_quality"
python3 "$BIN/check_quality.py" --strain "$STRAIN" --meta "$META" \
    --analysis-length 64 --output "${EVENT}_H1_quality.json"
have "${EVENT}_H1_quality.json"
python3 -c "import json;d=json.load(open('${EVENT}_H1_quality.json'));print('  pass =',d['pass'],'| failures =',d['failures'])"
ok "quality verdict written"

# ---------------------------------------------------------------- step 3
step "Step 3: condition_strain"
python3 "$BIN/condition_strain.py" --strain "$STRAIN" --meta "$META" \
    --quality "${EVENT}_H1_quality.json" \
    --analysis-length 64 --psd-length 256 --offsource-length 64 \
    --n-offsource 3 \
    --output "${EVENT}_H1_conditioned.hdf5" --report "${EVENT}_H1_conditioning.json"
have "${EVENT}_H1_conditioned.hdf5"; ok "conditioned segment + PSD + off-source spans"

ok "H1 conditioned"

# L1 is a genuine second realization (fixture) or a real fetch -- never a copy.
python3 "$BIN/check_quality.py" --strain "${EVENT}_L1_strain.hdf5" \
    --meta "${EVENT}_L1_fetch.json" --analysis-length 64 \
    --output "${EVENT}_L1_quality.json" >/dev/null
python3 "$BIN/condition_strain.py" --strain "${EVENT}_L1_strain.hdf5" \
    --meta "${EVENT}_L1_fetch.json" --quality "${EVENT}_L1_quality.json" \
    --analysis-length 64 --psd-length 256 --offsource-length 64 --n-offsource 3 \
    --output "${EVENT}_L1_conditioned.hdf5" \
    --report "${EVENT}_L1_conditioning.json" >/dev/null
ok "second detector conditioned"

# ---------------------------------------------------------------- step 4
# Exercises the SAME three-stage path the workflow runs (bank_psd -> bank_row
# per q -> merge_bank), not the serial build_bank.py, so a mismatch between a
# wrapper's arguments and the generator shows up here rather than in the DAG.
step "Step 4a: bank_psd (harmonic-mean placement PSD)"
python3 "$BIN/bank_psd.py" \
    --psd "${EVENT}_H1_conditioned.hdf5" --psd "${EVENT}_L1_conditioned.hdf5" \
    --output bank_psd.hdf5 --report bank_psd.json >/dev/null \
    || die "bank_psd failed"
have bank_psd.hdf5
python3 -c "
import json;r=json.load(open('bank_psd.json'))
print('  segments =',r['n_segments'],'| bins =',r['n_bins'],'| df =',round(r['df_hz'],6))"
ok "placement PSD"

step "Step 4b: bank_row x 2 (rows are independent -- this is the fan-out)"
N_Q=2; Q_MAX=4
for i in $(seq 0 $((N_Q - 1))); do
    q=$(python3 -c "import math;n=$N_Q;print(1.0 if n<2 else math.exp(math.log($Q_MAX)*$i/(n-1)))")
    python3 "$BIN/bank_row.py" --psd bank_psd.hdf5 --q "$q" --row "$i" \
        --mc-min 25 --mc-max 40 --minimal-match 0.90 \
        --allow-fallback-waveform \
        --output "bank_row_0${i}.json" >/dev/null \
        || die "bank_row $i failed"
    have "bank_row_0${i}.json"
done
python3 -c "
import glob,json
for f in sorted(glob.glob('bank_row_*.json')):
    r=json.load(open(f))
    print(f'  row {r[\"row\"]} q={r[\"q\"]:.3f}: {r[\"n_templates\"]} templates in {r[\"build_seconds\"]}s')"
ok "rows placed"

step "Step 4c: merge_bank (concatenate + B6 adequacy gate)"
python3 "$BIN/merge_bank.py" \
    $(for i in $(seq 0 $((N_Q - 1))); do printf -- "--row bank_row_0%d.json " "$i"; done) \
    --psd bank_psd.hdf5 --psd-report bank_psd.json \
    --reference "$REF" --tier-a "$EVENT" --adequacy-tol 0.5 \
    --output bank.hdf5 --manifest bank_manifest.json \
    || die "merge_bank failed"
have bank.hdf5
python3 -c "
import json;m=json.load(open('bank_manifest.json'))
print('  templates =',m['n_templates'],'| measured min match =',m['minimal_match_measured_min'])
print('  fallback waveform =',m['waveform_fallback_used'])
print('  cpu total =',m['build_seconds_total_cpu'],'s | slowest row =',m['build_seconds_slowest_row'],'s')"
ok "bank merged"

# ---------------------------------------------------------------- step 5
step "Step 5: matched_filter (both detectors)"
for det in H1 L1; do
    python3 "$BIN/matched_filter.py" \
        --conditioned "${EVENT}_${det}_conditioned.hdf5" --bank bank.hdf5 \
        --chunk 0 --n-chunks 1 --snr-threshold 5.5 \
        --output "${EVENT}_${det}_triggers_000.json" >/dev/null
    have "${EVENT}_${det}_triggers_000.json"
done
python3 -c "
import json
for d in ('H1','L1'):
    b=json.load(open('${EVENT}_'+d+'_triggers_000.json'))
    print(f'  {d}: {b[\"n_triggers\"]} triggers')"
ok "on-source filtering"

step "Step 5b: matched_filter on an off-source span (--dataset, SPEC V5a path)"
for det in H1 L1; do
    python3 "$BIN/matched_filter.py" \
        --conditioned "${EVENT}_${det}_conditioned.hdf5" --bank bank.hdf5 \
        --dataset offsource/span_00 --chunk 0 --n-chunks 1 \
        --output "off_${det}_000.json" >/dev/null
done
ok "off-source filtering through the identical code path"

# ---------------------------------------------------------------- step 6
step "Step 6: coincidence"
python3 "$BIN/coincidence.py" --event "$EVENT" \
    --triggers "${EVENT}_H1_triggers_000.json" \
    --triggers "${EVENT}_L1_triggers_000.json" \
    --output "${EVENT}_coincidence.json"
have "${EVENT}_coincidence.json"
python3 "$BIN/coincidence.py" --event "${EVENT}_offsource_span00" \
    --triggers off_H1_000.json --triggers off_L1_000.json \
    --output off_s00_coincidence.json >/dev/null
ok "coincidence, on- and off-source"

# ---------------------------------------------------------------- step 7
step "Step 7: estimate_params"
python3 "$BIN/estimate_params.py" --coincidence "${EVENT}_coincidence.json" \
    --reference "$REF" --event "$EVENT" --output "${EVENT}_estimates.json"
have "${EVENT}_estimates.json"
ok "point estimates (no FAR, no intervals — SPEC N1/N3)"

# ---------------------------------------------------------------- step 8
step "Step 8: spectrogram"
python3 "$BIN/spectrogram.py" --conditioned "${EVENT}_H1_conditioned.hdf5" \
    --output "${EVENT}_H1_spectrogram.png" --report "${EVENT}_H1_spectrogram.json"
have "${EVENT}_H1_spectrogram.png"; ok "PNG rendered headless"

# ---------------------------------------------------------------- step 9
step "Step 9: validate_event"
python3 "$BIN/validate_event.py" --estimates "${EVENT}_estimates.json" \
    --coincidence "${EVENT}_coincidence.json" --bank-manifest bank_manifest.json \
    --reference "$REF" --event "$EVENT" --tier A \
    --output "${EVENT}_validation.json"
have "${EVENT}_validation.json"
python3 -c "
import json;v=json.load(open('${EVENT}_validation.json'))
print('  verdict =',v['verdict'])
for c in v['criteria']: print(f\"    {'PASS' if c['pass'] else 'FAIL'}  {c['criterion']}\")"
ok "validation criteria evaluated"

# ---------------------------------------------------------------- step 10
step "Step 10: controls"
python3 "$BIN/control_offsource.py" --offsource off_s00_coincidence.json \
    --onsource "${EVENT}_coincidence.json" --event "$EVENT" \
    --output control_offsource.json
have control_offsource.json
python3 "$BIN/control_injection.py" --conditioned "${EVENT}_H1_conditioned.hdf5" \
    --bank bank.hdf5 --span offsource/span_00 --chirp-mass 31.0 \
    --target-snr 20 --crop-length 64 --output control_injection_00.json
have control_injection_00.json
python3 -c "
import json
for f in ('control_offsource.json','control_injection_00.json'):
    d=json.load(open(f)); print(f\"  {d['criterion']}: pass={d['pass']} {d.get('failure_reason','')}\")"
ok "both hard gates ran"

# ---------------------------------------------------------------- step 11
step "Step 11: aggregate_report"
python3 "$BIN/aggregate_report.py" \
    --validation "${EVENT}_validation.json" \
    --control control_offsource.json --control control_injection_00.json \
    --conditioning "${EVENT}_H1_conditioning.json" \
    --spectrogram "${EVENT}_H1_spectrogram.json" \
    --bank-manifest bank_manifest.json \
    --output-md report.md --output-json report.json
have report.md; have report.json
echo; sed -n '1,8p' report.md | sed 's/^/  /'
ok "report generated"

# ---------------------------------------------------------------- summary
step "Summary"
ls -la "$OUT" | sed 's/^/  /'
VERDICT=$(python3 -c "import json;print(json.load(open('report.json'))['verdict'])")
echo
echo -e "${GREEN}All 11 steps completed.${NC} Report verdict: ${VERDICT}"
if $OFFLINE; then
    echo -e "${YEL}This was an --offline run on synthetic noise: the plumbing is"
    echo -e "verified, the physics is not.${NC}"
fi
$KEEP || { echo "(pass --keep to retain $OUT)"; }
