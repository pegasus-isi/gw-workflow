#!/usr/bin/env python3
"""Shared signal-processing core for the gw-workflow.

Implemented from the equations and stated parameters in references/, per
CLAUDE.md §1.3: every non-obvious constant cites its paper and section.

Sources (references/):
  Usman 2016, CQG 33, 215004        — the PyCBC search pipeline
  Abbott 2016, PRL 116, 061102      — GW150914
  Abbott 2019, PRX 9, 031040        — GWTC-1
  Abbott 2021, SoftwareX 13, 100658 — GWOSC open data

This module is a *support file*, not a transformation: it is registered in the
Replica Catalog and staged into each job's working directory, so wrappers load
it with os.getcwd() on sys.path rather than a __file__-relative import
(PEGASUS.md, "Finding Support Files").
"""

from __future__ import annotations

import logging
import numpy as np
from scipy import signal

logger = logging.getLogger(__name__)

# --- Physical constants (SI) -------------------------------------------------
G = 6.67430e-11
C = 299792458.0
MSUN = 1.988409870698051e30
MSUN_SEC = G * MSUN / C**3          # geometrized solar mass, seconds

# --- Pipeline constants, each traceable to a paper --------------------------
PSD_BLOCK_S = 16.0      # Usman §3.3: inverse PSD is truncated to 16 s, so 16 s
                        # estimation blocks follow naturally.
INV_PSD_TRUNC_S = 16.0  # Usman §3.2: bounds the filter impulse response.
GATE_TAPER_ALPHA = 0.5  # Usman §3.2: Tukey window shape parameter.
GATE_ZERO_S = 1.0       # Usman §3.2: ~1 s zeroed around a short transient.
CLUSTER_WINDOW_S = 1.0  # Usman §3.3: 1 s clustering windows.
COINC_WINDOW_MS = 15.0  # Usman §3.5: 10 ms light travel + 5 ms timing.
LIGHT_TRAVEL_MS = 10.0  # Usman §3.5: max H1-L1 inter-site travel time.


# ===========================================================================
# Frequency-domain helpers
# ===========================================================================
def rfft_norm(x: np.ndarray, dt: float) -> np.ndarray:
    """Discrete approximation to the continuous Fourier transform."""
    return np.fft.rfft(x) * dt


def freqs(n: int, dt: float) -> np.ndarray:
    return np.fft.rfftfreq(n, d=dt)


def band_mask(f: np.ndarray, f_low: float, f_high: float) -> np.ndarray:
    return (f >= f_low) & (f <= f_high)


# ===========================================================================
# High-pass  (SoftwareX §2.2)
# ===========================================================================
def highpass(strain, fs, f_cut, order=4):
    """Zero-phase high-pass below the analysis band.

    SoftwareX §2.2: GWOSC strain is high-passed at 8 Hz "to avoid a number of
    digital signal processing problems related to spectral dynamic range and
    floating point precision limitation that may occur downstream when searching
    in the data". Re-applying it here means the pipeline enforces that condition
    instead of assuming it.

    It matters more than it looks. A large sub-band component leaks across the
    whole spectrum in an unwindowed segment FFT, so the data's in-band power can
    sit orders of magnitude above the windowed Welch PSD the matched filter
    divides by -- and rho scales as 1/sqrt(PSD), so the SNR explodes. A synthetic
    segment without this high-pass produced network SNR ~5.7e6 in testing.
    """
    sos = signal.butter(order, f_cut, btype="highpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, strain)


# ===========================================================================
# PSD estimation and inverse-spectrum truncation  (Usman §3.2-3.3)
# ===========================================================================
def estimate_psd(strain, fs, block_s=PSD_BLOCK_S):
    """Median-averaged Welch PSD over `block_s` blocks, 50% overlap.

    Usman §3.3 decouples PSD estimation from the filter block so many more
    averages are available (127 blocks from a 2048 s stretch vs ~15 overlapping
    blocks in the older scheme), cutting the estimator variance.

    The median average (rather than mean) is the standard robust choice — it
    keeps a loud transient inside one block from inflating the whole estimate.
    scipy applies the median-to-mean bias correction internally.
    """
    nperseg = int(block_s * fs)
    if len(strain) < nperseg:
        raise ValueError(
            f"need >= {block_s}s ({nperseg} samples) to estimate a PSD, got {len(strain)}"
        )
    f, pxx = signal.welch(
        strain, fs=fs, window="hann", nperseg=nperseg, noverlap=nperseg // 2,
        average="median", detrend=False,
    )
    n_blocks = max(1, (len(strain) - nperseg) // (nperseg // 2) + 1)
    return f, pxx, n_blocks


def interpolate_psd(f_psd, psd, f_target):
    """Interpolate a PSD onto the filter's frequency resolution (Usman §3.3)."""
    out = np.interp(f_target, f_psd, psd, left=psd[0], right=psd[-1])
    # Guard the DC/Nyquist ends so 1/S stays finite.
    floor = np.max(psd) * 1e12
    out = np.where((out <= 0) | ~np.isfinite(out), floor, out)
    return out


def truncate_inverse_psd(psd, fs, n_time, trunc_s=INV_PSD_TRUNC_S):
    """Truncate 1/S(f) to `trunc_s` in the time domain.

    Usman §3.2: the matched filter's impulse response is the convolution of the
    template with the inverse FT of 1/S(f). Truncating the inverse PSD to 16 s
    bounds that response — which is what makes the corrupted-edge length
    predictable and gating effective. Without it, a glitch's influence smears
    over a much longer stretch.
    """
    half = int(trunc_s * fs / 2)
    if half < 4 or 2 * half >= n_time:
        # Segment too short to truncate meaningfully; leave the PSD alone rather
        # than silently distorting it.
        return psd.copy()

    # Truncate the AMPLITUDE 1/sqrt(S), not 1/S. This is what LAL's
    # SpectrumInvertTruncate does, and it is not cosmetic: the inverse transform
    # of 1/sqrt(S) is the whitening kernel and is sharply concentrated at t = 0,
    # whereas the transform of 1/S is not. Truncating 1/S directly left the
    # in-band PSD ~700x too small in testing, and because rho scales as
    # 1/sqrt(S), every SNR came out wildly inflated. Measured in-band ratio after
    # truncation: 0.998 on the amplitude vs 0.0015 on 1/S.
    inv_amp = 1.0 / np.sqrt(psd)
    kernel = np.fft.irfft(inv_amp, n=n_time)

    # Flat-top window: unity across the kept region, tapering only at the
    # |t| = trunc_s/2 boundary. The kernel is centred on t = 0, which in a
    # wrapped array means the first and last samples. A full Hann roll-off would
    # attenuate the kernel everywhere instead of only at the cut.
    tuk = signal.windows.tukey(2 * half, alpha=0.25)
    win = np.zeros(n_time)
    win[:half] = tuk[half:]      # t >= 0: flat, then taper out to +trunc_s/2
    win[-half:] = tuk[:half]     # t <  0: taper in from -trunc_s/2, then flat

    amp_trunc = np.abs(np.fft.rfft(kernel * win, n=n_time))
    positive = amp_trunc[amp_trunc > 0]
    if not positive.size:
        return psd.copy()
    amp_trunc[amp_trunc <= 0] = positive.min()
    return 1.0 / amp_trunc ** 2


def check_truncation(psd, psd_trunc, f, f_low, f_high, tol=0.2):
    """Invariant: truncation must preserve the in-band PSD level.

    Returns (ok, ratio). Getting this wrong does not raise an exception anywhere
    -- it silently rescales every SNR -- so the conditioning step asserts it
    explicitly rather than trusting the implementation.
    """
    m = band_mask(f, f_low, f_high)
    if not np.any(m):
        return True, 1.0
    ratio = float(np.median(psd_trunc[m] / psd[m]))
    return abs(ratio - 1.0) <= tol, ratio


def whiten(strain, psd_interp, fs):
    """Whiten in the frequency domain. Used for gating and the spectrogram."""
    n = len(strain)
    dt = 1.0 / fs
    sf = rfft_norm(strain, dt)
    wf = sf / np.sqrt(psd_interp / (2 * dt))
    return np.fft.irfft(wf, n=n) / dt


# ===========================================================================
# Gating  (Usman §3.2)
# ===========================================================================
def find_gate_times(white, fs, threshold, cluster_s=CLUSTER_WINDOW_S):
    """Locate short, loud transients in whitened strain.

    Usman §3.2: whiten, take |whitened strain|, threshold, and cluster the
    peaks. The threshold is a tuning choice the paper does not fix; SPEC.md §9.3
    records it as an open question to be set against the off-source control.
    """
    mag = np.abs(white)
    over = mag > threshold
    if not np.any(over):
        return []
    idx = np.flatnonzero(over)
    gap = int(cluster_s * fs)
    groups, start, prev = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - prev > gap:
            groups.append((start, prev))
            start = i
        prev = i
    groups.append((start, prev))
    # Report the peak sample of each group.
    return [int(a + np.argmax(mag[a : b + 1])) for a, b in groups]


def apply_gates(strain, fs, gate_idx, zero_s=GATE_ZERO_S, alpha=GATE_TAPER_ALPHA):
    """Zero `zero_s` around each gate time with a Tukey window.

    Usman §3.2: a Tukey window (alpha = 0.5) rolls the data smoothly to zero so
    the gate itself does not inject a discontinuity. One glitch left ungated
    produces SNR ~1000 triggers at its own time plus an SNR ~10 skirt spread
    over +/-8 s from the template and inverse-PSD impulse responses.
    """
    out = strain.copy()
    half = int(zero_s * fs / 2)
    for c in gate_idx:
        a, b = max(0, c - half), min(len(out), c + half)
        seg = b - a
        if seg <= 2:
            continue
        # Tukey with alpha: taper region rolls to zero, centre is fully zeroed.
        w = 1.0 - signal.windows.tukey(seg, alpha=alpha)
        out[a:b] *= w
    return out


# ===========================================================================
# Waveform generation  (SPEC M2)
# ===========================================================================
_PYCBC_OK = None


def _pycbc_available():
    global _PYCBC_OK
    if _PYCBC_OK is None:
        try:
            import pycbc.waveform  # noqa: F401
            _PYCBC_OK = True
        except Exception:
            _PYCBC_OK = False
    return _PYCBC_OK


def fd_waveform(m1, m2, f_low, df, f_high, approximant="IMRPhenomD"):
    """Frequency-domain template, normalised arbitrarily (SNR divides it out).

    SPEC M2: waveform generation uses an established library approximant.
    Re-deriving a waveform model is out of scope.

    Default is IMRPhenomD, not SEOBNRv4_ROM. SEOBNRv4_ROM is what the GWTC-1
    PyCBC bank used above 4 Msun and remains selectable, but LALSuite >= 7.25
    resolves its interpolation grids from $LAL_DATA_PATH
    (SEOBNRv4ROM_v3.0.hdf5), which ships in a separate multi-GB data release
    rather than the pip wheel:
        https://git.ligo.org/waveforms/software/lalsuite-waveform-data
        https://zenodo.org/records/14999310
    Without it the call raises and this function would silently degrade to the
    PN toy below. IMRPhenomD is fully analytic, needs no data files, and is a
    complete inspiral-merger-ringdown aligned-spin model -- appropriate for
    these masses and permitted by SPEC.md M2.

    Returns (h_of_f, used_fallback). When PyCBC is unavailable the restricted
    post-Newtonian inspiral below is substituted so run_manual.sh can smoke-test
    on a bare host -- it is NOT adequate for these masses (the merger and
    ringdown carry much of the in-band SNR, which is exactly why GWTC-1 used
    SEOBNRv4), so every consumer must propagate used_fallback and refuse to
    report a validated result from it.
    """
    n_f = int(f_high / df) + 1
    if _pycbc_available():
        from pycbc.waveform import get_fd_waveform
        try:
            hp, _ = get_fd_waveform(
                approximant=approximant, mass1=m1, mass2=m2,
                spin1z=0.0, spin2z=0.0,       # SPEC M3: non-spinning bank
                delta_f=df, f_lower=f_low, f_final=f_high,
            )
            h = np.asarray(hp.data)[:n_f] if len(hp) >= n_f else np.pad(
                np.asarray(hp.data), (0, n_f - len(hp))
            )
            return h.astype(np.complex128), False
        except Exception as exc:      # pragma: no cover - library-dependent
            logger.warning("pycbc waveform failed (%s); using PN fallback", exc)
    return _pn_inspiral(m1, m2, f_low, df, f_high, n_f), True


def _pn_inspiral(m1, m2, f_low, df, f_high, n_f):
    """Restricted post-Newtonian inspiral -- SMOKE TEST ONLY.

    Leading-order amplitude f^{-7/6} with 0PN phase, truncated at the ISCO
    frequency. Usman §2 notes the inspiral of a low-mass binary is well modelled
    by post-Newtonian approximations; for the BBH masses in SPEC.md it is not,
    which is why anything built on this path is marked used_fallback=True.
    """
    f = np.arange(n_f) * df
    mc = chirp_mass(m1, m2) * MSUN_SEC
    mt = (m1 + m2) * MSUN_SEC
    f_isco_hz = 1.0 / (6 ** 1.5 * np.pi * mt)
    ok = (f >= f_low) & (f <= min(f_high, f_isco_hz)) & (f > 0)
    h = np.zeros(n_f, dtype=np.complex128)
    x = np.pi * mc * f[ok]
    phase = -0.75 * (8 * x) ** (-5.0 / 3.0)
    h[ok] = (mc ** (5.0 / 6.0)) * f[ok] ** (-7.0 / 6.0) * np.exp(1j * phase)
    return h


# ===========================================================================
# Matched filter  (Usman §2, eq. 1-2)
# ===========================================================================
def sigma_sq(h, psd, df, mask):
    """<h|h> -- the template normalisation of Usman eq. (2)."""
    return float(4.0 * df * np.sum((np.abs(h[mask]) ** 2) / psd[mask]))


def matched_filter_snr(data_fd, h, psd, df, n_time, mask):
    """Complex matched-filter SNR time series.

    Usman eq. (1) maximises over the template's two orthogonal phases; taking
    the modulus of the complex correlation below is the equivalent closed form,
    so amplitude and phase (hence sky position and orientation) are maximised
    analytically rather than searched.

    We want the one-sided complex sum

        z(t_n) = 4 df  sum_{k=0}^{N/2}  q_k exp(2 pi i k n / N)

    which is a full-length complex inverse transform of the one-sided product
    zero-padded to N -- NOT a pair of real irfft calls. irfft interprets its
    input as a Hermitian spectrum and so returns
    q_0 + 2*sum_{k>=1} Re(q_k ...), i.e. exactly TWICE the wanted sum for real
    q. A zero-noise injection test (h filtered against itself, where the peak
    must equal sqrt(<h|h>)) recovers a ratio of 1.000 with this construction and
    2.000 with the irfft pair; run_manual.sh re-runs that check.
    """
    qh = np.zeros(len(data_fd), dtype=np.complex128)
    qh[mask] = data_fd[mask] * np.conj(h[mask]) / psd[mask]
    q_full = np.zeros(n_time, dtype=np.complex128)
    q_full[: len(qh)] = qh
    # numpy's ifft carries a 1/N, so multiply it back out.
    z = 4.0 * df * n_time * np.fft.ifft(q_full)
    s = np.sqrt(sigma_sq(h, psd, df, mask))
    if s <= 0:
        return np.zeros(n_time)
    return np.abs(z) / s


def overlap(h1, h2, psd, df, mask):
    """Normalised overlap in [0, 1] -- the bank's minimal-match measure."""
    n1, n2 = sigma_sq(h1, psd, df, mask), sigma_sq(h2, psd, df, mask)
    if n1 <= 0 or n2 <= 0:
        return 0.0
    x = 4.0 * df * np.sum(h1[mask] * np.conj(h2[mask]) / psd[mask])
    return float(np.abs(x) / np.sqrt(n1 * n2))


# ===========================================================================
# Trigger clustering  (Usman §3.3)
# ===========================================================================
def cluster_triggers(snr, fs, threshold, window_s=CLUSTER_WINDOW_S):
    """Threshold then keep local maxima of the SNR time series.

    Usman §3.3: divide into `window_s` windows, take each window's maximum, then
    keep it only if it exceeds both neighbouring windows -- so a real signal's
    single narrow peak survives and a broad glitch does not produce a comb.
    """
    w = max(1, int(window_s * fs))
    n_win = len(snr) // w
    if n_win == 0:
        return []
    trimmed = snr[: n_win * w].reshape(n_win, w)
    loc_idx = np.argmax(trimmed, axis=1)
    loc_max = trimmed[np.arange(n_win), loc_idx]
    out = []
    for i in range(n_win):
        if loc_max[i] < threshold:
            continue
        prev_v = loc_max[i - 1] if i > 0 else -np.inf
        next_v = loc_max[i + 1] if i < n_win - 1 else -np.inf
        if loc_max[i] >= prev_v and loc_max[i] >= next_v:
            out.append((int(i * w + loc_idx[i]), float(loc_max[i])))
    return out


# ===========================================================================
# Template bank placement  (Usman §3.1)
# ===========================================================================
def harmonic_mean_psd(psds):
    """Usman eq. (7): per-frequency harmonic mean of N power spectra."""
    stack = np.vstack(psds)
    return stack.shape[0] / np.sum(1.0 / stack, axis=0)


def place_bank_row(q, mc_min, mc_max, minimal_match, psd, freqs, df, mask,
                   f_low, f_high, approximant, max_templates=20000,
                   progress=None):
    """Place one constant-q row of the bank by walking chirp mass upward.

    Each next template is chosen by bisection so the MEASURED overlap with its
    neighbour equals `minimal_match`, so the bank reports a measured minimal
    match rather than asserting one (SPEC.md B4, V2d).

    Rows at fixed q are independent of one another, which is what lets the bank
    build fan out across q. Shared here so the serial and parallel paths cannot
    drift apart.

    Returns (templates, adjacent_overlaps, used_fallback, truncated).
    """
    used_fallback = [False]

    def wf(mc):
        m1, m2 = mass1_mass2(mc, q)
        h, fb = fd_waveform(m1, m2, f_low, df, f_high, approximant)
        if fb:
            used_fallback[0] = True
        if len(h) < len(freqs):
            h = np.pad(h, (0, len(freqs) - len(h)))
        return h[: len(freqs)]

    def ov(a, b):
        return overlap(wf(a), wf(b), psd, df, mask)

    templates, adjacent, mc = [], [], float(mc_min)
    while mc <= mc_max and len(templates) < max_templates:
        m1, m2 = mass1_mass2(mc, q)
        templates.append((float(mc), float(q), float(m1), float(m2)))
        if progress and len(templates) % progress == 0:
            logger.info("  q=%.3f: %d templates (Mc=%.3f)", q, len(templates), mc)

        lo, hi, guard = mc, mc * 1.5, 0
        while ov(mc, hi) > minimal_match and hi < mc_max * 2 and guard < 30:
            lo, hi = hi, hi * 1.5
            guard += 1
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            if ov(mc, mid) > minimal_match:
                lo = mid
            else:
                hi = mid
            if (hi - lo) / mc < 1e-4:
                break
        nxt = 0.5 * (lo + hi)
        if nxt <= mc * (1 + 1e-6):
            logger.warning("spacing collapsed at Mc=%.3f q=%.3f; ending row", mc, q)
            break
        adjacent.append(ov(mc, nxt))
        mc = nxt

    return (templates, adjacent, bool(used_fallback[0]),
            len(templates) >= max_templates)


def bank_adequacy(mcs, ref_events, tier_a, tol):
    """SPEC.md B6: a template within `tol` of every Tier A event's expected
    detector-frame chirp mass, plus the local spacing V2d needs."""
    by_name = {e["name"]: e for e in ref_events}
    mcs = np.asarray(mcs, dtype=float)
    srt = np.sort(mcs)
    records, failures = [], []
    for name in tier_a:
        ev = by_name.get(name)
        if ev is None:
            failures.append(f"{name} absent from the reference file")
            continue
        sp = ev["source_parameters"]
        mc_det = (sp["chirp_mass_source_msun"]["median"]
                  * (1.0 + sp["redshift"]["median"]))
        i = int(np.argmin(np.abs(mcs - mc_det)))
        frac = abs(mcs[i] - mc_det) / mc_det
        j = int(np.argmin(np.abs(srt - mc_det)))
        nbr = [srt[k] for k in (j - 1, j + 1) if 0 <= k < len(srt)]
        rec = {
            "event": name,
            "mc_detector_expected": round(mc_det, 4),
            "nearest_template_mc": round(float(mcs[i]), 4),
            "fractional_offset": round(float(frac), 5),
            "within_tolerance": bool(frac <= tol),
            "local_spacing_fractional": (
                round(float(min(abs(n - srt[j]) for n in nbr) / srt[j]), 5)
                if nbr else None),
        }
        records.append(rec)
        if not rec["within_tolerance"]:
            failures.append(
                f"{name}: nearest template Mc={mcs[i]:.3f} is {frac*100:.2f}% "
                f"from the expected {mc_det:.3f} (tolerance {tol*100:.0f}%)")
    return records, failures


# ===========================================================================
# Mass / frequency relations
# ===========================================================================
def chirp_mass(m1, m2):
    """Usman §2 / PRL 116 061102 eq. for the chirp mass."""
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def mass1_mass2(mc, q):
    """Invert (chirp mass, mass ratio q = m1/m2 >= 1) to component masses."""
    eta = q / (1.0 + q) ** 2
    mt = mc / eta ** 0.6
    m1 = mt * q / (1.0 + q)
    return m1, mt - m1


def f_isco(m_total_msun):
    """ISCO frequency in Hz -- the ~4400/M scale used throughout SPEC.md."""
    return 1.0 / (6 ** 1.5 * np.pi * m_total_msun * MSUN_SEC)


def detector_frame(mc_source, z):
    """Detector-frame chirp mass. SPEC.md §6.2 -- the conversion that must not
    be skipped: (1+z) runs 1.07-1.49 across this event set."""
    return mc_source * (1.0 + z)
