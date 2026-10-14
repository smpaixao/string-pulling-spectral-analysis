"""
=============================================================================
STRING PULLING KINEMATIC ANALYSIS — AUGMENTED VERSION
=============================================================================
Script  : string_pulling_AUGMENTED_v2.3_2026-09-17.py
Version : 2.3   (created 2026-09-17)
Based on: string_pulling_FINAL.py (Sonia's original, unchanged) and FINAL_methods.txt

Corrected and extended version of string_pulling_FINAL.py.
Input format is unchanged (raw DeepLabCut Excel files, one per animal, trials
separated by blank rows). Frame rates are read from data/frame_rates.xlsx
(and optionally data/frame_rates_per_video.xlsx).

VERSION HISTORY
  v1.0  2026-09-15  First augmented version (fixes 1-6 below, cycle measures,
                    sample entropy, 12-measure statistics; stops >= 2 s).
  v1.1  2026-09-15  Optional per-video frame rates (frame_rates_per_video.xlsx).
  v1.2  2026-09-15  Aligned with FINAL_methods.txt: excluded mice 22367, 26806,
                    26829; pauses = 10th-percentile threshold, >= 150 ms.
  v2.0  2026-09-16  Statistics reduced to the dominant (peak) frequency of both
                    hands on Y (vertical) and X (lateral); Shapiro-Wilk rule
                    chooses Welch t-test / F-test or Mann-Whitney U /
                    Brown-Forsythe; effect sizes with bootstrap 95% CI; Holm
                    across the four tests; X-axis spectra added; group colours
                    Control gray / Zic2 teal; outputs go to a versioned folder
                    with run_info.txt. Other measures are still saved per mouse
                    but not tested.
  v2.1  2026-09-17  Figures only: all four panels of group_spectra.png and of
                    group_metrics.png share one y-axis range, so panels can be
                    compared directly. No change to any analysis or statistic;
                    numerical outputs are identical to v2.0.
  v2.2  2026-09-17  Figures only: each figure is also written as PDF and SVG for
                    editing in Illustrator, with text kept as text rather than
                    outlines (pdf.fonttype 42, svg.fonttype 'none'). Formats are
                    set in FIG_FORMATS. No change to any analysis or statistic;
                    numerical outputs are identical to v2.0.
  v2.3  2026-09-17  Rhythmicity redefined and tested. It was the power in the
                    single peak bin over the band total, which on the 0.05 Hz
                    interpolated grid is set by the grid spacing rather than by
                    the movement (the true resolution from a 1.5 s window is
                    about 0.67 Hz, so neighbouring bins are not independent).
                    It is now the power within RHYTHM_HALFWIDTH_HZ of the peak
                    over the band total, which is resolution-independent and
                    reads as the fraction of movement power concentrated near
                    the dominant frequency. The four rhythmicity measures now
                    appear in group_statistics.csv and group_metrics.png
                    alongside the four dominant frequencies. Dominant-frequency
                    values and their statistics are unchanged from v2.0.

WHAT CHANGED COMPARED WITH string_pulling_FINAL.py
  [FIX 1] Spectra from different trials / frame rates are interpolated onto one
          common frequency grid before averaging (the original averaged arrays
          position by position although their frequency spacing differed).
  [FIX 2] Same Welch window length (PSD_WIN_S) for every trial and frame rate,
          so frequency resolution and rhythmicity do not depend on trial length.
  [FIX 3] Slow drift removed before spectral analysis (linear detrend +
          high-pass at HIGHPASS_HZ); peak search starts at FREQ_MIN = 1.5 Hz.
  [FIX 4] Pause detection follows FINAL_methods.txt: 10th-percentile adaptive
          velocity threshold, pauses >= 150 ms removed (the script used the
          25th percentile). STOP_PERCENTILE / STOP_MIN_S can be changed.
  [FIX 5] Data are never stitched across a stop; each continuous pulling
          segment is analysed separately.
  [FIX 6] Frame rates read from files, not typed into the script: one per video
          from frame_rates_per_video.xlsx, otherwise one per mouse from
          frame_rates.xlsx.
  [STATS] Group statistics for the dominant (peak) frequency of the preferred
          and non-preferred hand, on the Y (vertical) and X (lateral) axes.
          Rule used throughout the paper: Shapiro-Wilk in each group; if both
          groups are normal -> Welch t-test (typical value) and F-test (spread),
          otherwise -> Mann-Whitney U (typical value) and Brown-Forsythe (spread).
          Effect size with a bootstrap 95% CI; Holm-adjusted p across the four tests.
  [EXTRA] Cycle measures, sample entropy, inter-hand correlation and lateral
          variability are still calculated per mouse and saved in
          summary_per_animal.csv, but they are not tested or plotted.

All parameters are in the SETTINGS block. They are defensible defaults, not
validated values: choose them on biological grounds before comparing groups.

Run from the scripts/ folder:  python string_pulling_AUGMENTED_v2.0_2026-09-16.py
Requires: pandas numpy scipy matplotlib openpyxl
=============================================================================
"""

import os, sys, warnings, datetime, platform, scipy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal, stats

warnings.filterwarnings('ignore', category=RuntimeWarning)   # empty-slice means only

# =============================================================================
# SETTINGS
# =============================================================================
HERE        = os.path.dirname(os.path.abspath(__file__))
PROJECT     = os.path.dirname(HERE)
DATA_FOLDER = os.path.join(PROJECT, 'data')
SCRIPT_NAME  = 'string_pulling_AUGMENTED'
VERSION      = '2.3'
VERSION_DATE = '2026-09-17'
OUT_FOLDER  = os.path.join(PROJECT, 'results', f'augmented_v{VERSION}_{VERSION_DATE}')
FRAME_RATE_FILE = os.path.join(DATA_FOLDER, 'frame_rates.xlsx')
# Optional: one frame rate per video. Columns: animal | trial | fps
# ("trial" = position of the block in the mouse's Excel file, 1 = first block).
# Any video listed here overrides the per-mouse rate in frame_rates.xlsx.
VIDEO_RATE_FILE = os.path.join(DATA_FOLDER, 'frame_rates_per_video.xlsx')

GROUPS = [("Control", "#808080"), ("Zic2", "#1A9EA6")]   # Control gray, Zic2 teal

# Figure formats. 'png' to look at, 'pdf' / 'svg' to edit in Illustrator.   [v2.2]
FIG_FORMATS = ['png', 'pdf', 'svg']
matplotlib.rcParams['pdf.fonttype'] = 42      # keep text editable, not outlines
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# Animals excluded in FINAL_methods.txt (fewer than 40 pulling bouts)
EXCLUDED_ANIMALS = ['22367', '26806', '26829']

# Pauses, as described in FINAL_methods.txt                            [FIX 4]
STOP_MIN_S         = 0.15    # pauses at least this long are removed
STOP_PERCENTILE    = 10      # adaptive velocity threshold (percentile)
STOP_ABS_MIN_REF   = 10.0    # px/frame at 146 fps, as original
STOP_REF_FPS       = 146.0
STOP_SMOOTH_MS     = 150
MIN_SEGMENT_S      = 1.5     # shortest continuous pulling segment analysed

# Spectral analysis
PSD_WIN_S   = 1.5            # Welch window length (s)                  [FIX 2]
HIGHPASS_HZ = 1.0            # drift removal                           [FIX 3]
FREQ_MIN    = 1.5            # dominant-frequency search band (Hz)     [FIX 3]
FREQ_MAX    = 8.0
RHYTHM_HALFWIDTH_HZ = 0.5    # rhythmicity = power within +/- this of the peak [v2.3]
FREQ_GRID   = np.round(np.arange(0, 12.0001, 0.05), 4)                # [FIX 1]

# Cycle and entropy metrics
LOWPASS_HZ     = 15.0        # smoothing before cycle detection
CYCLE_PROM_SD  = 0.5         # minimum trough prominence (in SDs)
CYCLE_POINTS   = 50          # each cycle resampled to this many points
MIN_CYCLES     = 3
ENTROPY_FS     = 100.0       # common time base for sample entropy
ENTROPY_M      = 2
ENTROPY_R      = 0.2         # tolerance, x SD

# Statistics
NORMALITY_ALPHA  = 0.05      # Shapiro-Wilk p above this in BOTH groups -> normal
T_TEST_EQUAL_VAR = False     # False = Welch t-test, True = Student t-test
N_BOOT           = 10000     # bootstrap resamples for the 95% CI of the group difference

REQUIRED = ['x_preferred', 'y_preferred', 'x_non_preferred', 'y_non_preferred']

# =============================================================================
# LOADING
# =============================================================================
def load_frame_rates(path):
    fr = pd.read_excel(path)
    fr.columns = ['animal', 'group', 'fps']
    return {str(a): float(f) for a, f in zip(fr.animal, fr.fps)}

def load_video_rates(path):
    """{(animal, trial): fps} from the optional per-video file; {} if absent."""
    if not os.path.exists(path):
        return {}
    vr = pd.read_excel(path)
    vr.columns = [c.strip().lower() for c in vr.columns]
    missing = {'animal', 'trial', 'fps'} - set(vr.columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}')
    vr = vr.dropna(subset=['animal', 'trial', 'fps'])
    return {(str(int(a)), int(t)): float(f) for a, t, f in zip(vr.animal, vr.trial, vr.fps)}

def load_trials(filepath):
    """Split sheet into trials at blank rows (as original)."""
    df = pd.read_excel(filepath, sheet_name=0, header=0)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f'{filepath}: missing columns {missing}')
    blank = df[REQUIRED].isnull().all(axis=1).values
    trials, start = [], None
    for i, b in enumerate(np.append(blank, True)):
        if not b and start is None:
            start = i
        elif b and start is not None:
            t = df.iloc[start:i][REQUIRED].apply(pd.to_numeric, errors='coerce').dropna()
            if len(t): trials.append(t.reset_index(drop=True))
            start = None
    return trials

# =============================================================================
# STOPS AND SEGMENTS
# =============================================================================
def detect_pulling(yp, ynp, fs):
    """True = pulling. Low-velocity runs shorter than STOP_MIN_S count as pulling."""
    win = max(3, int(STOP_SMOOTH_MS / 1000 * fs)) | 1
    v = (np.convolve(np.abs(np.gradient(yp)),  np.ones(win) / win, mode='same') +
         np.convolve(np.abs(np.gradient(ynp)), np.ones(win) / win, mode='same'))
    thr = max(np.percentile(v, STOP_PERCENTILE), STOP_ABS_MIN_REF * STOP_REF_FPS / fs)
    pulling = v > thr
    for s, e in runs(~pulling):
        if (e - s) < STOP_MIN_S * fs:
            pulling[s:e] = True
    return pulling

def runs(mask):
    edges = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    return list(zip(np.where(edges == 1)[0], np.where(edges == -1)[0]))

# =============================================================================
# SIGNAL HELPERS
# =============================================================================
def zscore(a):
    s = np.std(a)
    return (a - np.mean(a)) / s if s > 0 else a - np.mean(a)

def bandpass(x, fs, lo=HIGHPASS_HZ, hi=LOWPASS_HZ):
    b, a = signal.butter(2, [lo / (fs / 2), min(hi / (fs / 2), 0.99)], btype='band')
    return signal.filtfilt(b, a, signal.detrend(x))

def welch_on_grid(x, fs):
    nper = int(round(PSD_WIN_S * fs))
    f, p = signal.welch(x, fs=fs, nperseg=nper, noverlap=nper // 2,
                        nfft=int(round(fs / 0.05)), window='hann', detrend='linear')
    return np.interp(FREQ_GRID, f, p)

def peak_metrics(psd):
    """Dominant frequency, and the share of band power within
    RHYTHM_HALFWIDTH_HZ of it (resolution-independent).            [v2.3]"""
    band = (FREQ_GRID >= FREQ_MIN) & (FREQ_GRID <= FREQ_MAX)
    f_band, p_band = FREQ_GRID[band], psd[band]
    i = int(np.argmax(p_band))
    f_peak = float(f_band[i])
    near = np.abs(f_band - f_peak) <= RHYTHM_HALFWIDTH_HZ
    return f_peak, float(p_band[near].sum() / p_band.sum())

def sample_entropy(x, m=ENTROPY_M, r=ENTROPY_R):
    x = zscore(np.asarray(x, float)); tol = r * np.std(x); n = len(x)
    if n < 10 * (m + 1): return np.nan
    def count(mm):
        T = np.lib.stride_tricks.sliding_window_view(x, mm)[:n - m]
        c = 0
        for i in range(len(T) - 1):
            c += np.sum(np.max(np.abs(T[i + 1:] - T[i]), axis=1) <= tol)
        return c
    B, A = count(m), count(m + 1)
    return float(-np.log(A / B)) if A > 0 and B > 0 else np.nan

def cycles(x_filt, fs):
    """Troughs of the filtered signal -> list of (start, end) cycle indices."""
    tr, _ = signal.find_peaks(-x_filt, distance=int(fs / FREQ_MAX),
                              prominence=CYCLE_PROM_SD * np.std(x_filt))
    out = []
    for s, e in zip(tr[:-1], tr[1:]):
        if 1 / FREQ_MAX <= (e - s) / fs <= 1 / FREQ_MIN:
            out.append((s, e))
    return out

def interhand_corr(yp, ynp, fs):
    """As original: correlation at lag 0 and most negative correlation within ±500 ms."""
    if len(yp) < 10: return np.nan, np.nan
    y1, y2 = zscore(signal.detrend(yp)), zscore(signal.detrend(ynp))
    corr = np.correlate(y1, y2, mode='full') / len(y1)
    lags = np.arange(-(len(y1) - 1), len(y1)); cen = len(y1) - 1; ml = int(0.5 * fs)
    cw, lw = corr[max(0, cen - ml):cen + ml + 1], lags[max(0, cen - ml):cen + ml + 1]
    return float(cw[np.argmin(np.abs(lw))]), float(cw.min())

# =============================================================================
# PER TRIAL / PER ANIMAL
# =============================================================================
def analyse_segment(seg, fs):
    """seg: dict of raw arrays for one continuous pulling segment."""
    r = {}
    yp, ynp = seg['yp'], seg['ynp']
    for key, sig in [('Y_nonpref', ynp), ('Y_pref', yp),
                     ('X_nonpref', seg['xnp']), ('X_pref', seg['xp'])]:
        r['psd_' + key] = welch_on_grid(zscore(bandpass(sig, fs, hi=fs / 2.2)), fs)
    r['corr_lag0'], r['corr_anti'] = interhand_corr(yp, ynp, fs)
    r['x_sd_nonpref'] = float(np.std(seg['xnp']))
    r['x_sd_pref']    = float(np.std(seg['xp']))

    f = bandpass(ynp, fs)
    cyc = cycles(f, fs)
    if len(cyc) >= MIN_CYCLES:
        dur = np.array([(e - s) / fs for s, e in cyc])
        shapes = np.array([zscore(np.interp(np.linspace(s, e, CYCLE_POINTS),
                                            np.arange(len(f)), f)) for s, e in cyc])
        tmpl = shapes.mean(axis=0)
        r['cycle_rate_Hz']   = float(1 / np.median(dur))
        r['cycle_dur_CV']    = float(np.std(dur) / np.mean(dur))
        r['cycle_amp_px']    = float(np.median([f[s:e + 1].max() - f[s:e + 1].min() for s, e in cyc]))
        r['cycle_shape_sim'] = float(np.mean([np.corrcoef(c, tmpl)[0, 1] for c in shapes]))
        r['n_cycles']        = len(cyc)
    t_new = np.arange(0, len(ynp) / fs, 1 / ENTROPY_FS)
    y100 = np.interp(t_new, np.arange(len(ynp)) / fs, bandpass(ynp, fs))
    r['sample_entropy'] = sample_entropy(y100)
    return r

SPECTRA = ['Y_nonpref', 'Y_pref', 'X_nonpref', 'X_pref']   # axis_hand

SCALARS = ['corr_lag0', 'corr_anti', 'x_sd_nonpref', 'x_sd_pref', 'cycle_rate_Hz',
           'cycle_dur_CV', 'cycle_amp_px', 'cycle_shape_sim', 'sample_entropy']

def process_animal(path, name, group, fs_mouse, video_rates=None):
    rows = []
    store = {k: [] for k in SPECTRA}
    n_trials, n_video_rates = 0, 0
    video_rates = video_rates or {}
    for ti, t in enumerate(load_trials(path)):
        n_trials += 1
        fs = video_rates.get((name, ti + 1), fs_mouse)      # per-video rate if given
        n_video_rates += (name, ti + 1) in video_rates
        yp, ynp = -t.y_preferred.values, -t.y_non_preferred.values
        xp, xnp = t.x_preferred.values, t.x_non_preferred.values
        pulling = detect_pulling(yp, ynp, fs)
        for si, (s, e) in enumerate(runs(pulling)):
            dur = (e - s) / fs
            if dur < max(MIN_SEGMENT_S, PSD_WIN_S): continue
            seg = dict(yp=yp[s:e], ynp=ynp[s:e], xp=xp[s:e], xnp=xnp[s:e])
            r = analyse_segment(seg, fs)
            for k in SPECTRA:
                psd = r.pop('psd_' + k); store[k].append((psd, dur))
                r['dom_Hz_' + k], r['rhythm_' + k] = peak_metrics(psd)
            rows.append(dict(animal=name, group=group, fps=fs, trial=ti + 1,
                             segment=si + 1, duration_s=round(dur, 2), **r))
    tr = pd.DataFrame(rows)
    unknown = [k for k in video_rates if k[0] == name and k[1] > n_trials]
    if unknown:
        print(f'WARNING: {name}: per-video rates given for trials {[k[1] for k in unknown]}, '
              f'but the file has only {n_trials} trials')
    m = dict(animal=name, group=group, fps=fs_mouse, n_trials_total=n_trials,
             n_trials_with_video_fps=n_video_rates,
             n_segments_used=len(tr),
             seconds_used=float(tr.duration_s.sum()) if len(tr) else 0.0)
    spectra = {}
    for key in SPECTRA:
        if store[key]:
            avg = np.average([p for p, _ in store[key]], axis=0, weights=[d for _, d in store[key]])
            m[f'dom_Hz_{key}'], m[f'rhythm_{key}'] = peak_metrics(avg)
            spectra[key] = avg
        else:
            m[f'dom_Hz_{key}'] = m[f'rhythm_{key}'] = np.nan; spectra[key] = None
    m['trial_SD_dom_Hz_Y_nonpref'] = float(tr.dom_Hz_Y_nonpref.std()) if len(tr) >= 2 else np.nan
    for c in SCALARS:
        m[c] = float(np.average(tr[c].dropna(), weights=tr.loc[tr[c].notna(), 'duration_s'])) \
               if len(tr) and c in tr and tr[c].notna().any() else np.nan
    return m, tr, spectra

# =============================================================================
# STATISTICS
# =============================================================================
METRICS = [  # column, label, axis unit  (the tested measures)        [v2.3]
    ('dom_Hz_Y_nonpref', 'Dominant frequency, Y (vertical), non-preferred hand (Hz)', 'Hz'),
    ('dom_Hz_Y_pref',    'Dominant frequency, Y (vertical), preferred hand (Hz)',     'Hz'),
    ('dom_Hz_X_nonpref', 'Dominant frequency, X (lateral), non-preferred hand (Hz)',  'Hz'),
    ('dom_Hz_X_pref',    'Dominant frequency, X (lateral), preferred hand (Hz)',      'Hz'),
    ('rhythm_Y_nonpref', 'Rhythmicity, Y (vertical), non-preferred hand',      'fraction'),
    ('rhythm_Y_pref',    'Rhythmicity, Y (vertical), preferred hand',          'fraction'),
    ('rhythm_X_nonpref', 'Rhythmicity, X (lateral), non-preferred hand',       'fraction'),
    ('rhythm_X_pref',    'Rhythmicity, X (lateral), preferred hand',           'fraction'),
]

def hedges_g(a, b):
    """Standardised mean difference (Zic2 - Control), small-sample corrected."""
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    d = (np.mean(b) - np.mean(a)) / sp if sp > 0 else np.nan
    return d * (1 - 3 / (4 * (na + nb) - 9))

def rank_biserial(a, b):
    """Rank-based effect size in [-1, 1]; > 0 means Zic2 values tend to be higher."""
    u_b = stats.mannwhitneyu(b, a, alternative='two-sided').statistic
    return 2 * u_b / (len(a) * len(b)) - 1

def bootstrap_ci(a, b, stat, n_boot=N_BOOT, seed=0):
    rng = np.random.default_rng(seed)
    diffs = [stat(rng.choice(b, len(b))) - stat(rng.choice(a, len(a))) for _ in range(n_boot)]
    return np.percentile(diffs, [2.5, 97.5])

def holm(p):
    p = np.asarray(p, float); order = np.argsort(p); adj = np.empty_like(p); run = 0
    for rank, i in enumerate(order):
        run = max(run, min(1, (len(p) - rank) * p[i])); adj[i] = run
    return adj

def compare_groups(a, b):
    """Shapiro-Wilk rule -> parametric or non-parametric tests for centre and spread."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    sw_a, sw_b = stats.shapiro(a).pvalue, stats.shapiro(b).pvalue
    normal = sw_a > NORMALITY_ALPHA and sw_b > NORMALITY_ALPHA
    out = dict(shapiro_p_ctrl=sw_a, shapiro_p_zic2=sw_b, normal='yes' if normal else 'no')
    if normal:
        out['centre_test'] = 'Welch t-test' if not T_TEST_EQUAL_VAR else 'Student t-test'
        out['p_centre'] = stats.ttest_ind(a, b, equal_var=T_TEST_EQUAL_VAR).pvalue
        out['effect_size'], out['effect_size_type'] = hedges_g(a, b), 'Hedges g'
        out['difference'], out['difference_type'] = np.mean(b) - np.mean(a), 'mean (Zic2 - Control)'
        lo, hi = bootstrap_ci(a, b, np.mean)
        out['spread_test'] = 'F-test'
        F = np.var(b, ddof=1) / np.var(a, ddof=1)
        out['p_spread'] = 2 * min(stats.f.cdf(F, len(b) - 1, len(a) - 1), stats.f.sf(F, len(b) - 1, len(a) - 1))
    else:
        ties = len(np.unique(np.r_[a, b])) < len(a) + len(b)
        method = 'asymptotic' if ties else 'exact'
        out['centre_test'] = f'Mann-Whitney U ({method})'
        out['p_centre'] = stats.mannwhitneyu(a, b, alternative='two-sided', method=method).pvalue
        out['effect_size'], out['effect_size_type'] = rank_biserial(a, b), 'rank-biserial r'
        out['difference'], out['difference_type'] = np.median(b) - np.median(a), 'median (Zic2 - Control)'
        lo, hi = bootstrap_ci(a, b, np.median)
        out['spread_test'] = 'Brown-Forsythe'
        out['p_spread'] = stats.levene(a, b, center='median').pvalue
    out['diff_ci95_low'], out['diff_ci95_high'] = lo, hi
    return out

def group_stats(df):
    g1, g2 = GROUPS[0][0], GROUPS[1][0]; rows = []
    for col, label, _unit in METRICS:
        a = df.loc[df.group == g1, col].dropna().values; b = df.loc[df.group == g2, col].dropna().values
        row = dict(metric=col, label=label,
                   n_ctrl=len(a), ctrl_mean=np.mean(a), ctrl_sd=np.std(a, ddof=1), ctrl_sem=stats.sem(a),
                   ctrl_median=np.median(a), ctrl_q1=np.percentile(a, 25), ctrl_q3=np.percentile(a, 75),
                   n_zic2=len(b), zic2_mean=np.mean(b), zic2_sd=np.std(b, ddof=1), zic2_sem=stats.sem(b),
                   zic2_median=np.median(b), zic2_q1=np.percentile(b, 25), zic2_q3=np.percentile(b, 75))
        row.update(compare_groups(a, b))
        rows.append(row)
    st = pd.DataFrame(rows)
    st['p_centre_holm'] = holm(st.p_centre.values)
    st['p_spread_holm'] = holm(st.p_spread.values)
    return st

# =============================================================================
# FIGURES
# =============================================================================
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'

def style(ax):
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(axis='y', color=GRID, lw=0.8)
    ax.set_axisbelow(True)

def save_figure(path):
    """Write the current figure once per entry in FIG_FORMATS."""      # [v2.2]
    stem = os.path.splitext(path)[0]
    for fmt in FIG_FORMATS:
        plt.savefig(f'{stem}.{fmt}', format=fmt, dpi=200)

def plot_spectra(spectra, path):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True)
    mask = FREQ_GRID <= 10
    panels = [('Y_pref', 'Y (vertical) · preferred hand'), ('Y_nonpref', 'Y (vertical) · non-preferred hand'),
              ('X_pref', 'X (lateral) · preferred hand'), ('X_nonpref', 'X (lateral) · non-preferred hand')]
    for ax, (key, title) in zip(axes.flat, panels):
        for grp, col in GROUPS:
            arr = np.array([s[key] for g, s in spectra if g == grp and s[key] is not None])
            avg, sem = arr.mean(0), arr.std(0, ddof=1) / np.sqrt(len(arr))
            ax.fill_between(FREQ_GRID[mask], (avg - sem)[mask], (avg + sem)[mask], color=col, alpha=.12, lw=0)
            ax.plot(FREQ_GRID[mask], avg[mask], color=col, lw=2, label=f'{grp} (n={len(arr)})')
        ax.axvspan(FREQ_MIN, FREQ_MAX, color=GRID, alpha=.35, lw=0)
        ax.set_title(title, color=INK, fontsize=10, loc='left')
        ax.set_ylabel('Power (z-scored position)', color=INK2, fontsize=9)
        style(ax)
    lo = min(ax.get_ylim()[0] for ax in axes.flat)          # common y range [v2.1]
    hi = max(ax.get_ylim()[1] for ax in axes.flat)
    for ax in axes.flat: ax.set_ylim(lo, hi)
    for ax in axes[1]: ax.set_xlabel('Frequency (Hz)', color=INK2, fontsize=9)
    axes[0, 1].legend(frameon=False, fontsize=8, labelcolor=INK2)
    plt.tight_layout(); save_figure(path); plt.close()

def plot_metrics(df, st, path):
    ncol = 4                                                          # [v2.3]
    nrow = int(np.ceil(len(METRICS) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.5 * nrow), squeeze=False)
    rng = np.random.default_rng(1)
    for ax, (col, label, unit), (_, s) in zip(axes.flat, METRICS, st.iterrows()):
        for i, (grp, c) in enumerate(GROUPS):
            v = df.loc[df.group == grp, col].dropna().values
            ax.scatter(i + rng.uniform(-.12, .12, len(v)), v, s=30, color=c,
                       edgecolor='white', linewidth=1, zorder=3)
            ax.hlines(np.median(v) if s.normal == 'no' else np.mean(v), i - .25, i + .25, color=INK, lw=2, zorder=4)
        ax.set_xticks([0, 1]); ax.set_xticklabels([g for g, _ in GROUPS], color=INK2, fontsize=9)
        ax.set_xlim(-.6, 1.6); ax.set_ylabel(unit, color=INK2, fontsize=9)
        ax.set_title(label.replace(' (Hz)', ''), fontsize=9, color=INK, loc='left', pad=26)
        line = 'mean' if s.normal == 'yes' else 'median'
        ax.text(0, 1.02, f"{s.centre_test}: p = {s.p_centre:.3f}  ·  {s.spread_test}: p = {s.p_spread:.3f}\n"
                         f"black line = {line}; Shapiro-Wilk p = {s.shapiro_p_ctrl:.2f} / {s.shapiro_p_zic2:.2f}",
                transform=ax.transAxes, ha='left', va='bottom', fontsize=7, color=MUTED)
        style(ax)
    for u in dict.fromkeys(u for _, _, u in METRICS):   # share y within a unit [v2.3]
        group = [ax for ax, (_, _, mu) in zip(axes.flat, METRICS) if mu == u]
        lo = min(ax.get_ylim()[0] for ax in group)
        hi = max(ax.get_ylim()[1] for ax in group)
        for ax in group: ax.set_ylim(lo, hi)
    for ax in axes.flat[len(METRICS):]: ax.set_visible(False)
    plt.tight_layout(); save_figure(path); plt.close()

# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    os.makedirs(OUT_FOLDER, exist_ok=True)
    run_started = datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %Z')
    print(f'{SCRIPT_NAME} v{VERSION} ({VERSION_DATE}) - run started {run_started}')
    fps = load_frame_rates(FRAME_RATE_FILE)
    video_rates = load_video_rates(VIDEO_RATE_FILE)
    print(f'Per-video frame rates: {len(video_rates)} videos from {VIDEO_RATE_FILE}'
          if video_rates else 'Per-video frame rates: none found, using one rate per mouse')
    animals, trials, spectra = [], [], []
    for grp, _ in GROUPS:
        gdir = os.path.join(DATA_FOLDER, grp)
        for f in sorted(x for x in os.listdir(gdir) if x.endswith('.xlsx')):
            name = f[:-5]
            if name in EXCLUDED_ANIMALS:
                print(f'[{grp:7}] {name}: excluded (EXCLUDED_ANIMALS)'); continue
            if name not in fps:
                print(f'WARNING: {name} has no frame rate in {FRAME_RATE_FILE} - skipped'); continue
            m, tr, sp = process_animal(os.path.join(gdir, f), name, grp, fps[name], video_rates)
            animals.append(m); trials.append(tr); spectra.append((grp, sp))
            print(f"[{grp:7}] {name}: {m['n_segments_used']:2d} segments, "
                  f"{m['seconds_used']:5.1f} s, dominant Y {m['dom_Hz_Y_nonpref']:.2f} / {m['dom_Hz_Y_pref']:.2f} Hz, "
                  f"X {m['dom_Hz_X_nonpref']:.2f} / {m['dom_Hz_X_pref']:.2f} Hz (non-pref / pref)")
    df = pd.DataFrame(animals); df.insert(0, 'script_version', VERSION)
    df.to_csv(os.path.join(OUT_FOLDER, 'summary_per_animal.csv'), index=False)
    seg = pd.concat(trials); seg.insert(0, 'script_version', VERSION)
    seg.to_csv(os.path.join(OUT_FOLDER, 'per_segment.csv'), index=False)
    st = group_stats(df); st.insert(0, 'script_version', VERSION)
    st.to_csv(os.path.join(OUT_FOLDER, 'group_statistics.csv'), index=False, float_format='%.4f')
    with open(os.path.join(OUT_FOLDER, 'run_info.txt'), 'w') as fh:
        fh.write(f'Script        : {os.path.basename(__file__)}\n')
        fh.write(f'Version       : {VERSION} ({VERSION_DATE})\n')
        fh.write(f'Run started   : {run_started}\n')
        fh.write(f'Python        : {platform.python_version()}  numpy {np.__version__}  scipy {scipy.__version__}  pandas {pd.__version__}\n')
        fh.write(f'Mice analysed : {len(df)} ({(df.group == GROUPS[0][0]).sum()} {GROUPS[0][0]}, {(df.group == GROUPS[1][0]).sum()} {GROUPS[1][0]})\n')
        fh.write(f'Excluded      : {", ".join(EXCLUDED_ANIMALS)}\n')
        fh.write(f'Frame rates   : {os.path.relpath(FRAME_RATE_FILE, PROJECT)}' + (f' + {len(video_rates)} per-video rates' if video_rates else '') + '\n')
        for k in ['STOP_MIN_S', 'STOP_PERCENTILE', 'MIN_SEGMENT_S', 'PSD_WIN_S', 'HIGHPASS_HZ', 'FREQ_MIN',
                  'FREQ_MAX', 'LOWPASS_HZ', 'CYCLE_PROM_SD', 'ENTROPY_FS', 'ENTROPY_M', 'ENTROPY_R',
                  'NORMALITY_ALPHA', 'T_TEST_EQUAL_VAR', 'N_BOOT']:
            fh.write(f'{k:<14}: {globals()[k]}\n')
    plot_spectra(spectra, os.path.join(OUT_FOLDER, 'group_spectra.png'))
    plot_metrics(df, st, os.path.join(OUT_FOLDER, 'group_metrics.png'))
    pd.set_option('display.width', 200)
    show = st[['metric', 'n_ctrl', 'ctrl_mean', 'ctrl_median', 'n_zic2', 'zic2_mean', 'zic2_median',
               'shapiro_p_ctrl', 'shapiro_p_zic2', 'normal', 'centre_test', 'p_centre',
               'spread_test', 'p_spread', 'effect_size_type', 'effect_size']]
    print(show.round(3).to_string(index=False))
    print(f'\nOutputs written to {OUT_FOLDER} (see run_info.txt)')
