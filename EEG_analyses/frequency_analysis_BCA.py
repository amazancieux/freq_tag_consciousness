# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:47:16 2026

@author: Audrey Mazancieux

Frequency analyses based on baseline-corrected amplitudes (BCA).
"""

import os
import glob
import pickle
import mne
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mne.viz import plot_topomap

# =============================================================================

# Define parameters

ROOT_DIR = "C:/Users/Admin/Desktop/RESEARCH PROJECTS ANALYSES/freq_tag_consciousness"
EEG_DIR = 'EEG_analyses'
DATA_DIR = 'Data'
RESULT_DIR = 'Results_BCA'
BEHAV_DIR = 'Behaviour'
RESULT_PATH = os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR)

SUBJECTS = [3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            34, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52]

RESAMPLE_FREQ = 250      
FACE_FREQ = 1.2
IMAGE_FREQ = 6.0

# Epoch = whole sequence without the fade-in / fade-out periods
EVENT_ID = 10
EPOCH_TMIN = 1.667        # s
EPOCH_DURATION = 40.0     # s (must be an integer number of cycles of both frequencies)

# Baseline correction
N_NEIGHBORS = 10          # bins taken on each side
N_SKIP = 1                # immediately adjacent bins skipped on each side
EXCLUDE_EXTREMES = True   # drop local min and max among the neighbouring bins

# Harmonic selection
FACE_MAX_FREQ = 30.0      # highest candidate face harmonic (Hz)
IMAGE_MAX_FREQ = 48.0     # highest candidate base harmonic (Hz)
Z_THRESH = 1.64           
N_NONSIG_STOP = 2         # stop after this many consecutive non-significant harmonics

FACE_HARMONICS = None
IMAGE_HARMONICS = None

SPECTRUM_FMAX = 50.0      # spectra are kept up to this frequency (Hz)
MIN_TRIALS = 1            # conditions with fewer sequences are skipped

# ROIs from Quek & de Heering (2024)
ROIS = {'OCC': ['Oz', 'O1', 'O2'],                    # base
        'OT1': ['O1', 'PO3', 'PO7', 'P7', 'P9'],      # face, left
        'OT2': ['O2', 'PO4', 'PO8', 'P8', 'P10']}     # face, right
FACE_SELECTION_ROIS = ['OT1', 'OT2']   # ROIs used to select the face harmonics
IMAGE_SELECTION_ROIS = ['OCC']         # ROIs used to select the base harmonics


# =============================================================================

## Define functions

def amplitude_spectrum(data, sfreq):
    """Single-sided FFT amplitude spectrum along the last axis.
    Returns the frequencies and the amplitudes in the same unit as `data`
    """
    n_times = data.shape[-1]
    amp = np.abs(np.fft.rfft(data, axis=-1)) / n_times
    amp[..., 1:] *= 2
    if n_times % 2 == 0:
        amp[..., -1] /= 2      # Nyquist bin is not doubled
    freqs = np.fft.rfftfreq(n_times, d=1. / sfreq)
    return freqs, amp


def baseline_correct(amp, n_neighbors=N_NEIGHBORS, n_skip=N_SKIP,
                     exclude_extremes=EXCLUDE_EXTREMES):
    """Baseline-correct an amplitude spectrum (frequency = last axis).

    For each bin, the noise level is the mean of the `n_neighbors` bins on
    each side, skipping the `n_skip` closest bins on each side and, if
    `exclude_extremes`, the local minimum and maximum of those bins.

    Returns
    -------
    bca : amplitude - noise mean (baseline-corrected amplitude)
    snr : amplitude / noise mean
    z   : (amplitude - noise mean) / noise SD
    """
    amp = np.asarray(amp, dtype=float)
    n_freqs = amp.shape[-1]
    pad = n_skip + n_neighbors
    padded = np.pad(amp, [(0, 0)] * (amp.ndim - 1) + [(pad, pad)],
                    constant_values=np.nan)

    # offsets of the neighbouring bins, e.g. -11..-2 and +2..+11
    offsets = np.r_[-pad:-n_skip, n_skip + 1:pad + 1]
    idx = np.arange(n_freqs)[:, None] + pad + offsets[None, :]
    neigh = padded[..., idx]                     # (..., n_freqs, 2 * n_neighbors)

    if exclude_extremes:
        imax = np.nanargmax(neigh, axis=-1)[..., None]
        imin = np.nanargmin(neigh, axis=-1)[..., None]
        np.put_along_axis(neigh, imax, np.nan, axis=-1)
        np.put_along_axis(neigh, imin, np.nan, axis=-1)

    noise_mean = np.nanmean(neigh, axis=-1)
    noise_sd = np.nanstd(neigh, axis=-1, ddof=1)

    bca = amp - noise_mean
    snr = amp / noise_mean
    z = (amp - noise_mean) / noise_sd
    return bca, snr, z


def candidate_harmonics(f0, fmax, exclude_every=None):
    """Harmonics of f0 up to fmax. `exclude_every=5` removes every 5th harmonic
    (for f0 = 1.2 Hz: 6, 12, 18 Hz... which belong to the 6 Hz response)."""
    n = np.arange(1, int(np.floor(fmax / f0 + 1e-9)) + 1)
    if exclude_every:
        n = n[n % exclude_every != 0]
    return n * f0


def freq_to_bins(freqs, targets):
    """Exact FFT bin of each target frequency (error if not on a bin)."""
    df = freqs[1] - freqs[0]
    targets = np.asarray(targets, dtype=float)
    idx = np.round(targets / df).astype(int)
    if not np.allclose(freqs[idx], targets, atol=df / 100):
        raise ValueError('Some target frequencies do not fall on an FFT bin: '
                         'the epoch must contain an integer number of cycles.')
    return idx


def select_harmonics(z_values, z_thresh=Z_THRESH, n_nonsig_stop=N_NONSIG_STOP):
    """Keep all harmonics from the first one up to the last significant one,
    stopping once `n_nonsig_stop` consecutive harmonics are not significant.
    Returns a boolean mask over the candidate harmonics."""
    last_sig, n_nonsig = -1, 0
    for i, z in enumerate(z_values):
        if z > z_thresh:
            last_sig, n_nonsig = i, 0
        else:
            n_nonsig += 1
            if n_nonsig >= n_nonsig_stop:
                break
    mask = np.zeros(len(z_values), dtype=bool)
    mask[:last_sig + 1] = True
    return mask

def normalize_contrast(value):
    """Normalize a contrast label so that different spellings of the same value
    end up in the same condition ('1.50%', '1,5 %', 1.5 -> '1.5%').
 
    The numeric value is extracted and reformatted; anything that cannot be
    parsed as a number is returned unchanged (after stripping spaces).
    """
    text = str(value).strip()
    cleaned = text.replace('%', '').replace(',', '.').replace(' ', '')
    try:
        number = float(cleaned)
    except ValueError:
        return text        # e.g. a real label such as 'catch' or 'full'
    return f'{number:g}%'  # 1.5 -> '1.5%', 1.0 -> '1%'


def get_conditions(behav):
    """List of (analysis, contrast, level, sequence indices)."""
    conditions = []
    for contrast in sorted(behav['contrast'].unique()):
        c_mask = (behav['contrast'] == contrast).to_numpy()
        conditions.append(('contrast', contrast, 'all', np.where(c_mask)[0]))
        for column, analysis in [('pas_score', 'PAS'),
                                 ('accuracy', 'accuracy'),
                                 ('conf_score', 'confidence')]:
            for level in sorted(behav.loc[c_mask, column].dropna().unique()):
                idx = np.where(c_mask & (behav[column] == level).to_numpy())[0]
                conditions.append((analysis, contrast, level, idx))
    return conditions


def plot_topo_grid(results, analysis, key, info, title, fname):
    """Group-average topographies: one row per contrast, one column per level,
    shared colour scale."""
    res = [r for r in results if r['analysis'] == analysis]
    contrasts = sorted(set(r['contrast'] for r in res))
    levels = sorted(set(r['level'] for r in res))

    maps = {}
    for c in contrasts:
        for lvl in levels:
            vals = [r[key] for r in res if r['contrast'] == c and r['level'] == lvl]
            if vals:
                maps[(c, lvl)] = (np.mean(vals, axis=0), len(vals))

    vmin = min(m.min() for m, _ in maps.values())
    vmax = max(m.max() for m, _ in maps.values())

    fig, axes = plt.subplots(len(contrasts), len(levels), squeeze=False,
                             figsize=(2.6 * len(levels) + 1.2, 2.8 * len(contrasts)))
    im = None
    for i, c in enumerate(contrasts):
        for j, lvl in enumerate(levels):
            ax = axes[i, j]
            if (c, lvl) not in maps:
                ax.axis('off')
                continue
            m, n = maps[(c, lvl)]
            im, _ = plot_topomap(m, info, axes=ax, show=False, cmap='viridis',
                                 vlim=(vmin, vmax))
            lvl_txt = '' if analysis == 'contrast' else f' - {analysis} {lvl}'
            ax.set_title(f'{c}{lvl_txt}\n(n = {n})', fontsize=10)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, label='Summed BCA (µV)')
    fig.suptitle(title)
    fig.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close(fig)
 
    
# =============================================================================

## Define frequency analyses parameters

# define frequency bins 
n_times = int(round(EPOCH_DURATION * RESAMPLE_FREQ))
freqs_full = np.fft.rfftfreq(n_times, d=1. / RESAMPLE_FREQ)
keep = freqs_full <= SPECTRUM_FMAX
freqs = freqs_full[keep]
df = freqs[1] - freqs[0]

# define number of harmonics and bins
face_cand = candidate_harmonics(FACE_FREQ, FACE_MAX_FREQ,
                                exclude_every=int(round(IMAGE_FREQ / FACE_FREQ)))
image_cand = candidate_harmonics(IMAGE_FREQ, IMAGE_MAX_FREQ)
face_bins = freq_to_bins(freqs, face_cand)
image_bins = freq_to_bins(freqs, image_cand)


# =============================================================================

## Analyses for each subject 

results = []

for subject in SUBJECTS:

    print(f'sub-{subject}')

    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*preproc.fif"))[0]
    raw = mne.io.read_raw_fif(data_file, verbose='error')
    if raw.info['sfreq'] != RESAMPLE_FREQ:
        raise ValueError(f"sub-{subject}: sampling rate is {raw.info['sfreq']} Hz")

    events = mne.find_events(raw, shortest_event=1, verbose='error')
    epochs = mne.Epochs(raw, events, event_id=EVENT_ID, tmin=EPOCH_TMIN,
                        tmax=EPOCH_TMIN + EPOCH_DURATION, baseline=None,
                        picks='eeg', preload=True, verbose='error')

    X = epochs.get_data(units='uV')[..., :n_times]
    if X.shape[-1] != n_times:
        raise ValueError(f'sub-{subject}: epochs shorter than {EPOCH_DURATION} s')
        
    ch_names, topo_info = epochs.ch_names, epochs.info
    
    # -------------------------------------------------------------------------   

    # get behaviour
    behav_file = glob.glob(os.path.join(ROOT_DIR, BEHAV_DIR, 'Data', f'sub-{subject}', 'Freq*preproc.csv'))[0]
    behav = pd.read_csv(behav_file)

    # harmonize contrast labels ('1.50%' and '1.5%' are the same condition)
    raw_contrast = behav['contrast'].copy()
    behav['contrast'] = raw_contrast.map(normalize_contrast)
    relabelled = {str(old): new for old, new in zip(raw_contrast, behav['contrast']) if str(old) != new}

    for analysis, contrast, level, idx in get_conditions(behav):

        if len(idx) < MIN_TRIALS:
            continue

        evoked = X[idx].mean(axis=0)                    # time-domain average
        _, amp = amplitude_spectrum(evoked, RESAMPLE_FREQ)
        amp = amp[:, keep]
        bca, _, _ = baseline_correct(amp)

        res = {'subject': f'sub-{subject}', 'analysis': analysis,
               'contrast': contrast, 'level': level, 'n_trials': len(idx),
               'bca_face': bca[:, face_bins],           
               'bca_image': bca[:, image_bins]}
        if analysis == 'contrast':
            res['amp'] = amp.astype(np.float32)         
        results.append(res)

    
# =============================================================================
    
## Select harmonic on the grand average (all subjects, both contrasts)

roi_idx = {roi: [ch_names.index(ch) for ch in chans] for roi, chans in ROIS.items()}
ga_amp = np.mean([r['amp'] for r in results if r['analysis'] == 'contrast'], axis=0)

def choose_harmonics(roi_names, cand, bins, manual, label):
    chans = sorted(set(i for roi in roi_names for i in roi_idx[roi]))
    _, _, z = baseline_correct(ga_amp[chans].mean(axis=0))
    z = z[bins]
    if manual is None:
        mask = select_harmonics(z)
    else:
        mask = np.array([np.any(np.isclose(f, manual)) for f in cand])
    print(f'\n{label} harmonics (z on grand average, ROI {"+".join(roi_names)}):')
    for f, zz, m in zip(cand, z, mask):
        print(f'  {f:5.1f} Hz  z = {zz:6.2f}  {"selected" if m else ""}')
    return mask, z


face_mask, face_z = choose_harmonics(FACE_SELECTION_ROIS, face_cand, face_bins, FACE_HARMONICS, 'Face (1.2 Hz)')
image_mask, image_z = choose_harmonics(IMAGE_SELECTION_ROIS, image_cand, image_bins, IMAGE_HARMONICS, 'Image (6 Hz)')


# =============================================================================

## Sum BCA per subject, condition and ROI  

rows = []
for r in results:
    r['face_sum'] = r['bca_face'][:, face_mask].sum(axis=1)     
    r['image_sum'] = r['bca_image'][:, image_mask].sum(axis=1)
    for roi, idx in roi_idx.items():
        rows.append({'Subject': r['subject'],
                     'Analysis': r['analysis'],
                     'Contrast': r['contrast'],
                     'Level': r['level'],
                     'N_trials': r['n_trials'],
                     'ROI': roi,
                     'BCA_1_2Hz': r['face_sum'][idx].mean(),
                     'BCA_6Hz': r['image_sum'][idx].mean()})

bca_df = pd.DataFrame(rows)
bca_df.to_csv(os.path.join(RESULT_PATH, 'BCA_summed_all_sub.csv'), index=False)

with open(os.path.join(RESULT_PATH, 'bca_all_sub.p'), 'wb') as f:
    pickle.dump({'results': results, 'freqs': freqs, 'ch_names': ch_names,
                 'face_harmonics': face_cand, 'face_mask': face_mask, 'face_z': face_z,
                 'image_harmonics': image_cand, 'image_mask': image_mask,
                 'image_z': image_z, 'rois': ROIS}, f)

# =============================================================================

## Plots

# grand-average BCA spectrum per ROI and contrast
contrasts_all = sorted(set(r['contrast'] for r in results))

fig, axes = plt.subplots(len(ROIS), 1, figsize=(11, 2.8 * len(ROIS)), sharex=True)
for ax, (roi, idx) in zip(np.atleast_1d(axes), roi_idx.items()):
    for f0 in face_cand[face_mask]:
        ax.axvline(f0, color='tab:red', alpha=0.15, lw=1)
    for f0 in image_cand[image_mask]:
        ax.axvline(f0, color='grey', alpha=0.3, lw=1)
    for c in contrasts_all:
        spectra = [baseline_correct(r['amp'][idx].mean(axis=0))[0]
                   for r in results if r['analysis'] == 'contrast' and r['contrast'] == c]
        ax.plot(freqs, np.mean(spectra, axis=0), lw=0.8, label=c)
    ax.set_title(f'ROI {roi}')
    ax.set_ylabel('BCA (µV)')
ax.set_xlabel('Frequency (Hz)')
ax.set_xlim(0, SPECTRUM_FMAX)
np.atleast_1d(axes)[0].legend(title='Contrast')
fig.tight_layout()
fig.savefig(os.path.join(RESULT_PATH, 'BCA_spectra_grand_average.png'), dpi=300)
plt.close(fig)

# topographies of summed BCA
for analysis in ['contrast', 'PAS', 'accuracy', 'confidence']:
    for key, label, tag in [('face_sum', 'Faces (1.2 Hz harmonics)', 'faces'),
                            ('image_sum', 'Images (6 Hz harmonics)', 'images')]:
        plot_topo_grid(results, analysis, key, topo_info,
                       f'Summed BCA - {label}',
                       os.path.join(RESULT_PATH, f'Topomaps_BCA_{tag}_{analysis}.png'))

