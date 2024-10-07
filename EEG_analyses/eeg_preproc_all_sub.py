# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 10:35:02 2024

@author: Audrey Mazancieux

Preprocessing for EEG data of the frequency tagging task.
"""

import os
import glob
import mne
import numpy as np
import pickle
from mne.channels import make_standard_montage
from pyprep.find_noisy_channels import NoisyChannels

# =============================================================================

# Import EEG data
ROOT_DIR = "C:/Users/Admin/Desktop/RESEARCH PROJECTS ANALYSES/freq_tag_consciousness"
EEG_DIR = 'EEG_analyses'
DATA_DIR = 'Data'
RESULT_DIR = 'Results'
BEHAV_DIR = 'Behaviour'

SUBJECTS = [3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51]
N_STIM_PER_SEQ = 240
RESAMPLE_FREQ = 250

info_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}


# =============================================================================

## Preprocess data for all subjects

for subject in SUBJECTS : 
    
    # import EEG data
    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*.bdf"))[0]
    raw_data = mne.io.read_raw_bdf(data_file, stim_channel="Status", preload=True)
    print(raw_data.info)
    
    # downsampling
    events = mne.find_events(raw_data, shortest_event = 1)
    info_all_sub[f'sub-{subject}']['events'] = events
    raw_data, events = raw_data.resample(RESAMPLE_FREQ, events=events)
    
    # set the 32 system BioSemi channel positions on the data
    raw_data.pick_types(meg=False, eeg=True)
    montage = make_standard_montage('biosemi64')
    new_names = dict(zip(raw_data.ch_names, montage.ch_names))
    raw_data.rename_channels(new_names)
    raw_data.set_montage(montage)
    
    # apply bandpass filter
    raw_data.filter(l_freq=0.1, h_freq=30, fir_design="firwin", verbose=False)
        
    # get bad electrods and interpolate    
    noisy= NoisyChannels(raw_data)
    noisy.find_all_bads()
    bad_chs = noisy.get_bads()
    info_all_sub[f'sub-{subject}']['bad_chs'] = bad_chs
    raw_data.info['bads'] = bad_chs
    raw_data = raw_data.interpolate_bads()
        
    # perform Independent Component Analysis (ICA) on the copied data
    ica = mne.preprocessing.ICA(random_state=21, max_iter='auto')
    ica.fit(raw_data)

    # set the channel types for EOG channels
    raw_data.set_channel_types({'Fp1': 'eog', 'Fp2': 'eog'})

    # find Independent Components that match the EOG pattern
    eog_indices, eog_scores = ica.find_bads_eog(raw_data, measure='correlation', threshold=0.5) 
    ica.exclude = eog_indices
    info_all_sub[f'sub-{subject}']['bad_ica'] = eog_indices

    # set the channel types back to EEG channels
    raw_data.set_channel_types({'Fp1': 'eeg', 'Fp2': 'eeg'})

    # if there are EOG components to exclude, plot diagnostics and save the figures
    if len(eog_indices) > 0:
        ica.apply(raw_data)
        fig = ica.plot_components(picks=np.arange(ica.n_components_), show = False) # create figures of the components with the one removed
        fig.savefig(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f'sub-{subject}_ICA.png'))

    del ica
    
    # re-referencing to a robust average
    raw_data.set_eeg_reference('average') 

    # save the preprocessed raw data
    raw_data.save(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f'sub-{subject}_FreqTag_preproc.fif'), overwrite=True)
    
    
## Save output

with open(os.path.join(ROOT_DIR, EEG_DIR, 'Results', 'eeg_info_all_subjects.pickle'), 'wb') as f:
    pickle.dump(info_all_sub, f)     
    
