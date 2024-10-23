# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 10:35:02 2024

@author: Audrey Mazancieux

Preprocessing for EEG data of the frequency tagging task.
"""

import os
import glob
import mne
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

SUBJECTS = [3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52]
N_STIM_PER_SEQ = 240
RESAMPLE_FREQ = 250

bads_electrods = {f'sub-{sub}': {} for sub in SUBJECTS}


# =============================================================================

## Preprocess data for all subjects

for subject in SUBJECTS : 
    
    # import EEG data
    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*.bdf"))[0]
    raw_data = mne.io.read_raw_bdf(data_file, stim_channel="Status", preload=True)
    print(raw_data.info)
    
    # downsampling
    events = mne.find_events(raw_data, shortest_event = 1)
    raw_data, events = raw_data.resample(RESAMPLE_FREQ, events=events)
    
    # set the 32 system BioSemi channel positions on the data
    #raw_data.pick_types(meg=False, eeg=True)
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
    bads_electrods[f'sub-{subject}'] = bad_chs
    raw_data.info['bads'] = bad_chs
    raw_data = raw_data.interpolate_bads()
           
    # re-referencing to a robust average
    raw_data.set_eeg_reference('average') 

    # save the preprocessed raw data
    raw_data.save(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f'sub-{subject}_FreqTag_preproc.fif'), overwrite=True)
    
    
## Save output

with open(os.path.join(ROOT_DIR, EEG_DIR, 'Results', 'bad_electodes_all_subjects.pickle'), 'wb') as f:
    pickle.dump(bads_electrods, f)     
    
