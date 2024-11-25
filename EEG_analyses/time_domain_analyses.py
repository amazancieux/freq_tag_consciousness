# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:58:09 2024

@author: Audrey Mazancieux

Time-domain analyses: epoch data according to contrast and behavioural data.
Plot and perform permutation cluster analyses. 
"""

import os
import glob
import mne
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function

def trigger_in_sequence(event, events, sequence_trigger):
    """
    Check if a given event is within a given sequence.
    
    Parameters:
    -----------
    event : array
        The specific event to check, e.g., [sample_index, 0, trigger_value].
    events : array
        Array of all events, where each row is [sample_index, 0, trigger_value].
    sequence_trigger : int
        The trigger value for the sequence to match. 
        
    """
    # extract event times and triggers
    event_samples = events[:, 0]
    event_triggers = events[:, 2]
    
    # find all sequence triggers
    sequence_indices = np.where(event_triggers == sequence_trigger)[0]
    
    # loop through sequence triggers to find the correct window
    for i, seq_idx in enumerate(sequence_indices):
        seq_sample = event_samples[seq_idx]
        
        # define end of the sequence as the next sequence trigger or the end of events
        if i + 1 < len(sequence_indices):
            next_seq_sample = event_samples[sequence_indices[i + 1]]
        else:
            next_seq_sample = np.inf  # no next sequence; assume sequence extends to infinity
        
        # check if the event is within this sequence window
        if seq_sample <= event[0] < next_seq_sample:
            result = True
        else:
            result = False
                     
    return result  

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
IMAGE_FREQ = 6


# =============================================================================

for subject in SUBJECTS : 

    # import preprocessed data 
    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*preproc.fif"))[0]
    data = mne.io.read_raw_fif(data_file, preload=True)
      
    # use notch filter for images frequencies and harmonics 
    freqs_to_remove = [IMAGE_FREQ, IMAGE_FREQ*2, IMAGE_FREQ*3, IMAGE_FREQ*4, IMAGE_FREQ*5, IMAGE_FREQ*6]
    data.notch_filter(freqs=freqs_to_remove, method='spectrum_fit') 
        
    # get trigger and only keep start seq trigger
    events = pd.DataFrame(mne.find_events(data, shortest_event = 1))
    event10 = events[(events[2] == 10)]
    
    # get contrast from behaviour data
    behav_file = glob.glob(os.path.join(ROOT_DIR, BEHAV_DIR, 'Data', f'sub-{subject}', 'Freq*preproc.csv'))[0]
    behav = pd.read_csv(behav_file)   
        
    # add trigger 11 for contrast 1.5% and 10 for 1%
    event10['contrast'] = np.array(behav['contrast'])
    event10[2] = np.where(event10['contrast'] == '1%', 10, 11)
    del event10['contrast']
    
    # merge with other trigger
    events_new = events.merge(event10, on=0, how='left', suffixes=('_df1', '_df2'))
    events_new['new_trigger'] = events_new['2_df2'].combine_first(events_new['2_df1'])
    events[2] = events_new['new_trigger']
    
    # change trigger within sequences if sequence is 11 (contrast 1.5%)
    events = np.array(events)
    events_updated = events.copy()

    for i, event in enumerate(events):
        if trigger_in_sequence(event, events, sequence_trigger=11):
            if event[2] == 20: 
                events_updated[i, 2] = 201  
            if event[2] == 21: 
                events_updated[i, 2] = 211 
            if event[2] == 31: 
                events_updated[i, 2] = 311  
            if event[2] == 32: 
                events_updated[i, 2] = 321 
            if event[2] == 33: 
                events_updated[i, 2] = 331  
            if event[2] == 34: 
                events_updated[i, 2] = 341 
            if event[2] == 40: 
                events_updated[i, 2] = 401  
            if event[2] == 41: 
                events_updated[i, 2] = 411 
            if event[2] == 51: 
                events_updated[i, 2] = 511  
            if event[2] == 52: 
                events_updated[i, 2] = 521 
            if event[2] == 53: 
                events_updated[i, 2] = 531  
            if event[2] == 54: 
                events_updated[i, 2] = 541 
       
    # get unique trigger
    event_id = np.unique(events_updated)
    event_id = event_id[np.where((event_id >= 20) & (event_id <= 541))]
    event_id = [int(x) for x in event_id]
    
    # epoch data 
    events_updated = events_updated.astype(int)
    epochs = mne.Epochs(data, events_updated, event_id=event_id, tmin=-0.167, tmax=0.833, baseline=(-0.167, 0))  
    
    # save
    epochs.save(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f'epochs_clean_sub-{subject}.fif'))
        
        
# =============================================================================

## Define ROIs

# from Quek & de Heering (2024)
ROI_OT_1 = ['O1', 'PO3', 'PO7', 'P7', 'P9'] # face left
ROI_OT_2 = ['O2', 'PO4', 'PO8', 'P8', 'P10'] # face right
roi_face = ROI_OT_1 + ROI_OT_2

# get epochs
epoch_files = [glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f'epochs_clean_sub-{subject}.fif')) for subject in SUBJECTS]
epoch_all = [mne.read_epochs(x[0]).pick(roi_face).apply_baseline((-0.167, 0)) for x in epoch_files]


## Grand average per stimuli

# trigger to analyses
conditions = [20, 201, 21, 211]

evokeds_stim_all = {}
for c in conditions:
    # average evoked responses for each participant
    evokeds_stim_all[c] = [
        sub[c].average()
        for sub in epoch_all]
    
# rename condition
evokeds_stim_all['non_face_1%'] = evokeds_stim_all[20]
evokeds_stim_all['non_face_1.5%'] = evokeds_stim_all[201]
evokeds_stim_all['face_1%'] = evokeds_stim_all[21]
evokeds_stim_all['face_1.5%'] = evokeds_stim_all[211]

del evokeds_stim_all[20]
del evokeds_stim_all[201]
del evokeds_stim_all[21]
del evokeds_stim_all[211]

# plot grand average per stimuli
mne.viz.plot_compare_evokeds(evokeds_stim_all,
                            linestyles=['solid', 'solid', 'solid', 'solid'],  
                             styles={'non_face_1%': {"linewidth": 4},  
                                     'non_face_1.5%': {"linewidth": 4},  
                                     'face_1%': {"linewidth": 4}, 
                                     'face_1.5%': {"linewidth": 4},  
                                     },      
                             colors=["#f26d41","#40b2d8", "#D84BBE", "#ac4bd8"])





