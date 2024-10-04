# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 10:35:02 2024

@author: Audrey Mazancieux

Preprocessing and first analyses for EEG data of the frequency tagging task.
"""

import os
import glob
import mne
import numpy as np
import pandas as pd
import seaborn as sns
import pickle
from mne.channels import make_standard_montage
from meegkit import ress
import scipy.signal as ss
from meegkit.utils import snr_spectrum
from pyprep.find_noisy_channels import NoisyChannels


# Define function
def find_nearest_index(array, target_value):
    '''Function which selects the index of value in a numpy array that is the 
    nearest to the target value
    '''
    array = np.asarray(array)
    return (np.abs(array - target_value)).argmin()


# =============================================================================

# Import EEG data
ROOT_DIR = "C:/Users/Admin/Desktop/RESEARCH PROJECTS ANALYSES/freq_tag_consciousness"
EEG_DIR = 'EEG_analyses'
DATA_DIR = 'Data'
RESULT_DIR = 'Results'
BEHAV_DIR = 'Behaviour'

SUBJECTS = [3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50]
N_STIM_PER_SEQ = 240
RESAMPLE_FREQ = 250
FACE_FREQ = 1.2
IMAGE_FREQ = 6

bads_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}

snr_contrast_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_pas_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_acc_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_conf_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}

epochs_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}


for subject in SUBJECTS : 
    
    # =========================================================================
    
    ## Get data 
    
    # import EEG data
    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*.bdf"))[0]
    raw_data = mne.io.read_raw_bdf(data_file, stim_channel="Status", preload=True)
    print(raw_data.info)
    
    # downsampling
    events = mne.find_events(raw_data, shortest_event = 1)
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
    bads_all_sub[f'sub-{subject}'] = bad_chs
    raw_data.info['bads'] = bad_chs
    raw_data = raw_data.interpolate_bads()
                                
    # re-referencing to a robust average
    raw_data.set_eeg_reference('average') 
    
    # epochs entire sequences (excluding fade periods)
    epochs = mne.Epochs(raw_data, events, event_id=10, tmin=1.667, tmax=41.667, baseline=None)
    epoch_data = epochs.get_data()
    
    n_trial = len(epoch_data)
     
    # =========================================================================
    # Frequency domain: data extraction 
    # =========================================================================
    
    ## Estimate SNR for each contrast
    
    # get behaviour
    behav_file = glob.glob(os.path.join(ROOT_DIR, BEHAV_DIR, 'Data', f'sub-{subject}', 'Freq*preproc.csv'))[0]
    behav = pd.read_csv(behav_file)
    
    contrasts = behav['contrast'].unique()
    snr_contrast_all_sub[f'sub-{subject}'] =  {contrast: [] for contrast in contrasts}
    
    for contrast in contrasts:
        
        idx = np.where(behav['contrast'] == contrast)[0]
        epoch_data_pas = epoch_data[idx, :, :].mean(axis=0)
        bins, psd = ss.welch(np.squeeze(epoch_data_pas.T), RESAMPLE_FREQ, nperseg=len(epoch_data_pas.T), axis=0)
        
        snr = snr_spectrum(psd, bins, skipbins=2, n_avg=10)
        snr_contrast_all_sub[f'sub-{subject}'][contrast].append(snr)
    
    
    ## Get indexes for frequencies of interest

    # faces
    idx_1_2hz = find_nearest_index(bins, FACE_FREQ)
    idx_2_4hz = find_nearest_index(bins, FACE_FREQ*2) 
    idx_3_6hz = find_nearest_index(bins, FACE_FREQ*3) 
    
    # images
    idx_6hz = find_nearest_index(bins, FACE_FREQ)
    idx_12hz = find_nearest_index(bins, FACE_FREQ*2) 
    idx_18hz = find_nearest_index(bins, FACE_FREQ*3) 
                                          

    ## Estimate SNR for mean epochs per PAS
    
    # Get snr for each PAS and of each contrast
    snr_pas_all_sub[f'sub-{subject}'] = {contrast: [] for contrast in contrasts}
    
    for contrast in contrasts:
   
        snr_pas_all_sub[f'sub-{subject}'][contrast] = {pas_num: [] for pas_num in range(1, 5)}
        for pas_num in range(1, 5):
            
            # get sequence indexes
            idx_pas = np.where((behav['contrast'] == contrast) & (behav['pas_score'] == pas_num))[0]

            if len(idx_pas) != 0:
                
                # only keep relevant epochs
                epoch_data_pas = epoch_data[idx_pas, :, :].mean(axis=0)
                
                # get psd 
                bins, psd = ss.welch(np.squeeze(epoch_data_pas.T), RESAMPLE_FREQ, nperseg=len(epoch_data_pas.T), axis=0)
                
                # get snr
                snr_pas = snr_spectrum(psd, bins, skipbins=2, n_avg=10)
                snr_pas_all_sub[f'sub-{subject}'][contrast][pas_num].append(snr_pas)  
                
    
    ## Estimate SNR for mean epochs per accuracy
    
    accuracy = behav['accuracy'].unique()
    snr_acc_all_sub[f'sub-{subject}'] = {contrast: [] for contrast in contrasts}
    
    for contrast in contrasts:
   
        snr_acc_all_sub[f'sub-{subject}'][contrast] = {acc: [] for acc in accuracy}
        for acc in accuracy:
            
            # get sequence indexes
            idx_acc = np.where((behav['contrast'] == contrast) & (behav['accuracy'] == acc))[0]
                    
            if len(idx_acc) != 0:
                
                # only keep relevant epochs
                epoch_data_acc = epoch_data[idx_acc, :, :].mean(axis=0)
                
                # get psd and snr
                bins, psd = ss.welch(np.squeeze(epoch_data_acc.T), RESAMPLE_FREQ, nperseg=len(epoch_data_pas.T), axis=0)
                
                # get snr
                snr_pas = snr_spectrum(psd, bins, skipbins=2, n_avg=10)
                snr_acc_all_sub[f'sub-{subject}'][contrast][acc].append(snr_pas) 
                
                     
    ## Estimate SNR for mean epochs per confidence rating
            
    snr_conf_all_sub[f'sub-{subject}'] = {contrast: [] for contrast in contrasts}
    
    for contrast in contrasts:
   
        snr_conf_all_sub[f'sub-{subject}'][contrast] = {conf: [] for conf in range(1, 5)}
        for conf in range(1, 5):
            
            # get sequence indexes
            idx_conf = np.where((behav['contrast'] == contrast) & (behav['conf_score'] == conf))[0]
                    
            if len(idx_conf) != 0:
                
                # only keep relevant epochs
                epoch_data_conf = epoch_data[idx_conf, :, :].mean(axis=0)
                
                # get psd 
                bins, psd = ss.welch(np.squeeze(epoch_data_conf.T), RESAMPLE_FREQ, nperseg=len(epoch_data_pas.T), axis=0)
                
                # get snr
                snr_pas = snr_spectrum(psd, bins, skipbins=2, n_avg=10)
                snr_conf_all_sub[f'sub-{subject}'][contrast][conf].append(snr_pas) 
    
    # =========================================================================
    # Time domain: data extraction 
    # =========================================================================
    
    # get trigger ids
    event_id = np.unique(events[:,2])
    event_id = list(event_id[1:-1])
    
    # get epochs data         
    channels = np.asarray(raw_data.info['ch_names'])[0:64, ]
    epochs_behav = mne.Epochs(raw_data, events, event_id = event_id, picks=channels , tmin=-0.167, tmax=0.83333, baseline=(-0.167, 0))
        
    
    ## Get data and save 
    
    epochs_behav_data = epochs_behav.get_data()   
    n_epoch_per_seq = N_STIM_PER_SEQ + 3 # 3 is the number of questions (pas, gender categorisation, confidence)    
    n_seq = len(epochs_behav_data) / n_epoch_per_seq # check number of sequences
    print(f'Number of sequences = {n_seq}')
    
    # only keep relevany events
    events_clean = pd.DataFrame({'events': events[:,2]})
    events_clean = events_clean.drop(events_clean[events_clean.events == 10].index)
    events_clean = events_clean.drop(events_clean[events_clean.events == 65536].index)
        
    # get epoch per sequences and average
    epoch_data_all_seq = {'1%': {},
                          '1.5%': {}}
    for seq in range(0, int(n_seq)): 
        
        epoch_data_all_seq['1%'] = {'nonfaces': [], 'faces': [], 'responses': []}
        epoch_data_all_seq['1.5%'] = {'nonfaces': [], 'faces': [], 'responses': []}
        
        epoch_data_seq = epochs_behav_data[(seq*n_epoch_per_seq):(seq*n_epoch_per_seq + n_epoch_per_seq), :, :]
        event_seq = np.array(events_clean['events'])[(seq*n_epoch_per_seq):(seq*n_epoch_per_seq + n_epoch_per_seq)]
        nonfaces_idx = np.where(event_seq == 20)
        faces_idx = np.where(event_seq == 21)
        
        epoch_nonfaces = epoch_data_seq[nonfaces_idx, :, :].mean(axis=0)
        epoch_faces = epoch_data_seq[faces_idx, :, :].mean(axis=0)
      
        if np.array(behav['contrast'])[seq] == '1%':
            epoch_data_all_seq['1%']['nonfaces'].append(epoch_nonfaces.mean(axis=0))
            epoch_data_all_seq['1%']['faces'].append(epoch_faces.mean(axis=0))
            epoch_data_all_seq['1%']['responses'].append(epoch_data_seq[-3:])
        else:
            epoch_data_all_seq['1.5%']['nonfaces'].append(epoch_nonfaces.mean(axis=0))
            epoch_data_all_seq['1.5%']['faces'].append(epoch_faces.mean(axis=0))
            epoch_data_all_seq['1.5%']['responses'].append(epoch_data_seq[-3:])
                   
    # create data and info dictionnary 
    data_sub_dict = {'data': epoch_data_all_seq,
                     'channels': raw_data.info['ch_names'],
                     'time': epochs_behav.times}
    # save
    epochs_all_sub[f'sub-{subject}'] = data_sub_dict


# =============================================================================    
   
## Save outputs

bads_all_sub['channels'] = raw_data.info['ch_names']
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'bads_all_subjects.pickle'), 'wb') as f:
    pickle.dump(bads_all_sub, f) 

with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_contrast_all_subjects.pickle'), 'wb') as f:
    pickle.dump(snr_contrast_all_sub, f)

with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_pas_all_subjects.pickle'), 'wb') as f:
    pickle.dump(snr_pas_all_sub, f)
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_acc_all_subjects.pickle'), 'wb') as f:
    pickle.dump(snr_acc_all_sub, f)
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_conf_all_subjects.pickle'), 'wb') as f:
    pickle.dump(snr_conf_all_sub, f)
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'epochs_all_subjects.pickle'), 'wb') as f:
    pickle.dump(epochs_all_sub, f)
    
    