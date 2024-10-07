# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:21:21 2024

@author: Audrey Mazancieux

Frequency analyses: extract SNR at given frequencies and plots.
"""

import os
import mne
import glob
import numpy as np
import pandas as pd
import pickle
from mne.viz import plot_topomap
from matplotlib import pyplot as plt
import scipy.signal as ss
from meegkit.utils import snr_spectrum

# define function
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

SUBJECTS = [3, 14, 15, 17, 18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 37, 38, 39, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51]
N_STIM_PER_SEQ = 240
RESAMPLE_FREQ = 250
FACE_FREQ = 1.2
IMAGE_FREQ = 6

# =============================================================================

## Create snr dictionaries

snr_contrast_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_pas_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_acc_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}
snr_conf_all_sub = {f'sub-{sub}': {} for sub in SUBJECTS}


## Load infos

with open(os.path.join(ROOT_DIR, EEG_DIR, 'Results', 'eeg_info_all_subjects.pickle'), 'rb') as f:
    info_all_sub = pickle.load(f)    


for subject in SUBJECTS : 
    
    ## Epoch data
    
    # import preprocessed data 
    data_file = glob.glob(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, f'sub-{subject}', f"*{subject}*preproc.fif"))[0]
    data = mne.io.read_raw_fif(data_file)
            
    # epochs entire sequences (excluding fade periods)
    events = info_all_sub[f'sub-{subject}']['events']
    epochs = mne.Epochs(data, events, event_id=10, tmin=1.667, tmax=41.667, baseline=None)
    epoch_data = epochs.get_data()
    
    n_trial = len(epoch_data)
      
    
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
  
        
# =============================================================================

## Plot topographies 

# get faces indices
idx_1_2hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], FACE_FREQ)
idx_2_4hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], FACE_FREQ*2) 
idx_3_6hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], FACE_FREQ*3) 

# get images indices
idx_6hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], IMAGE_FREQ)
idx_12hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], IMAGE_FREQ*2) 
idx_18hz = find_nearest_index(snr_contrast_all_sub['frequency_bins'], IMAGE_FREQ*3) 

# get snr for each contrast and frequency
snr_1_2_C1 = []
snr_1_2_C2 = []
snr_6_C1 = []
snr_6_C2 = []

for sub in snr_pas_all_sub.keys():
    data = [snr_contrast_all_sub[sub][contrast][0][idx_1_2hz, :] for contrast in snr_contrast_all_sub[sub].keys()]
    snr_1_2_C1.append(data[0])
    snr_1_2_C2.append(data[1])
    data = [snr_contrast_all_sub[sub][contrast][0][idx_2_4hz, :] for contrast in snr_contrast_all_sub[sub].keys()]
    snr_6_C1.append(data[0])
    snr_6_C2.append(data[1])

snr_1_2_C1 = np.mean(snr_1_2_C1, axis=0)   
snr_1_2_C2 = np.mean(snr_1_2_C2, axis=0)   
snr_6_C1 = np.mean(snr_6_C1, axis=0)   
snr_6_C2 = np.mean(snr_6_C2, axis=0)    
 
# face topography 1.2 Hz at 1%
fig,(ax) = plt.subplots(figsize=(10, 10))
im, cm = plot_topomap(snr_1_2_C1, epochs.info, axes=ax)
ax_x_start = 0.92
ax_x_width = 0.02
ax_y_start = 0.25
ax_y_height = 0.45
cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
clb = fig.colorbar(im, cax=cbar_ax)
clb.ax.set_title('SNR', fontsize=18)
ax.title.set_text('SNR at 1.2Hz for 1% of contrast')
fig.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "Topo_1.2_HZ_SNR_1%.png"))

# face topography 1.2 Hz at 1.5%
fig,(ax) = plt.subplots(figsize=(10, 10))
im, cm = plot_topomap(snr_1_2_C2, epochs.info, axes=ax)
ax_x_start = 0.92
ax_x_width = 0.02
ax_y_start = 0.25
ax_y_height = 0.45
cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
clb = fig.colorbar(im, cax=cbar_ax)
clb.ax.set_title('SNR', fontsize=18)
ax.title.set_text('SNR at 1.2Hz for 1.5% of contrast')
fig.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "Topo_1.2_HZ_SNR_1.5%.png"))

# image topography 6 Hz at 1%
fig,(ax) = plt.subplots(figsize=(10, 10))
im, cm = plot_topomap(snr_6_C1, epochs.info, axes=ax)
ax_x_start = 0.92
ax_x_width = 0.02
ax_y_start = 0.25
ax_y_height = 0.45
cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
clb = fig.colorbar(im, cax=cbar_ax)
clb.ax.set_title('SNR', fontsize=18)
ax.title.set_text('SNR at 2.4Hz for 1% of contrast')
fig.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "Topo_2.4_HZ_SNR_1%.png"))

# image topography 6 Hz at 1.5%
fig,(ax) = plt.subplots(figsize=(10, 10))
im, cm = plot_topomap(snr_6_C2, epochs.info, axes=ax)
ax_x_start = 0.92
ax_x_width = 0.02
ax_y_start = 0.25
ax_y_height = 0.45
cbar_ax = fig.add_axes([ax_x_start, ax_y_start, ax_x_width, ax_y_height])
clb = fig.colorbar(im, cax=cbar_ax)
clb.ax.set_title('SNR', fontsize=18)
ax.title.set_text('SNR at 2.4Hz for 1.5% of contrast')
fig.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "Topo_2.4_HZ_SNR_1.5%.png"))
     

## Plot spectrograms for one participant

# get data for electrods for interest
SELECTED_ELECTRODS = ['Iz', 'Oz', 'O2', 'PO4', 'PO8', 'O1', 'PO3', 'PO7', 'P6', 'P8', 'P10', 'P5', 'P7', 'P9'] 
idx_ROIs = [data.info['ch_names'].index(electrod) for electrod in SELECTED_ELECTRODS]  

snr_C1 = []
snr_C2 = []
for sub in snr_pas_all_sub.keys():
    data = [snr_contrast_all_sub[sub][contrast][0][:, idx_ROIs] for contrast in snr_contrast_all_sub[sub].keys()]
    snr_C1.append(data[0])
    snr_C2.append(data[1])

snr_C1 = np.mean(snr_C1, axis=0)   
snr_C2 = np.mean(snr_C2, axis=0)   

# spectrogram for contrast 1%
f, ax = plt.subplots(figsize=(15, 8))
ax.plot(snr_contrast_all_sub['frequency_bins'], snr_C1.mean(axis=1), label="SNR")
ax.axvline(FACE_FREQ, ls=":", c="grey", zorder=0)
ax.axvline(FACE_FREQ*2, ls=":", c="grey", zorder=0)
ax.axvline(FACE_FREQ*3, ls=":", c="grey", zorder=0)
ax.axhline(1, ls=":", c="red", zorder=0)
ax.set_ylabel("SNR")
ax.set_xlabel("Frequency (Hz)")
ax.set_xlim([0, 20])
ax.set_ylim([0, 20])
ax.set_title('Spectrogram for te 1% contrast condition')
plt.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "spectrogram_1%.png"), dpi=300)

# spectrogram for contrast 1.5%
f, ax = plt.subplots(figsize=(15, 8))
ax.plot(snr_contrast_all_sub['frequency_bins'], snr_C2.mean(axis=1), label="SNR")
ax.axvline(FACE_FREQ, ls=":", c="grey", zorder=0)
ax.axvline(FACE_FREQ*2, ls=":", c="grey", zorder=0)
ax.axvline(FACE_FREQ*3, ls=":", c="grey", zorder=0)
ax.axhline(1, ls=":", c="red", zorder=0)
ax.set_ylabel("SNR")
ax.set_xlabel("Frequency (Hz)")
ax.set_xlim([0, 20])
ax.set_ylim([0, 20])
ax.set_title('Spectrogram for the 1.5% contrast condition')
plt.savefig(os.path.join(ROOT_DIR, EEG_DIR, RESULT_DIR, "spectrogram_1.5%.png"), dpi=300)


# =============================================================================
   
## Get SNR for each ROI

# from Quek & de Heering (2024)
ROI_OCC = ['Oz', 'O1', 'O2'] # for 6 Hz
ROI_OT_1 = ['O1', 'PO3', 'PO7', 'P7', 'P9'] # face left
ROI_OT_2 = ['O2', 'PO4', 'PO8', 'P8', 'P10'] # face right
                                                               
idx_ROI_OCC = [data.info['ch_names'].index(electrod) for electrod in ROI_OCC]
idx_ROI_OT_1 = [data.info['ch_names'].index(electrod) for electrod in ROI_OT_1]
idx_ROI_OT_2 = [data.info['ch_names'].index(electrod) for electrod in ROI_OT_2]

# idx_roi = [idx_ROI1, idx_ROI2, idx_ROI3, idx_ROI4, idx_ROI5]
idx_roi = [idx_ROI_OCC, idx_ROI_OT_1, idx_ROI_OT_2]
snr_roi_pas = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}


## Get snr per PAS at each frequency, subject, and ROI

for sub in snr_pas_all_sub.keys(): 
               
    for roi in range(len(idx_roi)): 
        contrasts = snr_pas_all_sub[sub].keys()
        
        for contrast in contrasts:    
               
            for pas_num in range(1, 5):
                
                if len(snr_pas_all_sub[sub][contrast][pas_num]) != 0:
                    
                    snr_1 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_1_2hz, idx_roi[roi]] 
                    snr_2 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_2_4hz, idx_roi[roi]]
                    snr_3 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_3_6hz, idx_roi[roi]]
                    snr_6 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_6hz, idx_roi[roi]]
                    snr_12 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_12hz, idx_roi[roi]]
                    snr_18 = snr_pas_all_sub[sub][contrast][pas_num][0][idx_18hz, idx_roi[roi]]
                     
                    new_data = [{'Subject': sub,
                                 'Contrast': contrast, 
                                 'PAS': pas_num,
                                 '1_2Hz': np.mean(snr_1, axis=0),
                                 '2_4Hz': np.mean(snr_2, axis=0),
                                 '3_6Hz': np.mean(snr_3, axis=0),
                                 '6Hz': np.mean(snr_6, axis=0),
                                 '12Hz': np.mean(snr_12, axis=0),
                                 '18Hz': np.mean(snr_18, axis=0)}]
                    
                    for row in new_data:
                        snr_roi_pas[f'roi{roi+1}'] = pd.concat([snr_roi_pas[f'roi{roi+1}'], pd.DataFrame(row, index=[0])], ignore_index=True)

# save data    
snr_roi_pas['roi1'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OCC_SNR_PAS_all_sub.csv'), index=False)  
snr_roi_pas['roi2'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT1_SNR_PAS_all_sub.csv'), index=False)  
snr_roi_pas['roi3'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT2_SNR_PAS_all_sub.csv'), index=False)  
                   

## Get snr per accuracy at each frequency, subject, and ROI

snr_roi_acc = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}
accuracy = snr_acc_all_sub['sub-3']['1%'].keys()
for sub in snr_pas_all_sub.keys(): 
               
    for roi in range(len(idx_roi)):   
        contrasts = snr_pas_all_sub[sub].keys()
        
        for contrast in contrasts:    
               
            for acc in accuracy:
                
                if len(snr_acc_all_sub[sub][contrast][acc]) != 0:
                    
                    snr_1 = snr_acc_all_sub[sub][contrast][acc][0][idx_1_2hz, idx_roi[roi]] 
                    snr_2 = snr_acc_all_sub[sub][contrast][acc][0][idx_2_4hz, idx_roi[roi]]
                    snr_3 = snr_acc_all_sub[sub][contrast][acc][0][idx_3_6hz, idx_roi[roi]]
                    snr_6 = snr_acc_all_sub[sub][contrast][acc][0][idx_6hz, idx_roi[roi]]
                    snr_12 = snr_acc_all_sub[sub][contrast][acc][0][idx_12hz, idx_roi[roi]]
                    snr_18 = snr_acc_all_sub[sub][contrast][acc][0][idx_18hz, idx_roi[roi]]
                     
                    new_data = [{'Subject': sub,
                                 'Contrast': contrast, 
                                 'Accuracy': acc,
                                 '1_2Hz': np.mean(snr_1, axis=0),
                                 '2_4Hz': np.mean(snr_2, axis=0),
                                 '3_6Hz': np.mean(snr_3, axis=0),
                                 '6Hz': np.mean(snr_6, axis=0),
                                 '12Hz': np.mean(snr_12, axis=0),
                                 '18Hz': np.mean(snr_18, axis=0)}]
                    
                    for row in new_data:
                        snr_roi_acc[f'roi{roi+1}'] = pd.concat([snr_roi_acc[f'roi{roi+1}'], pd.DataFrame(row, index=[0])], ignore_index=True)

# save data    
snr_roi_acc['roi1'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OCC_SNR_acc_all_sub.csv'), index=False)  
snr_roi_acc['roi2'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT1_SNR_acc_all_sub.csv'), index=False)  
snr_roi_acc['roi3'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT2_SNR_acc_all_sub.csv'), index=False)  

    
## Get snr per confidence rating at each frequency, subject, and ROI

snr_roi_conf = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}
for sub in snr_pas_all_sub.keys(): 
               
    for roi in range(len(idx_roi)):   
        contrasts = snr_pas_all_sub[sub].keys()
        
        for contrast in contrasts:    
               
            for conf in range(1, 5):
                
                if len(snr_conf_all_sub[sub][contrast][conf]) != 0:
                    
                    snr_1 = snr_conf_all_sub[sub][contrast][conf][0][idx_1_2hz, idx_roi[roi]] 
                    snr_2 = snr_conf_all_sub[sub][contrast][conf][0][idx_2_4hz, idx_roi[roi]]
                    snr_3 = snr_conf_all_sub[sub][contrast][conf][0][idx_3_6hz, idx_roi[roi]]
                    snr_6 = snr_conf_all_sub[sub][contrast][conf][0][idx_6hz, idx_roi[roi]]
                    snr_12 = snr_conf_all_sub[sub][contrast][conf][0][idx_12hz, idx_roi[roi]]
                    snr_18 = snr_conf_all_sub[sub][contrast][conf][0][idx_18hz, idx_roi[roi]]
                     
                    new_data = [{'Subject': sub,
                                 'Contrast': contrast, 
                                 'Confidence': conf,
                                 '1_2Hz': np.mean(snr_1, axis=0),
                                 '2_4Hz': np.mean(snr_2, axis=0),
                                 '3_6Hz': np.mean(snr_3, axis=0),
                                 '6Hz': np.mean(snr_6, axis=0),
                                 '12Hz': np.mean(snr_12, axis=0),
                                 '18Hz': np.mean(snr_18, axis=0)}]
                    
                    for row in new_data:
                        snr_roi_conf[f'roi{roi+1}'] = pd.concat([snr_roi_conf[f'roi{roi+1}'], pd.DataFrame(row, index=[0])], ignore_index=True)

# save data    
snr_roi_conf['roi1'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OCC_SNR_conf_all_sub.csv'), index=False)  
snr_roi_conf['roi2'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT1_SNR_conf_all_sub.csv'), index=False)  
snr_roi_conf['roi3'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_OT2_SNR_conf_all_sub.csv'), index=False)  
