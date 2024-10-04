# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 09:21:21 2024

@author: Audrey Mazancieux

Frequency analyses: extract SNR at given frequencies and plots
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import pickle
from mne.viz import plot_topomap
from matplotlib import pyplot as plt

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

FACE_FREQ = 1.2
IMAGE_FREQ = 6

# =============================================================================

## Load SNR data from preproc

with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'bads_all_subjects.pickle'), 'rb') as f:
    bads_all_sub = pickle.load(f) 
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_contrast_all_subjects.pickle'), 'rb') as f:
    snr_contrast_all_sub = pickle.load(f)

with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_pas_all_subjects.pickle'), 'rb') as f:
    snr_pas_all_sub = pickle.load(f) 
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_acc_all_subjects.pickle'), 'rb') as f:
    snr_acc_all_sub = pickle.load(f) 
    
with open(os.path.join(ROOT_DIR, EEG_DIR, DATA_DIR, 'SNR_conf_all_subjects.pickle'), 'rb') as f:
    snr_conf_all_sub = pickle.load(f) 
   
    
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
im, cm = plot_topomap(snr_2_4_C1, epochs.info, axes=ax)
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
im, cm = plot_topomap(snr_2_4_C2, epochs.info, axes=ax)
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
idx_ROIs = [bads_all_sub['channels'].index(electrod) for electrod in SELECTED_ELECTRODS]  

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

# from preprint
# ROI1 = ['Iz', 'Oz'] 
# ROI2 = ['O2', 'PO4', 'PO8']
# ROI3 = ['O1', 'PO3', 'PO7'] 
# ROI4 = ['P6', 'P8', 'P10'] 
# ROI5 = ['P5', 'P7', 'P9']

# from Quek & de Heering (2024)
ROI_OCC = ['Iz', 'Oz', 'O1', 'O2'] # for 6 Hz
ROI_OT_1 = ['O1', 'PO3', 'PO7', 'P7', 'P9'] # face left
ROI_OT_2 = ['O2', 'PO4', 'PO8', 'P8', 'P10'] # face right
                                                               
# idx_ROI1 = [bads_all_sub['channels'].index(electrod) for electrod in ROI1]
# idx_ROI2 = [bads_all_sub['channels'].index(electrod) for electrod in ROI2]
# idx_ROI3 = [bads_all_sub['channels'].index(electrod) for electrod in ROI3]
# idx_ROI4 = [bads_all_sub['channels'].index(electrod) for electrod in ROI4]
# idx_ROI5 = [bads_all_sub['channels'].index(electrod) for electrod in ROI5]  

idx_ROI_OCC = [bads_all_sub['channels'].index(electrod) for electrod in ROI_OCC]
idx_ROI_OT_1 = [bads_all_sub['channels'].index(electrod) for electrod in ROI_OT_1]
idx_ROI_OT_2 = [bads_all_sub['channels'].index(electrod) for electrod in ROI_OT_2]

# idx_roi = [idx_ROI1, idx_ROI2, idx_ROI3, idx_ROI4, idx_ROI5]
idx_roi = [idx_ROI_OCC, idx_ROI_OT_1, idx_ROI_OT_2]
snr_roi_pas = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}


## Get snr per PAS at each frequency, subject, and ROI

for subject in snr_pas_all_sub.keys(): 
               
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
                     
                    new_data = [{'Subject': subject,
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
# snr_roi_pas['roi4'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_4_SNR_PAS_all_sub.csv'), index=False)  
# snr_roi_pas['roi5'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_5_SNR_PAS_all_sub.csv'), index=False)  
                      

## Get snr per accuracy at each frequency, subject, and ROI

snr_roi_acc = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}
accuracy = snr_acc_all_sub['sub-3']['1%'].keys()
for subject in snr_pas_all_sub.keys(): 
               
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
                     
                    new_data = [{'Subject': subject,
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
# snr_roi_acc['roi4'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_4_SNR_acc_all_sub.csv'), index=False)  
# snr_roi_acc['roi5'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_5_SNR_acc_all_sub.csv'), index=False) 

    
## Get snr per confidence rating at each frequency, subject, and ROI

snr_roi_conf = {f'roi{roi+1}': pd.DataFrame(columns=['Subject', 'Contrast', 'PAS', '1_2Hz', '2_4Hz', '3_6Hz', '6Hz', '12Hz', '18Hz']) for roi in range(len(idx_roi))}
for subject in snr_pas_all_sub.keys(): 
               
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
                     
                    new_data = [{'Subject': subject,
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
# snr_roi_conf['roi4'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_4_SNR_conf_all_sub.csv'), index=False)  
# snr_roi_conf['roi5'].to_csv(os.path.join(ROOT_DIR, EEG_DIR, 'roi_5_SNR_conf_all_sub.csv'), index=False)  

