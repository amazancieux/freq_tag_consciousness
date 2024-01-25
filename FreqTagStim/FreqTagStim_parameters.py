# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:17:29 2024

This script correspond to parameters that can eb changed for the frequency 
tagging task.

@author: Audrey Mazancieux

"""

# Other devices

IsEEG               = 0

# Stimuli caracteristics

CONTRASTS           = [0.5, 1, 1.5] # sublinal, threashold, and supraliminal contrasts
NSEQ_PER_CONTRATS   = 1

# Presentation caracteristics

REFRESH_RATE        = 60 # in Hz
STIM_FREQ           = 6 # in Hz
FADE                = 1.6667 # in sec
NB_STIM_SEQ         = 20
NB_TARGET_FIXCROSS  = 6 # mean number of detections in the attentional task
DUR_CHANGE_FIXCROSS = 3 # duration of change for fix cross (in number of stim)

# Constants for colours

BLACK               = (0, 0, 0)
GRAY_B              = (119, 136, 153)
BLUE                = (0, 0, 255)
DARK_GREY           = (96, 96, 96)
WHITE               = (255, 255, 255)

EXP_FONT            = "Arial"