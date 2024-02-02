# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:17:29 2024

This script correspond to parameters that can eb changed for the frequency 
tagging task.

@author: Audrey Mazancieux

"""

# Main infos

SESSION             = 1 # 1 or 2 
IsEEG               = 0

# Stimuli caracteristics

CONTRASTS           = [3] # sublinal, threashold, and supraliminal contrasts
NSEQ_PER_CONTRATS   = 1

# Presentation caracteristics

REFRESH_RATE        = 60 # in Hz
STIM_FREQ           = 6 # in Hz
FADE                = 1.6667 # in sec
NB_STIM_SEQ         = 240
NB_TARGET_FIXCROSS  = 6 # mean number of detections in the attentional task
DUR_CHANGE_FIXCROSS = 4 # duration of change for fix cross (in number of stim)

# Constants for colours

BLACK               = (0, 0, 0)
RED                 = (255, 0, 0)
BLUE                = (0, 0, 255)
DARK_GREY           = (96, 96, 96)
WHITE               = (255, 255, 255)

EXP_FONT            = "Arial"