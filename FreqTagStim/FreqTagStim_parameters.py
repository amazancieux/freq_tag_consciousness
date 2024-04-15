# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:17:29 2024

This script correspond to parameters that can eb changed for the frequency 
tagging task.

@author: Audrey Mazancieux

"""

# Stimuli caracteristics

CONTRASTS            = [1] # one or more contrasts for session 2
NSEQ_PER_CONTRATS    = 50 
NSEQ_NOREPORT        = 20
NSEQ_CONTRAST_SESS1  = 6 # even number for equal male and female seq
NSEQ_TOT_SESS1       = 50

# Presentation caracteristics

REFRESH_RATE         = 60 # in Hz
STIM_FREQ            = 6 # in Hz
FADE                 = 1.6667 # in sec
NB_STIM_SEQ          = 240
NB_TARGET_FIXCROSS   = 6 # mean number of detections in the attentional task
DUR_CHANGE_FIXCROSS  = 4 # duration of change for fix cross (in number of stim)

# Constants for colours

BLACK                = (0, 0, 0)
RED                  = (255, 0, 0)
BLUE                 = (0, 0, 255)
DARK_GREY            = (96, 96, 96)
WHITE                = (255, 255, 255)

EXP_FONT             = "Arial"

# EEG triggers in hexadecimal and its corresponding ASCII decimal

TRIGGER_SEQ = b'\x0A' # DEC 10

TRIGGER_STIM = b'\x14' # DEC 20
TRIGGER_FACE = b'\x15' # DEC 21

TRIGGER_PAS1 = b'\x1F' # DEC 31
TRIGGER_PAS2 = b'\x20' # DEC 32
TRIGGER_PAS3 = b'\x21' # DEC 33
TRIGGER_PAS4 = b'\x22' # DEC 34

TRIGGER_CORRECT = b'\x28' # DEC 40
TRIGGER_INCORRECT = b'\x29' # DEC 41

TRIGGER_CONF1 = b'\x33' # DEC 51
TRIGGER_CONF2 = b'\x34' # DEC 52
TRIGGER_CONF3 = b'\x35' # DEC 53
TRIGGER_CONF4 = b'\x36' # DEC 54
