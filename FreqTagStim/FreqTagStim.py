# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 10:52:01 2024

SESSION 2
-------------------------------------------------------------------------------
This task is a frequency tagging task which present sequences of picture stimuli 
degraded in several level of contrasts and with 2AFC and PAS questions at the 
end of each sequence. EEG option available. 
Presentation of stimuli corresponds to square wave (half on half off).

All modified parameters are in the FreqTagStim_parameters.py 

@author: Audrey Mazancieux

"""

import os
import glob
import random
import numpy as np
import csv
import pandas as pd
from expyriment import design, control, stimuli, misc
from expyriment.misc import data_preprocessing as dp
import serial

# =============================================================================

# Load task parameters
from FreqTagStim_parameters import *

# Define useful function

def get_subject(exp):
    """Function which gets the subject number and put every data results
    for this subject number in the same file."""
    subject = exp.subject
    try :
        os.mkdir(f'data/sub-{subject}')
    except FileExistsError :
        pass
    exp.data.rename(f'sub-{subject}/{exp.data.filename}')
    return subject


def create_vector_change_fixcross(vector, target_duration, nb_target):
    result = vector.copy()

    for i in range(nb_target):
        # Generate a random position
        position = random.randint(5, len(vector) - 5) # avoid the 5 first and last stim

        # Check if the position meets the minimum spacing requirement (target duration*2)
        for element in vector[position:position+target_duration*2]:
            if element != 0:
                break

        # Add 1 to the selected position according to target duration
        for t in range(target_duration):
            result[position+t] += 1

    return result


def save_csv_file(subject_id): 
    xpd_filename = glob.glob(os.path.join(f'./data/sub-{subject_id}',
                                    'Freq*.xpd'))[0]
    
    data = dp.read_datafile(xpd_filename)
    output_file = f'./data/sub-{subject_id}/FreqTagStim_sub-{subject_id}.csv'     
    dataframe = pd.DataFrame(data[0], columns=data[1])
    dataframe.to_csv(output_file, sep=',', index=False)
    
    return output_file


def scan():
    available = []
    for i in range(256) : 
        try : 
            s= serial.Serial("COM"+str(i))
            available.append((s.portstr))
            s.close()
        except serial.SerialException :
            pass
    return available


def send_trigger(trigger):
    try : 
         port.write(trigger)
    except:
        pass
# https://www.ascii-code.com/ to find the corresponding b'' decimal trigger its decimal to hexadecimal


# =============================================================================
# INITIALISATION
# =============================================================================  

# get info for EEG port 
print(scan())
com_input = input("The port COM is: ")
port = serial.Serial("COM"+com_input, baudrate = 115200)

# initialize EEG port
if IsEEG == 1:
    port = serial.Serial("COM4", baudrate=115200, timeout=1)  # old baudrate:   9600 115200

exp = design.Experiment(name="Contrast face experiment",
                        background_colour=DARK_GREY,
                        foreground_colour=misc.constants.C_BLACK)  

exp.add_experiment_info("ContrastFace")
control.defaults.window_mode = True # ask for switch to windows mode
control.initialize(exp)

# set refresh rate 
exp.screen.refresh_rate = REFRESH_RATE


# =============================================================================
# DEFINE USEFUL SCREENS
# =============================================================================

# Fixation cross screens
fixation_cross = stimuli.FixCross(size=(30, 30), line_width=2)
fixation_cross.preload()

fixation_cross_blue = stimuli.FixCross(size=(30, 30), line_width=2, colour=(BLUE))
fixation_cross_blue.preload()

# Ready screen
ready_screen = stimuli.TextLine("Ready? Appuyez sur S pour lancer la séquence", 
                                 position=(0,0),
                                 text_size=21,
                                 text_font=EXP_FONT)
ready_screen.preload()

# Decision screen
decision_screen = stimuli.TextLine("Avez-vous vu des visages de femmes ou d'hommes? Appuyez sur F ou H", 
                                   position=(0,0),
                                   text_size=21,
                                   text_font=EXP_FONT)
decision_screen.preload()

# PAS screen
pas_screen = stimuli.TextLine("PAS? Appuyez sur 0, 1, 2 ou 3", 
                              position=(0,0),
                              text_size=21,
                              text_font=EXP_FONT)
pas_screen.preload()


# =============================================================================
# GENERATE STIMULI LISTS AND SEQUENCES
# =============================================================================

# get path to images
root_path = os.getcwd()
path_to_images = os.path.join(root_path, "./Stimuli") 

# generate sequences with 1 face every 5th item 
sequence = [0] * NB_STIM_SEQ
for i in range(4, len(sequence), 5):
    sequence[i] = 1

# get male face stimuli
m_face_0_files = sorted(glob.glob(os.path.join(path_to_images, '0%', 'Face', 'male', '*.png')))
m_face_05_files = sorted(glob.glob(os.path.join(path_to_images, '0.5%', 'Face', 'male', '*.png')))
m_face_1_files = sorted(glob.glob(os.path.join(path_to_images, '1%', 'Face', 'male', '*.png')))
m_face_15_files = sorted(glob.glob(os.path.join(path_to_images, '1.5%', 'Face', 'male', '*.png')))
m_face_2_files = sorted(glob.glob(os.path.join(path_to_images, '2%', 'Face', 'male', '*.png')))
m_face_25_files = sorted(glob.glob(os.path.join(path_to_images, '2.5%', 'Face', 'male', '*.png')))
m_face_3_files = sorted(glob.glob(os.path.join(path_to_images, '3%', 'Face', 'male', '*.png')))
m_face_35_files = sorted(glob.glob(os.path.join(path_to_images, '3.5%', 'Face', 'male', '*.png')))

m_face_dict = {'Contrast_0': m_face_0_files,
               'Contrast_0.5': m_face_05_files,
             'Contrast_1': m_face_1_files,
             'Contrast_1.5': m_face_15_files,
             'Contrast_2': m_face_2_files,
             'Contrast_2.5': m_face_25_files,
             'Contrast_3': m_face_3_files,
             'Contrast_3.5': m_face_35_files}

# get female face stimuli
f_face_0_files = sorted(glob.glob(os.path.join(path_to_images, '0%', 'Face', 'female', '*.png')))
f_face_05_files = sorted(glob.glob(os.path.join(path_to_images, '0.5%', 'Face', 'female', '*.png')))
f_face_1_files = sorted(glob.glob(os.path.join(path_to_images, '1%', 'Face', 'female', '*.png')))
f_face_15_files = sorted(glob.glob(os.path.join(path_to_images, '1.5%', 'Face', 'female', '*.png')))
f_face_2_files = sorted(glob.glob(os.path.join(path_to_images, '2%', 'Face', 'female', '*.png')))
f_face_25_files = sorted(glob.glob(os.path.join(path_to_images, '2.5%', 'Face', 'female', '*.png')))
f_face_3_files = sorted(glob.glob(os.path.join(path_to_images, '3%', 'Face', 'female', '*.png')))
f_face_35_files = sorted(glob.glob(os.path.join(path_to_images, '3.5%', 'Face', 'female', '*.png')))

f_face_dict = {'Contrast_0': f_face_0_files,
               'Contrast_0.5': f_face_05_files,
             'Contrast_1': f_face_1_files,
             'Contrast_1.5': f_face_15_files,
             'Contrast_2': f_face_2_files,
             'Contrast_2.5': f_face_25_files,
             'Contrast_3': f_face_3_files,
             'Contrast_3.5': f_face_35_files}

# get item stimuli
item_0_files = sorted(glob.glob(os.path.join(path_to_images, '0%', 'NonFace', '*.png')))
item_05_files = sorted(glob.glob(os.path.join(path_to_images, '0.5%', 'NonFace', '*.png')))
item_1_files = sorted(glob.glob(os.path.join(path_to_images, '1%', 'NonFace', '*.png')))
item_15_files = sorted(glob.glob(os.path.join(path_to_images, '1.5%', 'NonFace', '*.png')))
item_2_files = sorted(glob.glob(os.path.join(path_to_images, '2%', 'NonFace', '*.png')))
item_25_files = sorted(glob.glob(os.path.join(path_to_images, '2.5%', 'NonFace', '*.png')))
item_3_files = sorted(glob.glob(os.path.join(path_to_images, '3%', 'NonFace', '*.png')))
item_35_files = sorted(glob.glob(os.path.join(path_to_images, '3.5%', 'NonFace', '*.png')))

item_dict = {'Contrast_0': item_0_files,
             'Contrast_0.5': item_05_files,
             'Contrast_1': item_1_files,
             'Contrast_1.5': item_15_files,
             'Contrast_2': item_2_files,
             'Contrast_2.5': item_25_files,
             'Contrast_3': item_3_files,
             'Contrast_3.5': item_35_files}

# create N sequences per selected contrast 
generated_seq = {'contrast_type': [],
                 'sequence': []}

for contrast in CONTRASTS:       
    for n_seq in range(0, NSEQ_PER_CONTRATS):        
        seq = []
        for stim in sequence: 
                 
            if stim == 0:    
                # select a random index
                random_index = random.randrange(len(item_dict[f'Contrast_{contrast}']))
                # get item
                item = item_dict[f'Contrast_{contrast}'][random_index]
            
            else: 
                # select a random index
                random_index = random.randrange(len(item_dict[f'Contrast_{contrast}']))
                # get face
                if (n_seq % 2) == 0: 
                    face_type = 'male'
                    item = m_face_dict[f'Contrast_{contrast}'][random_index]
                else: 
                    face_type = 'female'
                    item = f_face_dict[f'Contrast_{contrast}'][random_index]
            
            # add it to seq
            seq.append(item)
                
        # add seq to block
        generated_seq['contrast_type'].append(f'Contrast_{contrast}_{face_type}')
        generated_seq['sequence'].append(seq)
            
        
# generate random order of sequence for this session  
shuffle_seq = list(range(0, len(generated_seq['sequence'])))  
random.shuffle(shuffle_seq)
       
# =============================================================================
# OTHER PARAMETERS
# =============================================================================

# Basic parameters from stimuli presentation
onset_stim = 1000/STIM_FREQ # one stimulus each onset (in ms)
t_one_frame = 1000/60 # time of one frame (in ms)
nb_frame_stim = round(onset_stim/t_one_frame) # nb of frame for stimulus + blanck presentation
nb_stim_fade = round(STIM_FREQ*FADE) # nb of stim for fade

# Computaing alpha values
# for experimental stimuli
alpha_exp = np.concatenate([np.zeros(round(nb_frame_stim/2)), np.repeat(100, round(nb_frame_stim/2))])

# for fade stimuli
alpha_values = np.linspace(5, 95, num=nb_frame_stim) # avoid extrem values (max=100)
alpha_fade_in = []
alpha_fade_out = []

for stim_fade in range(0, nb_frame_stim):
    alpha_stim = np.concatenate([np.zeros(round(nb_frame_stim/2)), 
                                 np.repeat(100-(alpha_values[stim_fade]), round(nb_frame_stim/2))])
    alpha_fade_in.append(alpha_stim)
    alpha_stim = np.concatenate([np.zeros(round(nb_frame_stim/2)), 
                                 np.repeat(alpha_values[stim_fade], round(nb_frame_stim/2))])
    alpha_fade_out.append(alpha_stim)

alpha_fade_in = np.vstack(alpha_fade_in) 
alpha_fade_out = np.vstack(alpha_fade_out) 


# =============================================================================
# START THE EXPERIMENT
# =============================================================================

# Start Experiment
control.start(skip_ready_screen=True)
  
# Get subject number    
subject = get_subject(exp)


# =============================================================================
# TASK DESIGN
# =============================================================================

# Define moments for fix cross change  
vector = [0] * NB_STIM_SEQ
fixcross_vector = create_vector_change_fixcross(vector=vector, 
                                                target_duration=DUR_CHANGE_FIXCROSS,
                                                nb_target=NB_TARGET_FIXCROSS)

# Loop to create and add blocks
for block_num in range(0, len(generated_seq['contrast_type'])):

    contrast_type = generated_seq['contrast_type'][block_num]
    block = design.Block(name=f'{contrast_type}') 
    
    # get total number of stim including fade period
    tot_stim_with_fade = NB_STIM_SEQ + nb_stim_fade*2
    
    # add trials to each block using generated seq
    for trial_num in range(len(tot_stim_with_fade)):
        
        if trial_num >= nb_stim_fade:             
            # get fade stim
            stimuli.Picture(item_dict['sequence'][block_num][]
            
        
        
            # get expe stim 
            stim = stimuli.Picture(generated_seq['sequence'][block_num][trial_num-nb_stim_fade])
            
            # add fix cross (black or blue depending on condition) 
            if fixcross_vector[trial_num] == 0:                  
                fixation_cross.plot(stim) 
                stim.preload()
                stim_name = generated_seq['sequence'][block_num][trial_num-nb_stim_fade]    
                trial = design.Trial()
                trial.set_factor("stim_name", stim_name)
                trial.add_stimulus(stim)         
            else:            
                fixation_cross_blue.plot(stim) 
                stim.preload()
                stim_name = generated_seq['sequence'][block_num][trial_num-nb_stim_fade]    
                trial = design.Trial()
                trial.set_factor("stim_name", stim_name)
                trial.add_stimulus(stim)  
                
        block.add_trial(trial)
    
    # add the block to the experiment
    exp.add_block(block)



# =============================================================================
# DEFINE THE VARIABLES THAT WILL BE SAVED 
# =============================================================================

fixcross_task = {'block': [],
                 'trial': [],
                 'frame': [],
                 'fixcross_resp': [],
                 'fixcross_time': [],
                 'fixcross_vector': []}
        

trial_variable_names = ["subject",
                        "contrast_type",
                        "trial",
                        "frame",
                        "stimulus",
                        "t_on_frame",
                        "t_off_frame"]
    
block_variable_names = ["t_decision_on",
                        "decision_resp",
                        "rt_decision",
                        "t_decision_off",
                        "t_pas_on",
                        "pas_resp",
                        "pas_rt",
                        "t_pas_off"]

all_variable_names = trial_variable_names + block_variable_names

# =============================================================================
# START THE TASK
# =============================================================================

# Initialize frame count
frame = 0

# Display the sequence
exp.keyboard.check()    
for block in exp.blocks: 
    
    # display ready screen
    ready_screen.present()
    exp.keyboard.wait_char("s")
    
    # fixation for a random duration between 2 and 5 sec
    fixation_cross.present() 
    exp.clock.wait(random.randint(2, 5))
            
    # fade-in
    
    # experimental trials        
    for t, trial in enumerate(block.trials):
        
        for alpha in alpha_exp:
            
            # update frame
            frame = frame+1 
            
            # space press for attentionnal task
            space_resp = exp.keyboard.check(misc.constants.K_SPACE)  
            if space_resp:
                space_time = exp.clock.time
            else: 
                space_resp = None
                space_time = None
            fixcross_task['block'].append(block.name)
            fixcross_task['trial'].append(trial.id+1)
            fixcross_task['frame'].append(frame)
            fixcross_task['fixcross_resp'].append(space_resp)
            fixcross_task['fixcross_time'].append(space_time)
            fixcross_task['fixcross_vector'].append(fixcross_vector[t])
            exp.keyboard.clear()
            
            if alpha == 0:    
                # present stim with fix cross
                t_on_frame = exp.clock.time
                trial.stimuli[0].present()
                t_off_frame = exp.clock.time
            
            else:
                # present only fix cross
                if fixcross_vector[t] == 0:  
                    fixation_cross.present()
                else: 
                    fixation_cross_blue.present()                
          
            # save trial data
            trial_data = [subject,
                          block.name,
                          trial.id+1,
                          trial.get_factor("stim_name"),
                          frame,
                          t_on_frame,
                          t_off_frame] + \
                [np.nan for _ in range(len(block_variable_names))]
        
            exp.data.add(trial_data)
    
    # fade out
              
    # decision    
    t_decision_on = exp.clock.time
    decision_screen.present()
    decision_resp, rt_decision = exp.keyboard.wait([misc.constants.K_h, misc.constants.K_f])
    t_decision_off = exp.clock.time
    
    # PAS response
    t_pas_on = exp.clock.time  
    pas_screen.present()
    pas_resp, pas_rt = exp.keyboard.wait([misc.constants.K_0, misc.constants.K_1,
                                          misc.constants.K_2, misc.constants.K_3])
    t_pas_off = exp.clock.time
      
    # save response data
    data_saved_table = [subject,
                        block.name,
                        "None",
                        "None",
                        "None",
                        "None",
                        "None",
                        t_decision_on,
                        decision_resp,
                        rt_decision,
                        t_decision_off,
                        t_pas_on,
                        pas_resp,
                        pas_rt,
                        t_pas_off]
    
    exp.data.add(data_saved_table)
        
exp.data.add_variable_names(all_variable_names)

# =============================================================================
# SAVE DATA AND END EXPERIMENT
# =============================================================================

control.end(goodbye_text="Fin de l'expérience. Merci pour votre participation !",
                goodbye_delay=2000)

# fixcross in dataframe and save
fixcross_dataframe = pd.DataFrame(fixcross_task)
file = f'./data/sub-{subject}/fixCross_sub-{subject}.csv' 
fixcross_dataframe.to_csv(file, sep=',', index=False)

# convert sequences data 
csv_filename = save_csv_file(subject_id=subject)




