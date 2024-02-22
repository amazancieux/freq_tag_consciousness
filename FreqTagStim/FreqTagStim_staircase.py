# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:09:19 2024

SESSION 1
-------------------------------------------------------------------------------
This task is staircase version of a frequency tagging task which present 
sequences of picture stimuli degraded in several level of contrasts and with 
2AFC and PAS questions at the end of each sequence. 
Presentation of stimuli corresponds to square wave (half on half off).

All modified parameters are in the FreqTagStim_parameters.py 

@author: Audrey Mazancieux

"""

import os
import glob
import random
import numpy as np
import pandas as pd
from expyriment import design, control, stimuli, misc
from expyriment.misc import data_preprocessing as dp
import FreqTagStim_parameters as param

# =============================================================================

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
        position = random.randint(15, len(vector) - 15) # avoid the 15 first and last stim

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
                                    '*staircase*.xpd'))[0]
    
    data = dp.read_datafile(xpd_filename)
    output_file = f'./data/sub-{subject_id}/FreqTagStim_staircase_sub-{subject_id}.csv'     
    dataframe = pd.DataFrame(data[0], columns=data[1])
    dataframe.to_csv(output_file, sep=',', index=False)
    
    return output_file


# =============================================================================
# INITIALISATION
# =============================================================================  

exp = design.Experiment(name="Contrast face experiment staircase",
                        background_colour=param.DARK_GREY,
                        foreground_colour=misc.constants.param.C_BLACK)  

exp.add_experiment_info("ContrastFace_staircase")
control.defaults.window_mode = True # ask for switch to windows mode
control.initialize(exp)

# set refresh rate 
exp.screen.refresh_rate = param.REFRESH_RATE


# =============================================================================
# DEFINE USEFUL SCREENS
# =============================================================================

# Fixation cross screens
fixation_cross_white = stimuli.FixCross(size=(20, 20), line_width=2, colour=(param.WHITE))
fixation_cross_white.preload()

fixation_cross_blue = stimuli.FixCross(size=(20, 20), line_width=2, colour=(param.BLUE))
fixation_cross_blue.preload()

# Ready screen
ready_screen = stimuli.TextLine("Prêt.e? Appuyez sur ESPACE pour lancer la séquence", 
                               position=(0,0),
                               text_size=21,
                               text_font=param.EXP_FONT)
ready_screen.preload()

# Decision screen
decision_screen = stimuli.TextBox("Avez-vous vu des visages de femmes ou d'hommes? Appuyez sur F ou H", 
                                  size=(600, 600),
                                  position=(0,0),
                                  text_size=21,
                                  text_font=param.EXP_FONT)
decision_screen.preload()

# PAS screen
pas_screen = stimuli.TextBox("Impression des visages? Appuyez sur 1 (aucune impression des visages), 2 (bref un aperçu des visages), 3 (une expérience presque claire des visages) ou 4 (une expérience claire des visages)", 
                              size=(600, 600),
                              position=(0,0),
                              text_size=21,
                              text_font=param.EXP_FONT)
pas_screen.preload()

# Confident screen
conf_screen = stimuli.TextBox("Quel est votre confidence dans votre réponse? Appuyez sur 1 pour 50% (réponse au hasard), 2 pour 60%, 3 pour 70%, 4 pour 80%, 5 pour 90% ou 6 pour 100% (très confiant.e)", 
                              size=(600, 600),
                              position=(0,0),
                              text_size=21,
                              text_font=param.EXP_FONT)
conf_screen.preload()


# =============================================================================
# GENERATE STIMULI LISTS AND SEQUENCES
# =============================================================================

# Basic parameters from stimuli presentation
onset_stim = 1000/param.STIM_FREQ # one stimulus each onset (in ms)
t_one_frame = 1000/60 # time of one frame (in ms)
nb_frame_stim = round(onset_stim/t_one_frame) # nb of frame for stimulus + blanck presentation
nb_stim_fade = round(param.param.STIM_FREQ*param.FADE) # nb of stim for fade

# get path to images
root_path = os.getcwd()
path_to_images = os.path.join(root_path, "./Stimuli") 

# generate sequences with 1 face every 5th item 
sequence = [0] * (param.NB_STIM_SEQ+nb_stim_fade*2)
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
m_face_4_files = sorted(glob.glob(os.path.join(path_to_images, '4%', 'Face', 'male', '*.png')))

m_face_dict = {'Contrast_0.0': m_face_0_files,
               'Contrast_0.5': m_face_05_files,
             'Contrast_1.0': m_face_1_files,
             'Contrast_1.5': m_face_15_files,
             'Contrast_2.0': m_face_2_files,
             'Contrast_2.5': m_face_25_files,
             'Contrast_3.0': m_face_3_files,
             'Contrast_3.5': m_face_35_files,
              'Contrast_4.0': m_face_4_files}

# get female face stimuli
f_face_0_files = sorted(glob.glob(os.path.join(path_to_images, '0%', 'Face', 'female', '*.png')))
f_face_05_files = sorted(glob.glob(os.path.join(path_to_images, '0.5%', 'Face', 'female', '*.png')))
f_face_1_files = sorted(glob.glob(os.path.join(path_to_images, '1%', 'Face', 'female', '*.png')))
f_face_15_files = sorted(glob.glob(os.path.join(path_to_images, '1.5%', 'Face', 'female', '*.png')))
f_face_2_files = sorted(glob.glob(os.path.join(path_to_images, '2%', 'Face', 'female', '*.png')))
f_face_25_files = sorted(glob.glob(os.path.join(path_to_images, '2.5%', 'Face', 'female', '*.png')))
f_face_3_files = sorted(glob.glob(os.path.join(path_to_images, '3%', 'Face', 'female', '*.png')))
f_face_35_files = sorted(glob.glob(os.path.join(path_to_images, '3.5%', 'Face', 'female', '*.png')))
f_face_4_files = sorted(glob.glob(os.path.join(path_to_images, '4%', 'Face', 'female', '*.png')))

f_face_dict = {'Contrast_0.0': f_face_0_files,
               'Contrast_0.5': f_face_05_files,
             'Contrast_1.0': f_face_1_files,
             'Contrast_1.5': f_face_15_files,
             'Contrast_2.0': f_face_2_files,
             'Contrast_2.5': f_face_25_files,
             'Contrast_3.0': f_face_3_files,
             'Contrast_3.5': f_face_35_files,
             'Contrast_4.0': f_face_4_files}

# get item stimuli
item_0_files = sorted(glob.glob(os.path.join(path_to_images, '0%', 'NonFace', '*.png')))
item_05_files = sorted(glob.glob(os.path.join(path_to_images, '0.5%', 'NonFace', '*.png')))
item_1_files = sorted(glob.glob(os.path.join(path_to_images, '1%', 'NonFace', '*.png')))
item_15_files = sorted(glob.glob(os.path.join(path_to_images, '1.5%', 'NonFace', '*.png')))
item_2_files = sorted(glob.glob(os.path.join(path_to_images, '2%', 'NonFace', '*.png')))
item_25_files = sorted(glob.glob(os.path.join(path_to_images, '2.5%', 'NonFace', '*.png')))
item_3_files = sorted(glob.glob(os.path.join(path_to_images, '3%', 'NonFace', '*.png')))
item_35_files = sorted(glob.glob(os.path.join(path_to_images, '3.5%', 'NonFace', '*.png')))
item_4_files = sorted(glob.glob(os.path.join(path_to_images, '4%', 'NonFace', '*.png')))

item_dict = {'Contrast_0.0': item_0_files,
             'Contrast_0.5': item_05_files,
             'Contrast_1.0': item_1_files,
             'Contrast_1.5': item_15_files,
             'Contrast_2.0': item_2_files,
             'Contrast_2.5': item_25_files,
             'Contrast_3.0': item_3_files,
             'Contrast_3.5': item_35_files,
             'Contrast_4.0': item_4_files}

# create N sequences per selected contrast 
generated_seq = {'contrast_type': [],
                 'sequence': []}

for contrast in item_dict.keys():       
    for n_seq in range(0, param.NSEQ_CONTRAST_SESS1):        
        seq = []
        for stim in sequence: 
                 
            if stim == 0:    
                # select a random index
                random_index = random.randrange(len(item_dict[contrast]))
                # get item
                item = item_dict[contrast][random_index]
            
            else:                
                if (n_seq % 2) == 0: 
                    # select a random index
                    random_index = random.randrange(len(m_face_dict[contrast]))
                    # get face
                    face_type = 'male'
                    item = m_face_dict[contrast][random_index]
                else: 
                    # select a random index
                    random_index = random.randrange(len(f_face_dict[contrast]))
                    # get face
                    face_type = 'female'
                    item = f_face_dict[contrast][random_index]
            
            # add it to seq
            seq.append(item)
                
        # add seq to block
        generated_seq['contrast_type'].append(f'{contrast}_{face_type}')
        generated_seq['sequence'].append(seq)
            
       
# =============================================================================
# OTHER PARAMETERS
# =============================================================================

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
vector = [0] * (param.NB_STIM_SEQ+nb_stim_fade*2)
fixcross_vector = create_vector_change_fixcross(vector=vector, 
                                                target_duration=param.DUR_CHANGE_FIXCROSS,
                                                nb_target=param.NB_TARGET_FIXCROSS)

# Loop to create and add blocks
for block_num in range(len(generated_seq['contrast_type'])):

    contrast_type = generated_seq['contrast_type'][block_num]
    block = design.Block(name=f'{contrast_type}') 
    
    
    # add trials to each block using generated seq
    for trial_num in range(len(generated_seq['sequence'][block_num])):
        
        trial = design.Trial()
        
        if trial_num < nb_stim_fade:             
            # get fade in stim without noise
            stim1 = stimuli.Picture(generated_seq['sequence'][block_num][trial_num])
            fixation_cross_white.plot(stim1)
            stim1.preload()
            # get fade in stim with noise
            stim2 = stimuli.Picture(generated_seq['sequence'][block_num][trial_num])                   
            stim2.add_noise(grain_size=1, percentage=alpha_fade_in[trial_num, -1], colour=param.DARK_GREY)      
            fixation_cross_white.plot(stim2) 
            stim2.preload()  
            
            trial.add_stimulus(stim1)
            trial.add_stimulus(stim2)
            stim_name = 'fade-in'
            
        elif trial_num > nb_stim_fade+param.NB_STIM_SEQ:
            # get fade out stim without noise
            stim1 = stimuli.Picture(generated_seq['sequence'][block_num][trial_num])
            fixation_cross_white.plot(stim1)
            stim1.preload()
            # get fade out stim with noise
            stim2 = stimuli.Picture(generated_seq['sequence'][block_num][trial_num])
            stim2.add_noise(grain_size=1, percentage=alpha_fade_out[trial_num-(nb_stim_fade+param.NB_STIM_SEQ), -1], colour=param.DARK_GREY)         
            fixation_cross_white.plot(stim2)    
            stim2.preload()  
            
            trial.add_stimulus(stim1)
            trial.add_stimulus(stim2)
            stim_name = 'fade-out'
        
        else:
            # get expe stim 
            stim = stimuli.Picture(generated_seq['sequence'][block_num][trial_num])               
            # add fix cross (red or blue depending on condition) 
            if fixcross_vector[trial_num] == 0:                  
                fixation_cross_white.plot(stim)         
            else:            
                fixation_cross_blue.plot(stim) 
            stim.preload()
            trial.add_stimulus(stim)            
            stim_name = generated_seq['sequence'][block_num][trial_num] 
                  
        trial.set_factor("stim_name", stim_name)
                
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
    
block_variable_names = ["t_pas_on",
                        "pas_resp",
                        "pas_rt",
                        "t_pas_off",
                        "t_decision_on",
                        "decision_resp",
                        "correct_response",
                        "rt_decision",
                        "t_decision_off",
                        "t_conf_on",
                        "conf_resp",
                        "conf_rt",
                        "t_conf_off"]

all_variable_names = trial_variable_names + block_variable_names

# =============================================================================
# START THE TASK
# =============================================================================

# Initialize frame count
frame = 0

# Display the sequence
exp.keyboard.check()    

for i_block in range(0, param.NSEQ_TOT_SESS1): 
    
    # staircase initial contrast
    if i_block == 0:
        c = 3
    
    # get relevant blocks and shuffle   
    index_contrast_m = [index for index, value in enumerate(generated_seq['contrast_type']) if value == f'Contrast_{c}_male']
    index_contrast_f = [index for index, value in enumerate(generated_seq['contrast_type']) if value == f'Contrast_{c}_female']
    index_all = index_contrast_m + index_contrast_f
    random.shuffle(index_all)

    # select one block 
    block = exp.blocks[index_all[0]]

    # display ready screen
    ready_screen.present()
    exp.keyboard.wait(misc.constants.K_SPACE)
    
    # fixation for a random duration between 2 and 5 sec
    fixation_cross_white.present() 
    exp.clock.wait(random.randint(2, 5))
    
    # loop for sequences        
    for t, trial in enumerate(block.trials):
        
        # Fade-in 
        if t < nb_stim_fade:
            for alpha in alpha_fade_in[t]:
     
                # update frame
                frame = frame+1 
                
                # present stim
                t_on_frame = exp.clock.time
                if alpha == 0:          
                    trial.stimuli[0].present()
                else:
                    trial.stimuli[1].present()                    
                t_off_frame = exp.clock
            
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
        
        # Fade-out
        if t > nb_stim_fade+param.NB_STIM_SEQ:
            for alpha in alpha_fade_out[t-(nb_stim_fade+param.NB_STIM_SEQ)]:
                     
                # update frame
                frame = frame+1 
                
                # present stim
                t_on_frame = exp.clock.time
                if alpha == 0:          
                    trial.stimuli[0].present()
                else:
                    trial.stimuli[1].present()                    
                t_off_frame = exp.clock.time
                
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
 
        # Experimental trials
        else: 
            for alpha in alpha_exp:
                
                # update frame
                t_on_frame = exp.clock.time
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
                    trial.stimuli[0].present()   
                    frame_name = trial.get_factor("stim_name")
                else:
                    # present only fix cross
                    frame_name = 'fix_cross'
                    if fixcross_vector[t] == 0:  
                        fixation_cross_white.present()
                    else: 
                        fixation_cross_blue.present() 
                        
                t_off_frame = exp.clock.time
          
                # save trial data
                trial_data = [subject,
                              block.name,
                              trial.id+1,
                              frame_name,
                              frame,
                              t_on_frame,
                              t_off_frame] + \
                    [np.nan for _ in range(len(block_variable_names))]
            
                exp.data.add(trial_data)

    # PAS response
    t_pas_on = exp.clock.time  
    pas_screen.present()
    pas_resp, pas_rt = exp.keyboard.wait([misc.constants.K_1, misc.constants.K_2,
                                          misc.constants.K_3, misc.constants.K_4])
    t_pas_off = exp.clock.time
        
    # decision    
    t_decision_on = exp.clock.time
    decision_screen.present()
    decision_resp, rt_decision = exp.keyboard.wait([misc.constants.K_h, misc.constants.K_f])
    t_decision_off = exp.clock.time
    
    # conf response
    t_conf_on = exp.clock.time  
    conf_screen.present()
    conf_resp, conf_rt = exp.keyboard.wait([misc.constants.K_1, misc.constants.K_2,
                                          misc.constants.K_3, misc.constants.K_4,
                                          misc.constants.K_5, misc.constants.K_6])
    t_conf_off = exp.clock.time
    
    # get correct resp
    if "f" in block.name:
        correct_response = 102 # ASCII for f key
    else: 
        correct_response = 104 # ASCII for m key
    
    # staircase
    accuracy = (correct_response == decision_resp)
    if accuracy == True: 
        if c == 0:
            c = 0
        else:
            c = c-0.5
    else: 
            c = c+0.5
      
    # save response data
    data_saved_table = [subject,
                        block.name,
                        "None",
                        "None",
                        "None",
                        "None",
                        "None",
                        t_pas_on,
                        pas_resp,
                        pas_rt,
                        t_pas_off,
                        t_decision_on,
                        decision_resp,
                        correct_response,
                        rt_decision,
                        t_decision_off,
                        t_conf_on,
                        conf_resp,
                        conf_rt,
                        t_conf_off]
    
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
