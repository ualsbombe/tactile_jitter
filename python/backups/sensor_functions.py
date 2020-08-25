#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:12:37 2020

@author: lau
"""

#%% IMPORTS

import mne
import numpy as np
from subprocess import check_output
from os.path import join, isfile
from sys import argv

#%% FIND RELEVANT FUNCTION AND USER



user = check_output('uname -n', shell=True)
if user == b'lau\n':
    project_path = '/home/lau/projects/cerebellar_clock'
    function_to_run = None
elif user[:6] == b'hyades':
    project_path = '/projects/MINDLAB2019_MEG-CerebellarClock'
    function_to_run = argv[1]
    
print(user)
print(function_to_run)
print(argv[:])

#%% PATHS

raw_path = join(project_path, 'raw')
scratch_path = join(project_path, 'scratch', 'tactile_jitter', 'MEG')
subjects_dir = join(project_path, 'scratch', 'tactile_jitter', 'freesurfer')

#%% FUNCTION DEFINITIONS
#%% GENERAL HELPER FUNCTIONS

def get_raw_path(subject, date, split_recording):
    full_date = date + '_000000'
    subject_path = join(raw_path, subject, full_date, 'MEG',
                    '001.tactile_jitter_raw', 'files')
                    
    return subject_path
                    
def get_save_path(subject, date):
    full_date = date + '_000000'
    save_path = join(scratch_path, subject, full_date)
    
    return save_path

def should_we_run(save_path, overwrite):
    run_function = overwrite or not isfile(save_path)
    if not run_function:
        print('not overwriting: ' + save_path)
    
    return run_function

def get_ico_string(ico):
    if ico == 3:
        ico_string = '1280'
    elif ico == 4:
        ico_string = '5120'
    elif ico == 5:
        ico_string = '20484'
        
    return ico_string

def hp_lp_string(h_freq, l_freq):
    if h_freq > 0 and l_freq > 0:
       prepend = 'hp_' + str(h_freq) + '_Hz_lp_' + str(l_freq) + '_Hz_'
    elif h_freq > 0 and l_freq == 0:
        prepend = 'lp_' + str(h_freq) + '_Hz_'
    elif l_freq > 0 and h_freq == 0:
       prepend = 'hp_' + str(l_freq) + '_Hz_'
    else:
       prepend = ''
    
    return prepend

       

#%% MEEG FUNCTIONS


def notch_filter(subject, date, split_recording, input_file, output_file,
                 overwrite):
    raw_path = get_raw_path(subject, date, split_recording)
    load_path = join(raw_path, input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        if not split_recording:
            raw = mne.io.read_raw_fif(load_path, preload=True)
        else:
            raws = []
            full_date = date + '_000000' 
            first_path = join(project_path, 'raw', subject, full_date, 'MEG',
                    '001.tactile_jitter_raw_1', 'files',
                    'tactile_jitter_raw_1.fif')
            raws.append(mne.io.read_raw_fif(first_path, preload=True))
            second_path = join(project_path, 'raw', subject, full_date, 'MEG',
                    '002.tactile_jitter_raw_2', 'files',
                    'tactile_jitter_raw_2.fif')
            raws.append(mne.io.read_raw_fif(second_path, preload=True))
            raw = mne.concatenate_raws(raws)
        
        raw.notch_filter(freqs=np.arange(50, 251, 50))
        raw.save(save_path, overwrite=overwrite)
        
def find_events(subject, date, input_file, output_file, overwrite):
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        raw = mne.io.read_raw_fif(load_path)
        events = mne.find_events(raw, 'STI101', min_duration=0.002)
        
        mne.write_events(save_path, events)
        
def high_and_lowpass(subject, date, input_file, output_file, overwrite, l_freq,
                     h_freq):
    load_path = join(get_save_path(subject, date), input_file)
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        raw = mne.io.read_raw_fif(load_path, preload=True)
        raw.filter(l_freq, h_freq)
        raw.save(save_path, overwrite=overwrite)
        
def epoch_data(subject, date, input_file, output_file, overwrite, h_freq,
               l_freq, events_file):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    events_path = join(get_save_path(subject, date), events_file)
    if should_we_run(save_path, overwrite):
        raw = mne.io.read_raw_fif(load_path, preload=True)
        events = mne.read_events(events_path)
        event_id = dict(s1=1, s2=3, s3=5,
                        s4_0=23, s5_0=25, s6_0=27,
                        s4_5_be=71, s5_5_be=73, s6_5_be=75, 
                        s4_5_af=87, s5_5_af=89, s6_5_af=91,
                        s4_15_be=81, s5_15_be=83, s6_15_be=85, 
                        s4_15_af=97, s5_15_af=99, s6_15_af=101,
                        o0=18, o5=28, o15=38,
                        n1_0=48, n2_0=50, n3_0=52, n4_0=54, n5_0=56,
                        n1_5=58, n2_5=60, n3_5=62, n4_5=64, n5_5=66,
                        n1_15=68, n2_15=70, n3_15=72, n4_15=74, n5_15=76
                        )
        tmin = -0.200
        tmax = 0.600
        baseline = (None, 0)
        decim = 1
        
        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline,
                            proj=False, decim=decim)
        epochs.save(save_path, overwrite=overwrite)
        
def get_evokeds(subject, date, input_file, output_file, overwrite, h_freq,
                l_freq):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        epochs = mne.read_epochs(load_path, proj=False)
        evokeds = []
        for event in epochs.event_id:
            evokeds.append(epochs[event].average())
        mne.write_evokeds(save_path, evokeds)
    
        
def watershed(subject):
    mne.bem.make_watershed_bem(subject, subjects_dir)
    
def source_space_and_bem_model_and_bem_solution(subject):
    ## source space
    print('lala')
    spacing = 'oct6'
    filename = subject + '-' + spacing[:3] + '-' + spacing[-1] + '-src.fif'
    bem_folder = join(subjects_dir, subject, 'bem')
    fullpath = join(bem_folder, filename)
    
    src = mne.source_space.setup_source_space(subject, spacing,
                                          subjects_dir=subjects_dir,
                                          n_jobs=-1, surface='white')
    
    mne.source_space.write_source_spaces(fullpath, src)
    
    ## bem model
    ico = 4
    ico_string = get_ico_string(ico)
    filename = subject + '-' + ico_string + '-bem.fif'
    fullpath = join(bem_folder, filename)
    
    surfaces = mne.bem.make_bem_model(subject, ico, conductivity=[0.3],
                                          subjects_dir=subjects_dir)
    mne.bem.write_bem_surfaces(fullpath, surfaces)
       
    ## bem solution
    filename = subject + '-' + ico_string + '-bem-sol.fif'
    fullpath = join(bem_folder, filename)
        
    bem = mne.bem.make_bem_solution(surfaces)
    mne.bem.write_bem_solution(fullpath, bem)
    
def forward_solution(subject, date, input_file, output_file, overwrite,
                     spacing, ico):
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    bem_folder = join(subjects_dir, subject, 'bem')

    if should_we_run(save_path, overwrite):
        ico_string = get_ico_string(ico)
        info = mne.io.read_info(load_path)
        trans = join(bem_folder, subject + '-trans.fif')
        src = join(bem_folder, subject + '-' + spacing[:3] + '-' + \
                   spacing[-1] + '-src.fif')
        bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')
        
            
        fwd = mne.make_forward_solution(info, trans, src, bem)
        mne.write_forward_solution(save_path, fwd, overwrite)
        
def estimate_covariance(subject, date, input_file, output_file, overwrite,
                        l_freq, h_freq):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        epochs = mne.read_epochs(load_path, proj=False)
        epochs.del_proj()
        cov = mne.compute_covariance(epochs, tmax=0)
        cov.save(save_path)

def inverse_operator(subject, date, input_file_1, input_file_2, input_file_3,
                     output_file, overwrite, h_freq, l_freq):
    load_path_1 = join(get_save_path(subject, date), input_file_1)
    load_path_2 = join(get_save_path(subject, date), input_file_2)
    input_file_3 = hp_lp_string(h_freq, l_freq) + input_file_3
    load_path_3 = join(get_save_path(subject, date), input_file_3)
   
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    save_path = join(get_save_path(subject, date), output_file)
    
    if should_we_run(save_path, overwrite):
        info = mne.io.read_info(load_path_1)
        fwd = mne.read_forward_solution(load_path_2)
        cov = mne.read_cov(load_path_3)
        
        inv = mne.minimum_norm.make_inverse_operator(info, fwd, cov)
        mne.minimum_norm.write_inverse_operator(save_path, inv)
 

def minimum_norm_estimate(subject, date, input_file_1, input_file_2,
                          output_file, overwrite, h_freq, l_freq):
    input_file_1 = hp_lp_string(h_freq, l_freq) + input_file_1
    input_file_2 = hp_lp_string(h_freq, l_freq) + input_file_2
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    org_name = output_file
    load_path_1 = join(get_save_path(subject, date), input_file_1)
    load_path_2 = join(get_save_path(subject, date), input_file_2)
    
    evokeds = mne.read_evokeds(load_path_1, proj=False)
    inv = mne.minimum_norm.read_inverse_operator(load_path_2)
        
    for evoked in evokeds:
        event = evoked.comment
        output_file = ''
        output_file += event + '_' + org_name
        save_path = join(get_save_path(subject, date), 'stcs', output_file)
        if should_we_run(save_path, overwrite):
            stc = mne.minimum_norm.apply_inverse(evoked, inv)
            stc.save(save_path)


    
#%% FUNCTION CALLS


if function_to_run == 'notch_filter':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    notch_filter(subject, date, split_recording, input_file, output_file,
                 overwrite)
    
if function_to_run == 'find_events':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    find_events(subject, date, input_file, output_file, overwrite)
    
if function_to_run == 'high_and_lowpass':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])
    high_and_lowpass(subject, date, input_file, output_file, overwrite,
                     h_freq, l_freq)
    
if function_to_run == 'epoch_data':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])
    events_file = argv[9]
    epoch_data(subject, date, input_file, output_file, overwrite,
               h_freq, l_freq, events_file)

if function_to_run == 'get_evokeds':  
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])
    get_evokeds(subject, date, input_file, output_file, overwrite,
                     h_freq, l_freq)    
    
if function_to_run == 'watershed':
    subject = argv[2]
    watershed(subject)

if function_to_run == 'source_space_and_bem_model_and_bem_solution':
    subject = argv[2]
    source_space_and_bem_model_and_bem_solution(subject)    
    
if function_to_run == 'forward_solution':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    spacing = argv[7]
    ico = int(argv[8])
    forward_solution(subject, date, input_file, output_file, overwrite,
                     spacing, ico)
    
if function_to_run == 'estimate_covariance':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])
    estimate_covariance(subject, date, input_file, output_file, overwrite,
                        h_freq, l_freq)
    
if function_to_run == 'inverse_operator':
    subject = argv[2]
    date = argv[3]
    input_file_1 = argv[4]
    input_file_2 = argv[5]
    input_file_3 = argv[6]
    output_file = argv[7]
    overwrite = bool(int(argv[8]))
    h_freq = int(argv[9])
    l_freq = int(argv[10])
  
    inverse_operator(subject, date, input_file_1, input_file_2, input_file_3,
                     output_file, overwrite, h_freq, l_freq)
    
if function_to_run == 'minimum_norm_estimate':
    subject = argv[2]
    date = argv[3]
    input_file_1 = argv[4]
    input_file_2 = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freq = int(argv[8])
    l_freq = int(argv[9])
    minimum_norm_estimate(subject, date, input_file_1, input_file_2,
                          output_file, overwrite, h_freq, l_freq)   