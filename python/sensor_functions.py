#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:12:37 2020

@author: lau
"""

#%% IMPORTS

import mne
import numpy as np
from scipy import stats as stats
from scipy import linalg
from subprocess import check_output
from os.path import join, isfile, dirname
from os import listdir
from sys import argv
import matplotlib.pyplot as plt
import warnings

#%% FIND RELEVANT FUNCTION AND USER

user = check_output('uname -n', shell=True)
if user == b'lau\n':
    project_path = '/home/lau/projects/cerebellar_clock'
    function_to_run = None
elif user[:6] == b'hyades':
    project_path = '/projects/MINDLAB2019_MEG-CerebellarClock'
    if len(argv) > 1: ## arguments has been supplied from bash
        function_to_run = argv[1]
    else:
        function_to_run = None

print(user)
# print(function_to_run)
# print(argv[:])

#%% SUBJECTS

subjects = [
    '0001', '0002', '0004', '0005', '0006', '0007',
    '0008', '0009', '0010', '0011', '0012', '0013',
    '0014', '0015', '0016', '0017', '0018', '0019',
    '0020', '0021', '0022', '0023', '0024', '0025',
    '0026', '0027', '0028', '0029', '0030', '0031',
    ]
dates = [
    '20200124', '20200121', '20200121', '20200128', '20200120', '20200120',
    '20200120', '20200121', '20200122', '20200122', '20200122', '20200124',
    '20200127', '20200128', '20200128', '20200129', '20200129', '20200129',
    '20200131', '20200131', '20200131', '20200203', '20200203', '20200203',
    '20200204', '20200204', '20200204', '20200205', '20200205', '20200205'
         ]

split_recordings = [
                    0, 0, 0, 0, 1, 0,
                    0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0
                    ]

n_subjects = len(subjects)
proj_name = 'MINDLAB2019_MEG-CerebellarClock'

all_events = []
collapsed_events = ['s1', 's2', 's3', 's4_0', 's5_0', 's6_0',
                    'o0', 'o5', 'o15', 'ns', 's4_5', 's4_15',
                    's5_5', 's5_15', 's6_5', 's6_15']

#%% PATHS

raw_path = join(project_path, 'raw')
scratch_path = join(project_path, 'scratch', 'tactile_jitter', 'MEG')
subjects_dir = join(project_path, 'scratch', 'tactile_jitter', 'freesurfer')
figures_path = join(project_path, 'scratch', 'tactile_jitter', 'figures')

#%% FUNCTION DEFINITIONS
#%% GENERAL HELPER FUNCTIONS

def get_raw_path(subject, date, split_recording):
    full_date = date + '_000000'
    subject_path = join(raw_path, subject, full_date, 'MEG',
                    '001.tactile_jitter_raw', 'files')

    return subject_path

def get_save_path(subject, date):
    if date != '':
        full_date = date + '_000000'
    else:
        full_date = ''
    save_path = join(scratch_path, subject, full_date)

    return save_path

def get_fig_path(subject, date):
    if date != '':
        full_date = date + '_000000'
    else:
        full_date = ''
    save_path = join(figures_path, subject, full_date)
    
    return save_path

def should_we_run(save_path, overwrite, suppress_print=False):
    run_function = overwrite or not isfile(save_path)
    if not run_function:
        if not suppress_print:
            print('not overwriting: ' + save_path)
    if overwrite and isfile(save_path):
        print('Overwriting: ' + save_path)

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
        prepend = 'hp_' + str(h_freq) + '_Hz_'
    elif l_freq > 0 and h_freq == 0:
       prepend = 'lp_' + str(l_freq) + '_Hz_'
    else:
       prepend = ''

    return prepend

def time_string(tmin, tmax):
    string = str(int(tmin * 1e3)) + '_ms_to_' + str(int(tmax * 1e3)) + '_ms_'
    return string

def hp_lp_filter(h_freq, l_freq, raw, n_jobs):
    if h_freq > 0 and l_freq > 0:
        raw.filter(h_freq, l_freq, n_jobs=n_jobs) ## bandpass
    elif h_freq > 0 and l_freq == 0:
        raw.filter(h_freq, None, n_jobs=n_jobs) ## highpass
    elif l_freq > 0 and h_freq == 0:
       raw.filter(None, l_freq, n_jobs=n_jobs)
    else:
       pass

    return raw

def handle_split_recording(subject, date, input_file):
    raws = []
    full_date = date + '_000000'
    if 'lp' in input_file or 'hp' in input_file:
        raw = mne.io.read_raw_fif(join(scratch_path, subject, full_date,
                input_file))
    elif 'sss' in input_file:
        slicer = 8

        first_path = join(scratch_path, subject, full_date,
                input_file[:-slicer] + '_1_sss.fif')
        raws.append(mne.io.read_raw_fif(first_path, preload=True))
        second_path = join(scratch_path, subject, full_date,
                input_file[:-slicer] + '_2_sss.fif')
        raws.append(mne.io.read_raw_fif(second_path, preload=True))
        raw = mne.concatenate_raws(raws)
    else:
        slicer = 4
        first_path = join(project_path, 'raw', subject, full_date, 'MEG',
                '001.tactile_jitter_raw_1', 'files',
                input_file[:-slicer] + '_1.fif')
        raws.append(mne.io.read_raw_fif(first_path, preload=True))
        second_path = join(project_path, 'raw', subject, full_date, 'MEG',
                '002.tactile_jitter_raw_2', 'files',
                input_file[:-slicer] + '_2.fif')
        raws.append(mne.io.read_raw_fif(second_path, preload=True))
        raw = mne.concatenate_raws(raws)

    return raw

def handle_split_recording_info(subject, date, input_file):
    full_date = date + '_000000'
    if 'lp' in input_file or 'hp' in input_file:
        info = mne.io.read_info(join(scratch_path, subject, full_date,
                input_file))
    elif 'sss' in input_file:
        slicer = 8

        first_path = join(scratch_path, subject, full_date,
                input_file[:-slicer] + '_1_sss.fif')
        info = mne.io.read_info(first_path)
        # second_path = join(scratch_path, subject, full_date,
        #         input_file[:-slicer] + '_2_sss.fif')
        # raws.append(mne.io.read_raw_fif(second_path, preload=True))
        # raw = mne.concatenate_raws(raws)
    else:
        slicer = 4
        first_path = join(project_path, 'raw', subject, full_date, 'MEG',
                '001.tactile_jitter_raw_1', 'files',
                input_file[:-slicer] + '_1.fif')
        info = mne.io.read_info(first_path)
        # second_path = join(project_path, 'raw', subject, full_date, 'MEG',
        #         '002.tactile_jitter_raw_2', 'files',
        #         input_file[:-slicer] + '_2.fif')
        # raws.append(mne.io.read_raw_fif(second_path, preload=True))
        # raw = mne.concatenate_raws(raws)

    return info

def handle_split_recording_fwd(subject, date, input_file):
    full_date = date + '_000000'
    if 'sss' in input_file:
        slicer = 8
        path = join(scratch_path, subject, full_date,
                    input_file[:-slicer] + '_1_sss.fif')
    else:
        slicer = 4
        path = join(project_path, 'raw', subject, full_date, 'MEG',
                    '001.tactile_jitter_raw_1', 'files',
                    input_file[:-slicer] + '_1.fif')
    return path

def collapse_conditions(epochs, ISI_type, baseline=False, tmin=-0.400):
    if baseline != 'no_baseline_Ss':
        to_collapse = {}
        to_collapse['ns']    = ['n1_0',  'n2_0',  'n3_0',  'n4_0',  'n5_0',
                                'n1_5',  'n2_5',  'n3_5',  'n4_5',  'n5_5',
                                'n1_15', 'n2_15', 'n3_15', 'n4_15', 'n5_15'
                                ]
        if ISI_type != 'local' and ISI_type != 'mean' and not baseline and \
            tmin == -0.400:
            to_collapse['s4_5']  = ['s4_5_be',  's4_5_af']
            to_collapse['s4_15'] = ['s4_15_be', 's4_15_af']
            to_collapse['s5_5']  = ['s5_5_be',  's5_5_af']
            to_collapse['s5_15'] = ['s5_15_be', 's5_15_af']
            to_collapse['s6_5']  = ['s6_5_be',  's6_5_af']
            to_collapse['s6_15'] = ['s6_15_be', 's6_15_af']
        new_triggers = [200, 33, 35, 37, 43, 45, 47]
    
        for cond_index, cond in enumerate(to_collapse):
                new_dict = {cond:new_triggers[cond_index]}
                mne.epochs.combine_event_ids(epochs, to_collapse[cond],
                                             new_dict,
                                             copy=False)
    return epochs

def freq_time_string_contrast(contrast, fmin, fmax, tmin, tmax):
    string = contrast[0] + '_' + contrast[1] + '_contrast_' + \
             str(fmin) + '_' + str(fmax) + '_Hz_' + \
             str(int(1e3 * tmin)) + '_' + str(int(1e3 * tmax)) + \
             '_ms_'
    return string

def freq_time_string_condition(condition, fmin, fmax, tmin, tmax):
    string = condition + '_' + str(fmin) + '_' + str(fmax) + '_Hz_' + \
             str(int(1e3 * tmin)) + '_' + str(int(1e3 * tmax)) + \
             '_ms_'
    return string

def spacing_contrast_string(spacing, contrast):
    string = contrast[0] + '_' + contrast[1] + '_contrast_' + str(spacing) + \
        '_mm_'
    return string

def spacing_condition_string(spacing, condition):
    string = condition + '_' + str(spacing) + '_mm_'
    return string

def spacing_condition_common_filter_string(spacing, condition, contrast):
    contrast_part = ''
    for part in contrast:
        contrast_part += part + '_'
    string = condition + '_' + str(spacing) + '_mm_common_filter_' + \
             contrast_part
    return string

def reg_and_weight_norm_string(reg, weight_norm):
    string = 'reg_' + str(reg) + '_weight_' + weight_norm + '_'
    return string

def ISI_and_baseline_string(ISI_type, baseline):
    string = ISI_type + '_ISI_' + baseline + '_'
    return string
    

def list_parser(string):

    List = []
    new_string = ''
    for character in string:
        if character != '[' and character != ']' and character != ',' and \
            character != ' ':
            new_string += character
        if character == ',' or character == ']':
            List.append(new_string)
            new_string = ''
    return List

def list_of_list_parser(string):
    list_of_lists = []
    first_level_parse = list_parser(string)
    sublist = [] # initialize
    for element in first_level_parse:
        if element == '':
            list_of_lists.append(sublist)
            sublist = []
        else:
            sublist.append(element)
    return list_of_lists

def compute_source_itc(stcs):
    n_trials = len(stcs)

    tmp = np.zeros(stcs[0].data.shape, dtype=np.complex)
    for stc in stcs:
        # divide by amplitude and sum angles
        tmp += stc.data / abs(stc.data)

    # take absolute value and normalize
    itc = abs(tmp) / n_trials

    return itc

def get_load_and_save_paths_for_common_filter_hilbert(subject, date,
                                                      spacing, contrast,
                                                      h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output,
                                                      reg, weight_norm,
                                                      channel_type,
                                                      grand_average,
                                                      statistics,
                                                      p_threshold,
                                                      tmin, tmax,
                                                      ISI_type,
                                                      figures=False,
                                                      stat_function=None,
                                                      proj=False,
                                                      baseline='no_baseline',
                                                      stat_tmin=None,
                                                      stat_tmax=None,
                                                      save_single_stcs=False):
    
    if grand_average and statistics:
        raise RuntimeError('"grand_average" and "statistics" cannot' + \
                           'both be true')
    if weight_norm is None:
        weight_norm = 'None'
        
    if channel_type == 'both':
        channel_type = 'mag_and_grad'
    if proj:
        proj = 'proj_'
    else:
        proj = ''
        
    ISI_string = ISI_type + '_ISI_'

    baseline_text = baseline + '_'

    if stat_tmin is not None and stat_tmax is not None:
        time_str = 'stat_' + time_string(stat_tmin, stat_tmax) + \
                   'range_' + time_string(tmin, tmax)
    else:
         time_str = time_string(tmin, tmax)
        
    
    input_file_1 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[0],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj + \
                                    org_name_input
    input_file_2 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[1],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj +  \
                                    org_name_input
    if len(contrast) == 2:
        input_file_3 = \
            ISI_string + baseline_text + channel_type + '_' + \
                reg_and_weight_norm_string(reg, weight_norm) + \
                spacing_contrast_string(spacing, contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                    time_str + \
                                        proj + org_name_input
    elif len(contrast) == 3:
        input_file_3 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[2],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj + \
                                    org_name_input
        
    output_file_1 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[0],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj + \
                                    org_name_output
    output_file_2 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[1],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj + \
                                    org_name_output
    if len(contrast) == 2:
        output_file_3 = \
            ISI_string + baseline_text + channel_type + '_' + \
                reg_and_weight_norm_string(reg, weight_norm) + \
                spacing_contrast_string(spacing, contrast) + \
                            hp_lp_string(h_freq, l_freq)  + \
                                    time_str + proj + \
                                        org_name_output
    elif len(contrast) == 3:
        output_file_3 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[2],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_str + proj + \
                                    org_name_output
                                        
                        
    if grand_average:
        subject = 'grand_averages'
        date = ''   
    if statistics:
        subject = 'grand_averages/statistics'
        date = ''
        input_file_1  = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        input_file_1
        input_file_2  = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        input_file_2
        input_file_3  = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        input_file_3
        output_file_1 = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        output_file_1
        output_file_2 = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        output_file_2
        output_file_3 = stat_function + '_' + 'p_' + str(p_threshold) + '_' + \
                        output_file_3
                
    load_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                       'common_filter', input_file_1)
    load_path_2 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                       'common_filter', input_file_2)
    load_path_3 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                       'common_filter', input_file_1)
    load_path_4 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                       'common_filter', input_file_2)
    if len(contrast) == 2:
        load_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', input_file_3)
        load_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', 'contrasts', input_file_3)
    elif len(contrast) == 3:
        load_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', input_file_3)
        load_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', input_file_3)
        
    if save_single_stcs:
        load_path_1 = join(get_save_path(subject, date), 'hilbert',
                           'connectivity', 'stcs',
                            input_file_1)
        load_path_2 = join(get_save_path(subject, date), 'hilbert',
                           'connectivity', 'stcs',
                           input_file_2)
        load_path_3 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', input_file_1)
        load_path_4 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', input_file_2)
        load_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', input_file_3)
        load_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', 'contrasts', input_file_3)
     
        
    if figures:                   
        save_path_1 = join(get_fig_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_1)
        save_path_2 = join(get_fig_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_2)
        save_path_3 = join(get_fig_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_1)
        save_path_4 = join(get_fig_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_2)
        if len(contrast) == 2:
            save_path_5 = join(get_fig_path(subject, date), 'hilbert', 'stcs',
                               'common_filter', 'contrasts', output_file_3)
            save_path_6 = join(get_fig_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', 'contrasts', output_file_3)
        elif len(contrast) == 3:
            save_path_5 = join(get_fig_path(subject, date), 'hilbert', 'stcs',
                               'common_filter', output_file_3)
            save_path_6 = join(get_fig_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', output_file_3)
    elif save_single_stcs:
        save_path_1 = join(get_save_path(subject, date), 'hilbert',
                           'connectivity', 'stcs', output_file_1)
        save_path_2 = join(get_save_path(subject, date), 'hilbert',
                           'connectivity', 'stcs', output_file_2)
        save_path_3 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_1)
        save_path_4 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_2)
        if len(contrast) == 2:
            save_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', output_file_3)
            save_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', 'contrasts', output_file_3)
        elif len(contrast) == 3:
            save_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                       'common_filter', output_file_3)
            save_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', output_file_3)
        
    else:
        save_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_1)
        save_path_2 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_2)
        save_path_3 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_1)
        save_path_4 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_2)
        if len(contrast) == 2:
            save_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', output_file_3)
            save_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', 'contrasts', output_file_3)
        elif len(contrast) == 3:
            save_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                       'common_filter', output_file_3)
            save_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                               'common_filter', output_file_3)
                
                
    
    load_paths = [load_path_1, load_path_2, load_path_3,
                  load_path_4, load_path_5, load_path_6]
    save_paths = [save_path_1, save_path_2, save_path_3,
                      save_path_4, save_path_5, save_path_6]
    paths = [load_paths, save_paths]
    
    return paths

def wilcoxon(x, y=None, zero_method="wilcox", correction=False,
             alternative="two-sided"):
    """
    Adapted by Lau to spit out z-values
    Calculate the Wilcoxon signed-rank test.

    The Wilcoxon signed-rank test tests the null hypothesis that two
    related paired samples come from the same distribution. In particular,
    it tests whether the distribution of the differences x - y is symmetric
    about zero. It is a non-parametric version of the paired T-test.

    Parameters
    ----------
    x : array_like
        Either the first set of measurements (in which case `y` is the second
        set of measurements), or the differences between two sets of
        measurements (in which case `y` is not to be specified.)  Must be
        one-dimensional.
    y : array_like, optional
        Either the second set of measurements (if `x` is the first set of
        measurements), or not specified (if `x` is the differences between
        two sets of measurements.)  Must be one-dimensional.
    zero_method : {'pratt', 'wilcox', 'zsplit'}, optional
        The following options are available (default is 'wilcox'):
     
          * 'pratt': Includes zero-differences in the ranking process,
            but drops the ranks of the zeros, see [4]_, (more conservative).
          * 'wilcox': Discards all zero-differences, the default.
          * 'zsplit': Includes zero-differences in the ranking process and 
            split the zero rank between positive and negative ones.
    correction : bool, optional
        If True, apply continuity correction by adjusting the Wilcoxon rank
        statistic by 0.5 towards the mean value when computing the
        z-statistic.  Default is False.
    alternative : {"two-sided", "greater", "less"}, optional
        The alternative hypothesis to be tested, see Notes. Default is
        "two-sided".

    Returns
    -------
    statistic : float
        If `alternative` is "two-sided", the sum of the ranks of the
        differences above or below zero, whichever is smaller.
        Otherwise the sum of the ranks of the differences above zero.
    pvalue : float
        The p-value for the test depending on `alternative`.

    See Also
    --------
    kruskal, mannwhitneyu

    Notes
    -----
    The test has been introduced in [4]_. Given n independent samples
    (xi, yi) from a bivariate distribution (i.e. paired samples),
    it computes the differences di = xi - yi. One assumption of the test
    is that the differences are symmetric, see [2]_.
    The two-sided test has the null hypothesis that the median of the
    differences is zero against the alternative that it is different from
    zero. The one-sided test has the null hypothesis that the median is 
    positive against the alternative that it is negative 
    (``alternative == 'less'``), or vice versa (``alternative == 'greater.'``).

    The test uses a normal approximation to derive the p-value (if
    ``zero_method == 'pratt'``, the approximation is adjusted as in [5]_).
    A typical rule is to require that n > 20 ([2]_, p. 383). For smaller n,
    exact tables can be used to find critical values.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test
    .. [2] Conover, W.J., Practical Nonparametric Statistics, 1971.
    .. [3] Pratt, J.W., Remarks on Zeros and Ties in the Wilcoxon Signed
       Rank Procedures, Journal of the American Statistical Association,
       Vol. 54, 1959, pp. 655-667. :doi:`10.1080/01621459.1959.10501526`
    .. [4] Wilcoxon, F., Individual Comparisons by Ranking Methods,
       Biometrics Bulletin, Vol. 1, 1945, pp. 80-83. :doi:`10.2307/3001968`
    .. [5] Cureton, E.E., The Normal Approximation to the Signed-Rank
       Sampling Distribution When Zero Differences are Present,
       Journal of the American Statistical Association, Vol. 62, 1967,
       pp. 1068-1069. :doi:`10.1080/01621459.1967.10500917`

    Examples
    --------
    In [4]_, the differences in height between cross- and self-fertilized
    corn plants is given as follows:

    >>> d = [6, 8, 14, 16, 23, 24, 28, 29, 41, -48, 49, 56, 60, -67, 75]

    Cross-fertilized plants appear to be be higher. To test the null
    hypothesis that there is no height difference, we can apply the
    two-sided test:

    >>> from scipy.stats import wilcoxon
    >>> w, p = wilcoxon(d)
    >>> w, p
    (24.0, 0.04088813291185591)

    Hence, we would reject the null hypothesis at a confidence level of 5%,
    concluding that there is a difference in height between the groups.
    To confirm that the median of the differences can be assumed to be
    positive, we use:

    >>> w, p = wilcoxon(d, alternative='greater')
    >>> w, p
    (96.0, 0.020444066455927955)

    This shows that the null hypothesis that the median is negative can be
    rejected at a confidence level of 5% in favor of the alternative that
    the median is greater than zero. The p-value based on the approximation
    is within the range of 0.019 and 0.054 given in [2]_.
    Note that the statistic changed to 96 in the one-sided case (the sum
    of ranks of positive differences) whereas it is 24 in the two-sided
    case (the minimum of sum of ranks above and below zero).

    """

    if zero_method not in ["wilcox", "pratt", "zsplit"]:
        raise ValueError("Zero method should be either 'wilcox' "
                         "or 'pratt' or 'zsplit'")

    if alternative not in ["two-sided", "less", "greater"]:
        raise ValueError("Alternative must be either 'two-sided', "
                         "'greater' or 'less'")

    if y is None:
        d = np.asarray(x)
        if d.ndim > 1:
            raise ValueError('Sample x must be one-dimensional.')
    else:
        x, y = map(np.asarray, (x, y))
        if x.ndim > 1 or y.ndim > 1:
            raise ValueError('Samples x and y must be one-dimensional.')
        if len(x) != len(y):
            raise ValueError('The samples x and y must have the same length.')
        d = x - y

    if zero_method in ["wilcox", "pratt"]:
        n_zero = np.sum(d == 0, axis=0)
        if n_zero == len(d):
            raise ValueError("zero_method 'wilcox' and 'pratt' do not work if "
                             "the x - y is zero for all elements.")

    if zero_method == "wilcox":
        # Keep all non-zero differences
        d = np.compress(np.not_equal(d, 0), d, axis=-1)

    count = len(d)
    if count < 10:
        warnings.warn("Sample size too small for normal approximation.")
        
    r = stats.rankdata(abs(d))
    r_plus = np.sum((d > 0) * r, axis=0) ## FIXME: real part is compared by np
    r_minus = np.sum((d < 0) * r, axis=0)

    if zero_method == "zsplit":
        r_zero = np.sum((d == 0) * r, axis=0)
        r_plus += r_zero / 2.
        r_minus += r_zero / 2.

    # return min for two-sided test, but r_plus for one-sided test
    # the literature is not consistent here
    # r_plus is more informative since r_plus + r_minus = count*(count+1)/2,
    # i.e. the sum of the ranks, so r_minus and the min can be inferred
    # (If alternative='pratt', r_plus + r_minus = count*(count+1)/2 - r_zero.)
    # [3] uses the r_plus for the one-sided test, keep min for two-sided test
    # to keep backwards compatibility
    if alternative == "two-sided":
        T = min(r_plus, r_minus)
    else:
        T = r_plus
    mn = count * (count + 1.) * 0.25
    se = count * (count + 1.) * (2. * count + 1.)

    if zero_method == "pratt":
        r = r[d != 0]
        # normal approximation needs to be adjusted, see Cureton (1967)
        mn -= n_zero * (n_zero + 1.) * 0.25
        se -= n_zero * (n_zero + 1.) * (2. * n_zero + 1.)

    replist, repnum = stats.find_repeats(r)
    if repnum.size != 0:
        # Correction for repeated elements.
        se -= 0.5 * (repnum * (repnum * repnum - 1)).sum()

    se = np.sqrt(se / 24)

    # apply continuity correction if applicable
    d = 0
    if correction:
        if alternative == "two-sided":
            d = 0.5 * np.sign(T - mn)
        elif alternative == "less":
            d = -0.5
        else:
            d = 0.5

    # compute statistic and p-value using normal approximation
    z = (T - mn - d) / se
    ## added by lau, necessary to get right sign of z, which is not 
    ## calculated in the scipy implementation
    if r_plus > r_minus:
        z = np.abs(z)

    return z

def get_max_vertices_per_roi(atlas, rois, subject, date, proj, h_freq, l_freq,
                             tmin, tmax, spacing, reg, weight_norm,
                             channel_type, contrast, ISI_type, baseline):
    
    from nilearn import image, datasets
    
    filepath = join(scratch_path, get_save_path(subject, date), 'hilbert',
                    'stcs', 'common_filter', 'contrasts')
    filename = ISI_and_baseline_string(ISI_type, baseline) + channel_type + \
               '_' + reg_and_weight_norm_string(reg, weight_norm) + \
               spacing_contrast_string(spacing, contrast) + \
               hp_lp_string(h_freq, l_freq) + time_string(tmin, tmax) + \
               'tactile_jitter_morph-vl.nii'
    img = image.load_img(join(filepath, filename))
    data = np.asanyarray(img.dataobj)
    
    src = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                      'volume-' + str(spacing) + 'mm-src.fif'))
    vertices = src[0]['vertno']

    if atlas == 'aal':
        atlas_dataset = datasets.fetch_atlas_aal()
        atlas_filepath = atlas_dataset['maps']
        atlas_img = image.load_img(atlas_filepath)
        labels = atlas_dataset['labels']
        indices = list(map(int, atlas_dataset['indices']))
    elif atlas == 'harvard_oxford':
        atlas_filename = 'cort_maxprob_thr0_1mm'
        atlas_datasets = datasets.fetch_atlas_harvard_oxford(atlas_filename)
        atlas_filepath = atlas_datasets.maps
        atlas_img = image.load_img(atlas_filepath)
        labels = atlas_dataset['labels']
        indices = list(map(int, np.unique(atlas_img.get_fdata())))
        
    ## resampling
    atlas_interpolated = image.resample_to_img(atlas_img, img, 'nearest')
    atlas_interpolated_data = np.asanyarray(atlas_interpolated.dataobj)
    
    ## find peak voxels
    n_labels = len(labels)
    
    if atlas == 'aal':
        ## identify max voxel
        max_coordinates = np.zeros((n_labels, 3))
        for label_index, label in enumerate(indices):
            this_mask = atlas_interpolated_data == label
            label_data = data[this_mask, :]
            argmax = np.argmax(np.max(label_data, axis=1))
            xs, ys, zs = np.where(this_mask)
            max_coordinates[label_index] = xs[argmax], ys[argmax], zs[argmax]
            
        ## find vertices
        stc_voxels = np.array(
            np.unravel_index(vertices, img.shape[:3], order='F')).T
        vert_nos = [None] * n_labels
        n_xs, n_ys, n_zs = src[0]['shape']
        
        for label_index, label in enumerate(labels):
            max_coordinate = max_coordinates[label_index]
            distances = np.linalg.norm(max_coordinate - stc_voxels, axis=1)
            vert_nos[label_index] = vertices[np.argmin(distances)]
            
    if atlas == 'harvard_oxford':
         stc_voxels = np.array(
            np.unravel_index(vertices, img.shape[:3], order='F')).T
         for label_index, label in enumerate(indices):
            if labels[label_index] not in rois:
                continue
            vert_nos[label_index] = []
            
            lefts = []
            rights = []
            this_mask = atlas_interpolated_data == label
            xs, ys, zs = np.where(this_mask)
            distances = np.zeros(len(stc_voxels))
            for coordinate_index, coordinate in enumerate(zip(xs, ys, zs)):
                for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
                    distances[stc_voxel_index] = \
                        np.linalg.norm(coordinate - stc_voxel)
                closest_vertex = vertices[np.argmin(distances)]
                coordinate = src[0]['rr'][closest_vertex]
                if coordinate[0] > 0:
                    rights.append(coordinate_index)
                    
                if coordinate[0] < 0:
                    lefts.append(coordinate_index)
            label_data = data[this_mask, :]
            left_data = label_data[lefts, :]
            right_data = label_data[rights, :]
            
            ## left
            argmax = np.argmax(np.max(left_data, axis=1))
            max_coordinate = xs[lefts][argmax], ys[lefts][argmax], \
                            zs[lefts][argmax]
            distances = np.zeros(len(stc_voxels))
            for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
                distances[stc_voxel_index] = \
                    np.linalg.norm(max_coordinate - stc_voxel)
            vert_nos[label_index].append(vertices[np.argmin(distances)])
            
            ## right
            argmax = np.argmax(np.max(right_data, axis=1))
            max_coordinate = xs[rights][argmax], ys[rights][argmax], \
                            zs[rights][argmax]
            distances = np.zeros(len(stc_voxels))
            for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
                distances[stc_voxel_index] = \
                    np.linalg.norm(max_coordinate - stc_voxel)
            vert_nos[label_index].append(vertices[np.argmin(distances)])

         
    ## map vert nos to rois
    roi_vertices = dict()
    for label_index, label in enumerate(labels):
        if label in rois:
            if atlas == 'aal':
                roi_vertices[label] = vert_nos[label_index]
            if atlas == 'harvard_oxford':
                roi_vertices[label + '_L'] = vert_nos[label_index][0]
                roi_vertices[label + '_R'] = vert_nos[label_index][1]
    
    return roi_vertices
              

#%% MEEG FUNCTIONS

def high_and_lowpass(subject, date, split_recording, input_file, output_file,
                     overwrite, h_freq, l_freq):
    load_path = join(get_save_path(subject, date), input_file)
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        if not split_recording:
            raw = mne.io.read_raw_fif(load_path, preload=True)
        else:
            raw = handle_split_recording(subject, date, input_file)
        raw = hp_lp_filter(h_freq, l_freq, raw)
        raw.save(save_path, overwrite=overwrite)

def find_events(subject, date, input_file, output_file, overwrite):
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    local_save_path = join(get_save_path(subject, date), 'local_' + \
                                                        output_file)
    mean_save_path = join(get_save_path(subject, date), 'mean_' + \
                                                        output_file)
    if should_we_run(save_path, overwrite) or \
            should_we_run(local_save_path, overwrite) or \
            should_we_run(mean_save_path, overwrite):
        raw = mne.io.read_raw_fif(load_path)
        events = mne.find_events(raw, 'STI101', min_duration=0.002)
        mne.write_events(save_path, events)
        
        ## find local and mean ISI events
        omission_indices = [28, 38]
        local_events = events.copy()
        mean_events = events.copy()
        
        for omission_index in omission_indices:
            
            indices = np.where(events[:, 2] == omission_index)[0]
            
            for index in indices:
                s3_index = index - 4
                s4_index = index - 3
                s5_index = index - 2
                s6_index = index - 1
                local_ISI = events[s6_index, 0] - events[s5_index, 0]
                mean_ISI = int(np.mean([events[s4_index, 0] - \
                                            events[s3_index, 0],
                                        events[s5_index, 0] - \
                                            events[s4_index, 0],
                                        events[s6_index, 0] - \
                                            events[s5_index, 0]]))
                local_sample = events[s6_index, 0] + local_ISI
                mean_sample  = events[s6_index, 0] + mean_ISI
                
                local_events[index, 0] = local_sample
                mean_events[index, 0] = mean_sample
        

        mne.write_events(local_save_path, local_events)
        mne.write_events(mean_save_path, mean_events)

        

def epoch_data(subject, date, split_recording, input_file, output_file,
               overwrite, h_freq, l_freq, events_file):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    raw_path = join(get_raw_path(subject, date, split_recording), input_file)
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    events_path = join(get_save_path(subject, date), events_file)
    if should_we_run(save_path, overwrite):
        if not split_recording and 'sss' in input_file:
            raw = mne.io.read_raw_fif(load_path, preload=True)
        elif not split_recording and 'sss' not in input_file:
            raw = mne.io.read_raw_fif(raw_path, preload=True)
        else:
            raw = handle_split_recording(subject, date, input_file)
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
        reject =  dict(grad=4000e-13,
                       mag=4e-12,
                       eog=250e-6)

        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline,
                            proj=False, decim=decim, reject=reject)
        epochs.save(save_path, overwrite=overwrite)

def get_evokeds(subject, date, input_file, output_file, overwrite, h_freq,
                l_freq, collapse_cond):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        epochs = mne.read_epochs(load_path, proj=False)
        epochs.del_proj()
        picks = mne.pick_types(epochs.info, meg=True, eog=True,
                                   ecg=True, emg=True, misc=True)
        if collapse_cond:
            epochs = collapse_conditions(epochs)
        evokeds = []
        for event in epochs.event_id:
            evokeds.append(epochs[event].average(picks=picks))
        mne.write_evokeds(save_path, evokeds)

def grand_average_evokeds(subjects, dates, input_file, output_file, overwrite,
                          h_freq, l_freq):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    save_path = join(get_save_path('grand_averages', ''), output_file)
    if should_we_run(save_path, overwrite):
        evokeds_dict = dict()
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            load_path = join(get_save_path(subject, date), input_file)
            evokeds = mne.read_evokeds(load_path)
            for evoked in evokeds:
                event = evoked.comment
                print('Loading event: ' + event + ' for subject: ' + subject)
                if subject_index == 0:
                    evokeds_dict[event] = []
                evokeds_dict[event].append(evoked)
            evokeds = None ## free up memory

        ga_evokeds = []
        for event in evokeds_dict:
            temp = evokeds_dict[event]
            ga = mne.grand_average(temp)
            ga.comment = event
            ga_evokeds.append(ga)

        mne.evoked.write_evokeds(save_path, ga_evokeds)

def watershed(subject, overwrite):
    mne.bem.make_watershed_bem(subject, subjects_dir, overwrite=overwrite)

def source_space_and_bem_model_and_bem_solution(subject):
    ## source space
    if subject == 'fsaverage':
        spacing = 'ico4'
    else:
        spacing = 'oct6'
    filename = subject + '-' + spacing[:3] + '-' + spacing[-1] + '-src.fif'
    bem_folder = join(subjects_dir, subject, 'bem')
    fullpath = join(bem_folder, filename)

    src = mne.source_space.setup_source_space(subject, spacing,
                                          subjects_dir=subjects_dir,
                                          n_jobs=-1, surface='white')

    mne.source_space.write_source_spaces(fullpath, src, overwrite=True)

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

def volumetric_source_space(subject, spacing):

    bem_folder = join(subjects_dir, subject, 'bem')

    ## bem model
    ico = 4
    ico_string = get_ico_string(ico)
    filename = subject + '-' + ico_string + '-bem.fif'
    fullpath = join(bem_folder, filename)

    src = mne.source_space.setup_volume_source_space(subject, pos=spacing,
                                                     bem=fullpath)
    filename = 'volume-' + str(spacing) +'mm-src.fif'
    fullpath = join(bem_folder, filename)
    mne.source_space.write_source_spaces(fullpath, src, overwrite=True)


def morph_to_fsaverage(subject, date, output_file, overwrite, ico):
    
    output_file = 'ico-' + str(ico) + '_' + output_file
    save_path = join(get_save_path(subject, date), output_file)

    spacing = 'oct6'
    filename = subject + '-' + spacing[:3] + '-' + spacing[-1] + '-src.fif'
    bem_folder = join(subjects_dir, subject, 'bem')
    fullpath = join(bem_folder, filename)
    if should_we_run(save_path, overwrite):
        src = mne.source_space.read_source_spaces(fullpath)
        morph = mne.compute_source_morph(src, subject,
                                         subjects_dir=subjects_dir,
                                         spacing=ico)

        morph.save(save_path, overwrite=overwrite)

def morph_to_fsaverage_volume(subject, date, output_file, overwrite,
                              spacing):
    output_file = str(spacing) + '_mm_' + output_file
    save_path = join(get_save_path(subject, date), output_file)

    bem_folder = join(subjects_dir, subject, 'bem')
    filename = 'volume-' + str(spacing) + 'mm-src.fif'
    src_to = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                         filename)) ## from mne 0.20
    if should_we_run(save_path, overwrite):
        src = join(bem_folder, filename)
        morph = mne.compute_source_morph(src, subject,
                                         subjects_dir=subjects_dir,
                                         src_to=src_to)

        morph.save(save_path, overwrite=overwrite)

def forward_solution(subject, date, split_recording, input_file, output_file,
                     overwrite, spacing, ico):
    if split_recording:
        load_path = handle_split_recording_fwd(subject, date, input_file)
    else:
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

        fwd = mne.make_forward_solution(info, trans, src, bem, n_jobs=-1)
        mne.write_forward_solution(save_path, fwd, overwrite)

def forward_solution_volumetric(subject, date, split_recording, input_file,
                                output_file, overwrite, spacing, ico):
    if split_recording:
        load_path = handle_split_recording_fwd(subject, date, input_file)
    else:
        load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    bem_folder = join(subjects_dir, subject, 'bem')


    if should_we_run(save_path, overwrite):
        ico_string = get_ico_string(ico)
        info = mne.io.read_info(load_path)
        trans = join(bem_folder, subject + '-trans.fif')
        src = join(bem_folder, 'volume-' + str(spacing) + 'mm-src.fif')
        bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')

        fwd = mne.make_forward_solution(info, trans, src, bem, n_jobs=1)
        mne.write_forward_solution(save_path, fwd, overwrite)

def estimate_covariance(subject, date, input_file, output_file, overwrite,
                        h_freq, l_freq):
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        epochs = mne.read_epochs(load_path, proj=False)
        epochs.del_proj()
        if 'sss' in input_file:
            rank = 'info'
        else:
            rank = None
        cov = mne.compute_covariance(epochs, tmax=0, rank=rank)
        cov.save(save_path)

def inverse_operator(subject, date, split_recording, input_file_1,
                     input_file_2, input_file_3, output_file, overwrite,
                     h_freq, l_freq):
    if split_recording:
        load_path_1 = handle_split_recording_fwd(subject, date, input_file_1)
    else:
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

def morph_minimum_norm_estimate(subject, date, input_file, input_file_morph,
                                output_file, overwrite, h_freq, l_freq,
                                collapse_cond, ico):
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    output_file = hp_lp_string(h_freq, l_freq) + output_file + '_ico-' + \
                 str(ico)
    input_file_morph = 'ico-' + str(ico) + '_' + input_file_morph
    org_name_input = input_file
    org_name_output = output_file
    load_path_morph = join(get_save_path(subject, date), input_file_morph)
    morph = mne.read_source_morph(load_path_morph)
    if collapse_cond:
        these_events = collapsed_events
    else:
        these_events = all_events
    for event in these_events:
        input_file = event + '_' + org_name_input
        output_file = event + '_' + org_name_output
        load_path = join(get_save_path(subject, date), 'stcs', input_file)
        save_path = join(get_save_path(subject, date), 'stcs', output_file)
        if should_we_run(save_path, overwrite):
            stc = mne.read_source_estimate(load_path)
            stc_morph = morph.apply(stc)
            stc_morph.save(save_path)

def grand_average_stc(subjects, dates, input_file, output_file, overwrite,
                      h_freq, l_freq, collapse_cond, ico):
    input_file = hp_lp_string(h_freq, l_freq) + input_file + '_ico-' + str(ico)
    output_file = hp_lp_string(h_freq, l_freq) + output_file + '_ico-' + \
                  str(ico)
    org_name_input = input_file
    org_name_output = output_file
    save_path = join(get_save_path('grand_averages', 'stcs'), output_file)
    if collapse_cond:
        these_events = collapsed_events
    else:
        these_events = all_events
    if should_we_run(save_path, overwrite):
        stcs_dict = dict()
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            for event in these_events:
                input_file = event + '_' + org_name_input
                load_path = join(get_save_path(subject, date), 'stcs',
                                 input_file)
                print('Loading: ' + event + ' for subject: ' + subject)
                stc = mne.read_source_estimate(load_path)
                if subject_index == 0:
                    stcs_dict[event] = []
                stcs_dict[event].append(stc)

        ga_stcs = dict()
        for event in stcs_dict:
            stcs = stcs_dict[event]
            for stc_index, stc in enumerate(stcs):
                if stc_index == 0:
                    temp = stc.copy()
                else:
                    temp._data += stc.data
            temp._data /= n_subjects
            ga_stcs[event] = temp

        for event in ga_stcs:
            output_file = event + '_' + org_name_output
            save_path = join(get_save_path('grand_averages', ''), 'stcs',
                             output_file)
            stc = ga_stcs[event]
            stc.save(save_path)
            
def statistics_stc(subjects, dates, input_file, output_file, overwrite, h_freq,
                   l_freq, contrasts, p_threshold, tmin, tmax, ico, n_jobs):
    input_file = hp_lp_string(h_freq, l_freq) + input_file + '_ico-' + str(ico)
    output_file = hp_lp_string(h_freq, l_freq) + output_file + '_ico-' + \
                  str(ico)
    org_name_input = input_file
    org_name_output = output_file
    n_subjects = len(subjects)
    
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'fsaverage-ico-' + \
                                                     str(ico) +'-src.fif'))
    
    ## gather data
    for contrast in contrasts:
        contrast_name = contrast[0] + '_vs_' + contrast[1]
        print('Running contrast: ' + contrast_name)
        output_file = contrast_name + '_' + org_name_output
        save_path = join(get_save_path('grand_averages', ''), 'statistics',
                         'stcs', output_file)
        stat_arrays = [None] * len(contrast)
        if should_we_run(save_path, overwrite):
            for event_index, event in enumerate(contrast):
                input_file = event + '_' + org_name_input
                for subject_index, subject in enumerate(subjects):
                    date = dates[subject_index]
                    load_path = join(get_save_path(subject, date), 'stcs',
                                     input_file)
                    print('Loading: ' + input_file + ' for subject: ' + subject)
                    stc = mne.read_source_estimate(load_path)
                    stc.crop(tmin, tmax)
                    if subject_index == 0:
                        n_times = stc.data.shape[1]
                        n_vertices = stc.data.shape[0]
                        stat_arrays[event_index] = np.zeros((n_subjects,
                                                             n_times,
                                                             n_vertices))
                    stat_arrays[event_index][subject_index, :, :] = stc.data.T
            ## do statistics
            connectivity = mne.spatio_temporal_src_connectivity(src, n_times)
            t_threshold = -stats.distributions.t.ppf(p_threshold / 2.0,
                                                     n_subjects - 1)
            t_obs, clusters, cluster_p_values, H0 = clu = \
                mne.stats.spatio_temporal_cluster_test(stat_arrays,
                    connectivity=connectivity, n_jobs=n_jobs,
                    threshold=t_threshold, verbose=True, seed=7,
                    n_permutations=1024)
            print('Saving contrast: ' + contrast_name)
            cluster_dict = dict()
            names = ['t_obs', 'clusters', 'cluster_p_values', 'H0']
            for name_index, name in enumerate(names):
                cluster_dict[name] = clu[name_index]
            np.save(save_path, cluster_dict)
                       


def filter_for_hilbert(subject, date, split_recording, input_file, 
                       output_file, overwrite, h_freqs, l_freqs, n_jobs):
    if 'sss' in input_file:
        load_path = join(get_save_path(subject, date), input_file)
    else:
        load_path = join(get_raw_path(subject, date, split_recording), 
                         input_file)
    org_name_output = output_file

    raw = None
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        output_file = hp_lp_string(h_freq, l_freq) + org_name_output
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            if raw is None:
                if not split_recording:
                    raw = mne.io.read_raw_fif(load_path, preload=True)
                else:
                    raw = handle_split_recording(subject, date, input_file)
                copy = raw.copy()
            else: raw = copy
            raw = hp_lp_filter(h_freq, l_freq, raw, n_jobs)
            raw.save(save_path, overwrite=overwrite)
            

def hilbert_transform(subject, date, input_file, output_file,
                      overwrite, h_freqs, l_freqs, baseline_tmin,
                      baseline_tmax, tmin, tmax, events_file, ISI_type):
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
        
    events_file = prepend  + events_file
    
    if baseline_tmin == 'None' and baseline_tmax == 'None':
        baseline = (None, None)
        baseline_text = 'full_baseline_'
    elif baseline_tmin == 'None' and baseline_tmax == '-0.050':
        baseline = (None, -0.050)
        baseline_text = ''
    elif baseline_tmin == 'no_baseline' and baseline_tmax == 'no_baseline':
        baseline = None
        baseline_text = 'no_baseline_'
    elif baseline_tmin == 'no_baseline_Ss' and \
         baseline_tmax == 'no_baseline_Ss':
         baseline = None
         baseline_text = 'no_baseline_Ss_'
             
    else:
        raise NameError('baseline not implemented')
        
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + baseline_text + \
                      hp_lp_string(h_freq, l_freq) + \
                      time_string(tmin, tmax) + org_name_output
        
        load_path = join(get_save_path(subject, date), 'hilbert', input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        events_path = join(get_save_path(subject, date), events_file)
        if should_we_run(save_path, overwrite):
            raw = None ## conserve memory
            epochs_hilbert = None 
            raw = mne.io.read_raw_fif(load_path, preload=True)
            raw.apply_hilbert(envelope=False)
            events = mne.read_events(events_path)
            if (ISI_type != 'normal' or baseline_text != '' or \
                tmin != -0.400) and baseline_text != 'no_baseline_Ss_':
                event_id = dict(o0=18, o5=28, o15=38,
                                n1_0=48, n2_0=50, n3_0=52, n4_0=54, n5_0=56,
                                n1_5=58, n2_5=60, n3_5=62, n4_5=64, n5_5=66,
                                n1_15=68, n2_15=70, n3_15=72, n4_15=74,
                                n5_15=76
                                )
            elif baseline_text == 'no_baseline_Ss_':
                event_id = dict(s1=1, s2=3, s3=5)
            else:
                event_id = dict(s1=1, s2=3, s3=5,
                                s4_0=23, s5_0=25, s6_0=27,
                                s4_5_be=71, s5_5_be=73, s6_5_be=75,
                                s4_5_af=87, s5_5_af=89, s6_5_af=91,
                                s4_15_be=81, s5_15_be=83, s6_15_be=85,
                                s4_15_af=97, s5_15_af=99, s6_15_af=101,
                                o0=18, o5=28, o15=38,
                                n1_0=48, n2_0=50, n3_0=52, n4_0=54, n5_0=56,
                                n1_5=58, n2_5=60, n3_5=62, n4_5=64, n5_5=66,
                                n1_15=68, n2_15=70, n3_15=72, n4_15=74,
                                n5_15=76
                                )
            # tmin = -0.400
            # tmax = 1.100
            decim = 1
            reject =  dict(grad=4000e-13,
                           mag=4e-12)
    
            epochs_hilbert = mne.Epochs(raw, events, event_id, tmin, tmax,
                                        baseline, proj=False, decim=decim,
                                        reject=reject)
            epochs_hilbert.save(save_path, overwrite=overwrite)
                            
                
def evoked_hilbert(subject, date, input_file, output_file, overwrite,
                   collapse_cond, proj, h_freqs, l_freqs, tmin, tmax, ISI_type,
                   full_baseline):
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
    if proj:
        proj_string = 'proj_'
    else:
        proj_string = ''

    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend +  hp_lp_string(h_freq, l_freq) + \
            time_string(tmin, tmax) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + \
            time_string(tmin, tmax) + proj_string + org_name_output
        load_path = join(get_save_path(subject, date), 'hilbert', input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            epochs_hilbert = mne.read_epochs(load_path, proj=False)
            if not proj:
                epochs_hilbert.del_proj()
            picks = mne.pick_types(epochs_hilbert.info, meg=True, eog=True,
                                   ecg=True, emg=True, misc=True)
            if collapse_cond:
                epochs_hilbert = collapse_conditions(epochs_hilbert, ISI_type,
                                                     full_baseline, tmin)
            evokeds = []
            for event in epochs_hilbert.event_id:
                print('Averaging event: ' + event)
                evoked = epochs_hilbert[event].average(picks=picks)
                evoked.apply_proj()
                evokeds.append(evoked)
            for evoked in evokeds:
                evoked._data = np.abs(evoked.data)
                
            mne.write_evokeds(save_path, evokeds)
        epochs_hilbert = None # memory control
        
def z_transform_wilcoxon_contrasts(subject, date, input_file,
                                   output_file, overwrite, collapse_cond, proj,
                                   h_freqs, l_freqs, tmin, tmax, contrasts,
                                   ISI_type, baseline, n_projs):
    
    org_name_input = input_file
    org_name_output = output_file
    
    prepend = ISI_type + '_' + baseline + '_'
    
    if proj:
        proj_string = 'proj_'
    else:
        proj_string = ''
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + \
            time_string(tmin, tmax) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + \
            time_string(tmin, tmax) + proj_string + org_name_output
        load_path = join(get_save_path(subject, date), 'hilbert',
                           input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            epochs_hilbert = mne.read_epochs(load_path, proj=False,
                                             preload=False)
            if proj:
                epochs_hilbert.info['projs'] = \
                    epochs_hilbert.info['projs'][:n_projs]
                print('Removed first ' + str(n_projs) + ' projections')
                epochs_hilbert.apply_proj()
            else:
                epochs_hilbert.del_proj()
            if collapse_cond:
                epochs_hilbert = collapse_conditions(epochs_hilbert, ISI_type,
                                                     baseline, tmin)
            z_transformed = []
            for contrast in contrasts:
                contrast_name = contrast[0] + '_vs_' + contrast[1]
                print('Running contrast: ' + contrast_name)
                contrast_1 = epochs_hilbert[contrast[0]]
                contrast_2 = epochs_hilbert[contrast[1]]
                contrast_1.load_data()
                contrast_2.load_data()
                contrast_1.pick_types(meg=True)
                contrast_2.pick_types(meg=True)
                mne.epochs.equalize_epoch_counts([contrast_1, contrast_2])
                data_1 = contrast_1.get_data()
                data_2 = contrast_2.get_data()

                n_channels = data_1.shape[1]
                n_samples = data_1.shape[2]
                zs = np.zeros((n_channels, n_samples))
            
                for channel_index in range(n_channels):
                    for sample_index in range(n_samples):
                        this_data_1 = data_1[:, channel_index, sample_index]
                        this_data_2 = data_2[:, channel_index, sample_index]
                        z = wilcoxon(abs(this_data_1), abs(this_data_2))
                        zs[channel_index, sample_index] = z
                z_transform = mne.EvokedArray(zs, info=contrast_1.info,
                                              tmin=contrast_1.tmin,
                                              comment=contrast_name,
                                              nave=this_data_1.shape[0])
                z_transformed.append(z_transform)
                
            mne.write_evokeds(save_path, z_transformed)
        epochs_hilbert = None ## memory control        

def grand_average_evoked_hilbert(subjects, dates, input_file, output_file,
                                 overwrite, h_freqs, l_freqs, tmin, tmax,
                                 ISI_type, baseline):
    org_name_input = input_file
    org_name_output = output_file
    
    prepend_input = ISI_type + '_' + baseline + '_'
    prepend_output = ISI_type + '_' + baseline + '_'
  
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend_input + hp_lp_string(h_freq, l_freq) + \
                    time_string(tmin, tmax) + org_name_input
        output_file = prepend_output + hp_lp_string(h_freq, l_freq) + \
                    time_string(tmin, tmax) + org_name_output
        save_path = join(get_save_path('grand_averages', ''),
                         'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            evokeds_dict = dict()
            for subject_index, subject in enumerate(subjects):
                date = dates[subject_index]
                load_path = join(get_save_path(subject, date), 'hilbert',
                                 input_file)
                evokeds = mne.read_evokeds(load_path)
                for evoked in evokeds:
                    event = evoked.comment
                    print('Loading event: ' + event + ' for subject: ' + \
                          subject)
                    if subject_index == 0:
                        evokeds_dict[event] = []
                    evokeds_dict[event].append(evoked)
                evokeds = None ## free up memory

        ga_evokeds = []
        for event in evokeds_dict:
            temp = evokeds_dict[event]
            ga = mne.grand_average(temp)
            ga.comment = event
            ga_evokeds.append(ga)

        mne.evoked.write_evokeds(save_path, ga_evokeds)
        
def statistics_evoked_hilbert(subjects, dates, input_file, output_file,
                              overwrite, h_freqs, l_freqs,
                              channel_type, contrasts, p_threshold,
                              tmin, tmax, stat_tmin, stat_tmax, ISI_type,
                              baseline, n_jobs):
    org_name_input = input_file
    org_name_output = output_file
    
    prepend = ISI_type + '_' + baseline + '_'

    
    prepend_input = prepend
    prepend_output = 't-test_p_' + str(p_threshold) + '_' + prepend + \
            channel_type + '_'
            
    stat_time_string = time_string(stat_tmin, stat_tmax)
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend_input + hp_lp_string(h_freq, l_freq) + \
                    time_string(tmin, tmax) + org_name_input
        evokeds_dict = dict()
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            load_path = join(get_save_path(subject, date), 'hilbert',
                             input_file)
            evokeds = mne.read_evokeds(load_path)
            for contrast in contrasts:
                contrast_name = contrast[0] + '_vs_' + contrast[1]
                for evoked in evokeds:
                    if evoked.comment == contrast_name:
                        print('Loading contrast: ' + contrast_name + \
                              ' for subject: ' + subject) 
                        if subject_index == 0:
                            evokeds_dict[contrast_name] = []
                        evoked.crop(stat_tmin, stat_tmax)
                        evoked.pick_types(channel_type)
                        evokeds_dict[contrast_name].append(evoked)
            evokeds = None ## free up memory
        for contrast in contrasts:
            contrast_name = contrast[0] + '_vs_' + contrast[1]
            contrast_string = contrast[0] + '_' + contrast[1] + \
                '_contrast_'
            output_file = prepend_output + contrast_string + \
                hp_lp_string(h_freq, l_freq) + 'stat_' + \
                    stat_time_string + 'range_' + time_string(tmin, tmax) + \
                        org_name_output
            save_path = join(get_save_path('grand_averages', ''),
                             'statistics', 'hilbert', output_file)
            if should_we_run(save_path, overwrite):
                print('Clustering for contrast: ' + contrast_name)
                evokeds = evokeds_dict[contrast_name]
                n_subjects = len(evokeds)
                n_channels, n_times = evokeds[0].data.shape
                stat_array = np.zeros((n_subjects, n_times, n_channels))
                for evoked_index, evoked in enumerate(evokeds):
                    stat_array[evoked_index, :, :] = evoked.data.T
                
                ## running the stats
                threshold = -stats.distributions.t.ppf(p_threshold / 2.0,
                                                             n_subjects)  
                t_obs, clusters, cluster_p_values, H0 = clu = \
                    mne.stats.spatio_temporal_cluster_1samp_test(
                        stat_array, threshold, n_permutations=1024,
                        n_jobs=n_jobs, seed=7)
                cluster_dict = dict()
                names = ['t_obs', 'clusters', 'cluster_p_values', 'H0']
                for name_index, name in enumerate(names):
                    cluster_dict[name] = clu[name_index]
                np.save(save_path, cluster_dict)
                
            
                
def beamformer_contrast_hilbert(subject, date, input_file_1, input_file_2,
                                output_file, overwrite,
                                collapse_cond, whiten, proj, h_freqs, l_freqs,
                                tmin, tmax,
                                spacing, reg, weight_norm, channel_type,
                                ISI_type, contrasts, baseline_type, n_projs,
                                save_single_stcs):

    
    org_name_input_1 = input_file_1
    org_name_input_2 = input_file_2                 
    org_name_output = output_file              
    bem_folder = join(subjects_dir, subject, 'bem')
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file_1 = hp_lp_string(h_freq, l_freq) + org_name_input_1
        input_file_2 = hp_lp_string(h_freq, l_freq) + \
                       time_string(tmin, tmax) + org_name_input_2
        # prepend = ISI_and_baseline_string(ISI_type, baseline_type)
        prepend = ISI_type + '_' + baseline_type + '_' #FIXME

        input_file_2 = prepend + input_file_2
        
        load_path_1 = join(get_save_path(subject, date), 'hilbert',
                           input_file_1)
        load_path_2 = join(get_save_path(subject, date), 'hilbert',
                           input_file_2)
    
        has_been_loaded = False
        channels_deleted = False
        for contrast_index, contrast in enumerate(contrasts):
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                      contrast[1])
            
                
            paths = \
            get_load_and_save_paths_for_common_filter_hilbert(subject, date,
                                                      spacing, contrast,
                                                      h_freq, l_freq,
                                                      org_name_input_1,
                                                      org_name_output,
                                                      reg, weight_norm,
                                                      channel_type,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=None,
                                                      tmin=tmin, tmax=tmax,
                                                      ISI_type=ISI_type,
                                                      proj=proj,
                                                      baseline=baseline_type,
                                          save_single_stcs=save_single_stcs)  
                
            save_paths = paths[1]                

            if should_we_run(save_paths[0], overwrite) or \
                should_we_run(save_paths[1], overwrite) or \
                should_we_run(save_paths[2], overwrite) or \
                should_we_run(save_paths[3], overwrite) or \
                should_we_run(save_paths[4], overwrite) or \
                should_we_run(save_paths[5], overwrite):
                    
                if not has_been_loaded:
                    raw = mne.io.read_raw_fif(load_path_1, preload=False)
                    epochs_hilbert = mne.read_epochs(load_path_2,
                                                     preload=False,
                                                     proj=False)
                    if channel_type != 'both':
                        if whiten:
                            raise RuntimeError('Only whiten if you use ' + \
                                               'channel_type="both"!')     
                        picks = mne.pick_types(epochs_hilbert.info,
                                               meg=channel_type)
                    else:
                        picks = mne.pick_types(epochs_hilbert.info,
                                               meg=True)
                                         
                    if collapse_cond:
                        epochs_hilbert = collapse_conditions(epochs_hilbert,
                                                             ISI_type,
                                                             baseline_type,
                                                             tmin=tmin)
                    
                events = epochs_hilbert.events
                baseline = epochs_hilbert.baseline
                reject = epochs_hilbert.reject
                if channel_type == 'mag' and not channels_deleted:
                    del reject['grad']
                    channels_deleted = True
                if channel_type == 'grad' and not channels_deleted:
                    del reject['mag']
                    channels_deleted = True
                event_ids = epochs_hilbert.event_id
                
                new_event_id = dict()
                for event in event_ids:
                    if event in contrast:
                        new_event_id[event] = event_ids[event]
                
                epochs_cov = mne.Epochs(raw, events, new_event_id, tmin, tmax,
                                        baseline, proj=False, reject=reject,
                                        preload=True, picks=picks)
                
                if proj:
                    epochs_hilbert.info['projs'] = \
                        epochs_hilbert.info['projs'][:n_projs]
                    epochs_cov.info['projs'] = \
                        epochs_cov.info['projs'][:n_projs]
                    print('Removed first ' + str(n_projs) + ' projections')
                    epochs_hilbert.apply_proj()
                    epochs_cov.apply_proj()
                else:
                    epochs_hilbert.del_proj()
                    epochs_cov.del_proj()

                
                if not has_been_loaded:
                    ico_string = '5120'
                    trans = join(bem_folder, subject + '-trans.fif')
                    src = join(bem_folder,
                               'volume-' + str(spacing) + 'mm-src.fif')
                    bem = join(bem_folder,
                               subject + '-' + ico_string + '-bem-sol.fif')
                    
                    fwd = mne.make_forward_solution(epochs_cov.info, trans,
                                                    src, bem)
                    has_been_loaded = True
                
                if 'sss' in input_file_1:
                    rank = 'info'
                else:
                    rank = None
                if baseline_type != 'prestim':
                    data_cov = mne.compute_covariance(epochs_cov, tmin=None,
                                                      tmax=None, rank=rank)
                else:
                    data_cov = mne.compute_covariance(epochs_cov, tmin=0,
                                                      tmax=None, rank=rank)
                if whiten: ## probably throws an error
                    noise_cov = mne.compute_covariance(epochs_cov, tmin=None,
                                                       tmax=baseline[1],
                                                       rank=rank)
                else:
                    noise_cov = None
                        
                if weight_norm == 'None':
                    weight_norm = None
                filters = mne.beamformer.make_lcmv(epochs_cov.info, fwd,
                                                   data_cov=data_cov,
                                                   noise_cov=noise_cov,
                                                   pick_ori='max-power',
                                                   weight_norm=weight_norm,
                                                   reg=reg)
                epochs_cov = None
                stcs_dict = dict()
                itcs_dict = dict()
                for event_index, event in enumerate(contrast):
                    these_epochs = epochs_hilbert[event].load_data()
                    these_epochs.pick(picks)
                    stcs = \
                mne.beamformer.apply_lcmv_epochs(these_epochs,
                                                 filters, max_ori_out='signed')
                
                    if save_single_stcs:
                        print('Saving single stcs')
                        for stc_index, stc in enumerate(stcs):
                            if event_index == 0:
                                this_save_path = save_paths[0]
                            if event_index == 1:
                                this_save_path = save_paths[1]
                            replace_string = 'ms_tactile_jitter'
                            full_path = this_save_path.replace(replace_string,
                                replace_string + '_trial_' + str(stc_index))
                            stc.save(full_path, ftype='h5')
                        stcs = None
                
                    
                    
                    if not save_single_stcs:
                        print('Running ITC')
                        itc = compute_source_itc(stcs)
                        itc_stc = stcs[0].copy()
                        itc_stc._data = itc
                        itcs_dict[event] = itc_stc
                        
                        print('Getting absolute values')
                        for stc in stcs:
                            stc._data = np.array(np.abs(stc.data),
                                                 dtype='float64')
                        print('Averaging response')    
                        stc_mean = stcs[0].copy()
                        mean_data = np.mean([stc.data for stc in stcs], axis=0)
                        stc_mean._data = mean_data
                        stcs_dict[event] = stc_mean
                    
                if not save_single_stcs:
                ## save_single conditions
                    for condition_index, condition in enumerate(contrast):
                        stc = stcs_dict[condition]
                        itc = itcs_dict[condition]
                        if condition_index == 0:
                            stc.save(save_paths[0])
                            itc.save(save_paths[2])
                        if condition_index == 1:
                            stc.save(save_paths[1])
                            itc.save(save_paths[3])
                        if len(contrast) == 3:
                            if condition_index == 2:
                                stc.save(save_paths[4])
                                itc.save(save_paths[5])
                            
                    if len(contrast) < 3:
                            
                        ## calculate ratio and save
                        for condition_index, condition in enumerate(contrast):
                            if condition_index == 0:
                                ratio_stc = stcs_dict[condition].copy()
                                diff_itc = itcs_dict[condition].copy()
                            if condition_index == 1:
                                ratio_stc._data = (ratio_stc.data - \
                                           stcs_dict[condition].data) / \
                                    (ratio_stc.data + \
                                     stcs_dict[condition].data)
                                diff_itc._data = diff_itc.data - \
                                    itcs_dict[condition].data
                                                      
                        ratio_stc.save(save_paths[4])
                        diff_itc.save(save_paths[5])
   
                ## memory control
                stcs = None
                stcs_dict = None
                itcs_dict = None
                ratio_stc = None
                diff_itc = None
def morph_beamformer_contrast_hilbert(subject, date, input_file,
                                      input_file_morph, output_file, overwrite,
                                      proj, h_freqs, l_freqs, tmin, tmax,
                                      spacing,
                                      reg, weight_norm, channel_type, ISI_type,
                                      contrasts, baseline, save_single_stcs):
    
    org_name_input = input_file    
    org_name_output = output_file
    input_file_morph = str(spacing) + '_mm_' + input_file_morph
    load_path_morph = join(get_save_path(subject, date), input_file_morph)
    morph = mne.read_source_morph(load_path_morph)
    
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])
            paths = \
    get_load_and_save_paths_for_common_filter_hilbert(subject, date, spacing,
                                                      contrast, h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output, reg,
                                                      weight_norm,
                                                      channel_type,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=None,
                                                      tmin=tmin, tmax=tmax,
                                                      ISI_type=ISI_type,
                                                      proj=proj,
                                                      baseline=baseline,
                                          save_single_stcs=save_single_stcs)
            load_paths = paths[0]
            save_paths = paths[1]
            load_path_index = 0
            
            for (load_path, save_path) in zip(load_paths, save_paths):
                if load_path_index == 2:
                    break
                if should_we_run(save_path, overwrite):
                    if save_single_stcs:
                        poststring = '-vl-stc.h5'
                        print(load_path)
                        prestring = load_path[:load_path.index(poststring)]
                        folder_name = dirname(prestring)
                        files = listdir(folder_name)
                        file_counter = 0
                        for File in files:
                            condition_string = \
                                prestring[(len(folder_name) + 1):]
                            if File.find(condition_string) != -1:
                                file_counter += 1
                        for file_index in range(file_counter):
                            trial_string = '_trial_' + str(file_index)
                            this_load_path = prestring + trial_string + \
                                            poststring
                            this_save_path = prestring + '_morph' + \
                                             trial_string + poststring
                            ## FIXME: add overwrite protection
                            stc = mne.read_source_estimate(this_load_path)
                            print('Applying morph no.: ' + str(file_index))
                            stc_morph = morph.apply(stc)
                            print('Saving morph no.: ' + str(file_index))
                            stc_morph.save(this_save_path, ftype='h5')
                
                    else: ## evoked stc
                        stc = mne.read_source_estimate(load_path)
                        stc_morph = morph.apply(stc)
                        stc_morph.save(save_path)     
                load_path_index += 1
            
def grand_average_beamformer_contrast_hilbert(subjects, dates, input_file,
                                              output_file, overwrite, proj,
                                              h_freqs,
                                              l_freqs, tmin, tmax, spacing,
                                              reg, weight_norm, channel_type,
                                              ISI_type, contrasts,
                                              baseline, save_std):
    org_name_input = input_file
    org_name_output = output_file
    
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])

            for subject_index, subject in enumerate(subjects):
                date = dates[subject_index]

                paths = \
                get_load_and_save_paths_for_common_filter_hilbert(subject,
                                                                  date, 
                                                                  spacing,
                                                    contrast, h_freq, l_freq,
                                                    org_name_input,
                                                    org_name_output, reg,
                                                    weight_norm,
                                                    channel_type,
                                                    grand_average=False,
                                                    statistics=False,
                                                    p_threshold=None,
                                                    tmin=tmin, tmax=tmax,
                                                    ISI_type=ISI_type,
                                                    proj=proj,
                                                    baseline=baseline)
                load_paths = paths[0]
                
                paths = \
                get_load_and_save_paths_for_common_filter_hilbert(subject,
                                                                  date, 
                                                                  spacing,
                                                    contrast, h_freq, l_freq,
                                                    org_name_input,
                                                    org_name_output, reg,
                                                    weight_norm,
                                                    channel_type,
                                                    grand_average=True,
                                                    statistics=False,
                                                    p_threshold=None,
                                                    tmin=tmin, tmax=tmax,
                                                    ISI_type=ISI_type,
                                                    proj=proj,
                                                    baseline=baseline)
                save_paths = paths[1]
                n_paths = len(save_paths)
                if subject_index == 0:
                    stcs_list = [list([]) for _ in range(n_paths)]
                    running_gas = [True] * n_paths                          
                for path_index, (load_path, save_path) in \
                    enumerate(zip(load_paths, save_paths)):
                    if not should_we_run(save_path, overwrite,
                                         suppress_print=True):
                        running_gas[path_index] = False
                    else:
                        print('Loading: ' + load_path)
                        stc = mne.read_source_estimate(load_path)
                        stcs_list[path_index].append(stc)
                    
            ## run ga and save          
            for path_index in range(n_paths):
                stcs = stcs_list[path_index]
                n_subjects = len(stcs)
                n_sources = stcs[0].data.shape[0]
                n_times = stcs[0].data.shape[1]
                if running_gas[path_index]:
                    save_path = save_paths[path_index]
                    print('Saving grand average: ' + save_path)
                    grand_average = sum(stcs) / n_subjects
                    grand_average.save(save_path)
                if save_std:
                    stcs_data = np.zeros((n_subjects, n_sources, n_times))
                    for subject_index in range(n_subjects):
                        stcs_data[subject_index, :, :] = \
                            stcs[subject_index].data
                        
                    std_data = np.std(stcs_data, axis=0)
                    std_save_path = save_path[:-7] + '_std' + \
                        save_path[-7:]
                    std_stc = stcs[0].copy()
                    std_stc._data = std_data
                    print('Saving std to: ' + std_save_path)
                    std_stc.save(std_save_path)
                else:
                    print('Not overwriting: ' + save_path)
                    
def statistics_beamformer_contrast_hilbert(subjects, dates, input_file,
                                           output_file, overwrite, proj,
                                           h_freqs,
                                           l_freqs, spacing, reg, weight_norm,
                                           channel_type, contrasts,
                                           p_threshold, tmin, tmax,
                                           stat_tmin, stat_tmax,
                                           Type, ISI_type, stat_function,
                                           baseline, n_jobs):
    org_name_input = input_file
    org_name_output = output_file
    
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'volume-' + str(spacing) + \
                                                     'mm-src.fif'))

    
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):

        stcs_dict = dict()
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            for contrast in contrasts:
                print('Loading contrast: ' + contrast[0] + ' versus ' + \
                      contrast[1] + ' for subject: ' + subject)
                contrast_name = contrast[0] + '_vs_' + contrast[1]
                
                paths = \
                get_load_and_save_paths_for_common_filter_hilbert(subject,
                                                                  date, 
                                                                  spacing,
                                                    contrast, h_freq, l_freq,
                                                    org_name_input,
                                                    org_name_output, reg,
                                                    weight_norm,
                                                    channel_type,
                                                    grand_average=False,
                                                    statistics=False,
                                                    p_threshold=p_threshold,
                                                    tmin=tmin,
                                                    tmax=tmax,
                                                    ISI_type=ISI_type,
                                                    proj=proj,
                                                    baseline=baseline)
                if Type == 'stc':
                    load_path = paths[0][4]
                elif Type == 'itc':
                    load_path = paths[0][5]
                    
                stc = mne.read_source_estimate(load_path)
                if subject_index == 0:
                    stcs_dict[contrast_name] = []
                
                stcs_dict[contrast_name].append(stc)
    
        stat_dict = dict()
        for contrast in contrasts:
            contrast_name = contrast[0] + '_vs_' + contrast[1]
            stcs = stcs_dict[contrast_name]
            stat_array = None
            
            paths = \
            get_load_and_save_paths_for_common_filter_hilbert(subject,
                                                              date, 
                                                              spacing,
                                                contrast, h_freq, l_freq,
                                                org_name_input,
                                                org_name_output, reg,
                                                weight_norm,
                                                channel_type,
                                                grand_average=False,
                                                statistics=True,
                                                p_threshold=p_threshold,
                                                tmin=tmin, tmax=tmax,
                                                ISI_type=ISI_type,
                                                stat_function=stat_function,
                                                proj=proj,
                                                baseline=baseline,
                                                stat_tmin=stat_tmin,
                                                stat_tmax=stat_tmax)
            if Type == 'stc':
                save_path = paths[1][4]
            elif Type == 'itc':
                save_path = paths[1][5]
                
            if should_we_run(save_path, overwrite):
                for stc in stcs:
                    temp = stc.copy()
                    temp._data = temp.data.astype('float32')
                    temp.crop(stat_tmin, stat_tmax)

                    if stat_array is None:
                        stat_array = temp.data.T
                        stat_array = np.expand_dims(stat_array, 0)
                    else:
                        temp = temp.data.T
                        temp = np.expand_dims(temp, 0)
                        stat_array = np.concatenate((stat_array, temp), 0)
                print('Clustering for contrast: ' + contrast_name)
                n_times = stat_array.shape[1]
                connectivity = mne.spatio_temporal_src_connectivity(src,
                                                                    n_times)
                n_subjects = len(subjects)
                if stat_function == 't-test':
                    threshold = -stats.distributions.t.ppf(p_threshold / 2.0,
                                                             n_subjects - 1)
                    stat_fun = mne.stats.ttest_1samp_no_p
                    buffer_size = 1000
                elif stat_function == 'wilcoxon':
                    threshold = \
                        -stats.distributions.norm.ppf(p_threshold / 2.0)
                    ## only exists in my git branch
                    stat_fun = mne.stats.wilcoxon
                    buffer_size = None
                        
                t_obs, clusters, cluster_p_values, H0 = clu = \
                    mne.stats.spatio_temporal_cluster_1samp_test(stat_array,
                        connectivity=connectivity, n_jobs=n_jobs,
                        threshold=threshold, verbose=True, max_step=1,
                        seed=7, n_permutations=1024, stat_fun=stat_fun,
                        buffer_size=buffer_size)
                stat_dict[contrast_name] = clu
                
                print('Saving contrast: ' + contrast_name)
    
                clu = stat_dict[contrast_name]
                cluster_dict = dict()
                names = ['t_obs', 'clusters', 'cluster_p_values', 'H0']
                for name_index, name in enumerate(names):
                    cluster_dict[name] = clu[name_index]
                np.save(save_path, cluster_dict)
                
def save_nifti(subject, date, input_file, output_file, overwrite, h_freqs,
               l_freqs, tmin, tmax, spacing, reg, weight_norm, channel_type,
               ISI_type, contrasts, baseline):
    
    org_name_input = input_file
    org_name_output = output_file
    
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'volume-' + str(spacing) + \
                                                     'mm-src.fif'))
            
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])
            paths = \
    get_load_and_save_paths_for_common_filter_hilbert(subject, date, spacing,
                                                      contrast, h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output, reg,
                                                      weight_norm,
                                                      channel_type,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=None,
                                                      tmin=tmin, tmax=tmax,
                                                      ISI_type=ISI_type,
                                                      baseline=baseline,
                                                      proj=False)
            load_paths = paths[0]
            save_paths = paths[1]
          
            for (load_path, save_path) in zip(load_paths, save_paths):
            
                if should_we_run(save_path, overwrite):
                    stc = mne.read_source_estimate(load_path)
                    print('Saving:' + save_path)
                    stc.save_as_volume(save_path, src)

def save_nifti_stat(input_file_1, input_file_2, output_file, overwrite, h_freq,
                    l_freq, spacing, reg, weight_norm, channel_type, contrasts,
                    p_threshold, tmin, tmax, stat_tmin, stat_tmax, Type, 
                    ISI_type, stat_function):
    
    org_name_input_1 = input_file_1
    org_name_input_2 = input_file_2
    org_name_output = output_file
    
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'volume-' + str(spacing) + \
                                                     'mm-src.fif'))

    for contrast in contrasts:
        print('Loading contrast: ' + contrast[0] + ' versus ' + contrast[1])
        
        ## load cluster
        paths = \
        get_load_and_save_paths_for_common_filter_hilbert('',
                                                          '', 
                                                          spacing,
                                            contrast, h_freq, l_freq,
                                            org_name_input_1,
                                            org_name_output, reg,
                                            weight_norm,
                                            channel_type,
                                            grand_average=False,
                                            statistics=True,
                                            p_threshold=p_threshold,
                                            tmin=stat_tmin,
                                            tmax=stat_tmax,
                                            ISI_type=ISI_type,
                                            stat_function=stat_function)
        if Type == 'stc':
            load_path = paths[0][4]
            save_path = paths[1][4]
        elif Type == 'itc':
            load_path = paths[0][5]
            save_path = paths[1][5]
            
        if should_we_run(save_path, overwrite):
            clu = np.load(load_path, allow_pickle=True).item()

            ## load stc
            paths = \
            get_load_and_save_paths_for_common_filter_hilbert(
                                                '',
                                                '', 
                                                spacing,
                                                contrast, h_freq, l_freq,
                                                org_name_input_2,
                                                org_name_output, reg,
                                                weight_norm,
                                                channel_type,
                                                grand_average=True,
                                                statistics=False,
                                                p_threshold=p_threshold,
                                                tmin=tmin,
                                                tmax=tmax,
                                                ISI_type=ISI_type)
            if Type == 'stc':
                load_path = paths[0][4]
            elif Type == 'itc':
                load_path = paths[0][5]
                
            stc = mne.read_source_estimate(load_path)
                
            # t_obs = clu['t_obs']
            # H0 = clu['H0']
            clusters = clu['clusters']
            cluster_p_values = clu['cluster_p_values']
            cluster_stc = stc.copy()
            cluster_stc.crop(stat_tmin, stat_tmax)
            org_data = cluster_stc.data
            cluster_stc._data = np.zeros(cluster_stc.data.shape)
            cluster_index = np.argmin(cluster_p_values) ## FIXME
                                                    ##  only plotting min
            assert (cluster_p_values[cluster_index] < p_threshold), \
                    'There are no p_values smaller than p_threshold (' + \
                        str(p_threshold) + ')'
                        
            cluster = clusters[cluster_index]
            cluster_stc._data[cluster[1], cluster[0]] = org_data[cluster[1],
                                                                 cluster[0]]
            print('Saving:' + save_path)
            cluster_stc.save_as_volume(save_path, src)
            

def envelope_correlations():
    pass
                      
            
def read_and_resave_epochs(subject, date, input_file, output_file, overwrite,
                           h_freqs, l_freqs, ISI_type, tmin, tmax,
                           full_baseline):
    org_name_input = input_file
    org_name_output = output_file
    if full_baseline:
        baseline_text = 'full_baseline_'
    else:
        baseline_text = ''
    prepend = ISI_type + '_' + baseline_text

    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend +  hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + \
                      time_string(tmin, tmax) + org_name_output
        load_path = join(get_save_path(subject, date), 'hilbert', input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            epochs_hilbert = mne.read_epochs(load_path, proj=False)
            epochs_hilbert.save(save_path, overwrite=overwrite)
            epochs_hilbert = None
             
def plot_alignment(subject, date, split_recording, input_file, output_file,
                   overwrite):
    mne.viz.set_3d_backend('pyvista')
    bem_folder = join(subjects_dir, subject, 'bem')
    trans = mne.read_trans(join(bem_folder, subject + '-trans.fif'))
    raw_path = get_raw_path(subject, date, '')
    load_path = join(raw_path, input_file)
    save_path = join(get_fig_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        if not split_recording:
            info = mne.io.read_info(load_path)
        else:
            info = handle_split_recording_info(subject, date, input_file)
        fig = mne.viz.plot_alignment(info, trans, subject, subjects_dir,
                                     dig=True, meg='sensors', show_axes=True,
                                     coord_frame='meg', surfaces='head-dense')
        mne.viz.set_3d_view(fig, 0, 90)
        fig.plotter.save_graphic(save_path)
        fig.plotter.close()
                
def plot_evoked_hilbert(subject, date, input_file, output_file, overwrite,
                         h_freqs, l_freqs):
    org_name_input = input_file
    org_name_output = output_file
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = hp_lp_string(h_freq, l_freq) + org_name_input
        load_path = join(get_save_path(subject, date), 'hilbert', input_file)
        evokeds = mne.read_evokeds(load_path, proj=False)
        
        for evoked in evokeds:
            output_file = evoked.comment + '_' + \
                hp_lp_string(h_freq, l_freq) + org_name_output       
            save_path = join(get_fig_path(subject, date), 'hilbert', 'evokeds',
                             output_file)
            if should_we_run(save_path, overwrite):
                print('Plotting: ' + evoked.comment)
                fig = evoked.plot()
                fig.savefig(save_path)
                
def plot_projs_evoked_hilbert(subject, date, input_file, output_file, 
                              overwrite, h_freqs, l_freqs, ISI_type):
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'

    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_input
        load_path = join(get_save_path(subject, date), 'hilbert', input_file)
        evokeds = mne.read_evokeds(load_path, proj=True)
        
        for evoked in evokeds[0:1]:
            output_file = prepend + \
                hp_lp_string(h_freq, l_freq) + org_name_output       
            save_path = join(get_fig_path(subject, date), 'hilbert', 'evokeds',
                             output_file)
            if should_we_run(save_path, overwrite):
                print('Plotting: ' + evoked.comment)
                print(save_path)
                fig = evoked.plot_projs_topomap()
                fig.savefig(save_path)                
                
def plot_all_together_evoked_hilbert(subjects, dates, input_file, output_file,
                                     overwrite, h_freqs, l_freqs,
                                     channel_type, ISI_type):
    org_name_input = input_file
    org_name_output = output_file
    n_subjects = len(subjects)
    prepend = ISI_type + '_'
    
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        plt.close('all')
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            input_file = prepend + hp_lp_string(h_freq, l_freq) + \
                         org_name_input
            load_path = join(get_save_path(subject, date), 'hilbert', 
                             input_file)
            evokeds = mne.read_evokeds(load_path, proj=False)
            n_evokeds = len(evokeds)
            if subject_index == 0:
                figs = [None] * n_evokeds
                axes_list = [None] * n_evokeds
            for evoked_index, evoked in enumerate(evokeds):
                output_file = prepend + channel_type + '_' + \
                              evoked.comment + '_' + \
                    hp_lp_string(h_freq, l_freq) + org_name_output       
                save_path = join(get_fig_path('all_together', ''), 'hilbert',
                                 'evokeds', output_file)
                if should_we_run(save_path, overwrite):
                    print('Adding subject: ' + subject + ' to plot')
                    if subject_index == 0:
                        n_rows = 5
                        n_cols = 6
                        fig, axes = plt.subplots(n_rows, n_cols)
                        fig.set_figheight(14.4)
                        fig.set_figwidth(25.6)
                        figs[evoked_index] = fig
                        axes_list[evoked_index] = axes
                    fig = figs[evoked_index]
                    axes = axes_list[evoked_index]
                    row_number = (subject_index // n_cols)
                    col_number = (subject_index % n_cols)
                    evoked.plot(picks=[channel_type], axes=axes[row_number,
                                                                col_number],
                                titles=dict(grad='', mag=''))
                    if (subject_index + 1) == n_subjects:
                        fig.savefig(save_path)
                        
def plot_beamformer_contrast_hilbert(subject, date, input_file, output_file,
                                     overwrite, h_freqs, l_freqs, tmin, tmax,
                                     spacing, reg, weight_norm, channel_type,
                                     ISI_type, contrasts,
                                     initial_time,
                                     initial_pos):
    org_name_input = input_file
    org_name_output = output_file
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'volume-' + str(spacing) + \
                                                     'mm-src.fif'))
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])
            paths = \
    get_load_and_save_paths_for_common_filter_hilbert(subject, date, spacing,
                                                      contrast, h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output, reg,
                                                      weight_norm,
                                                      channel_type,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=None,
                                                      tmin=tmin, tmax=tmax,
                                                      ISI_type=ISI_type,
                                                      figures=True)
            load_paths = paths[0]
            save_paths = paths[1]
            
            for (load_path, save_path) in zip(load_paths, save_paths):
            
                if should_we_run(save_path, overwrite):
                    stc = mne.read_source_estimate(load_path)
                    fig = stc.plot(src, src[0]['subject_his_id'],
                                   initial_time=initial_time,
                                   initial_pos=initial_pos)
                    
                    pos_string = str(int(initial_pos[0] * 1e3 )) + '_' + \
                                 str(int(initial_pos[1] * 1e3 )) + '_' + \
                                 str(int(initial_pos[2] * 1e3 )) + '_xyz_mm_'
                    save_path = str(initial_time) + '_ms_' + \
                                pos_string + save_path
                                
                    fig.savefig(save_path)
                
                
def plot_beamformer_time_courses_hilbert(subject, date, input_file,
                                         output_file, overwrite, proj,
                                         h_freqs, l_freqs, tmin, tmax, spacing,
                                         reg, weight_norm, channel_type,
                                         ISI_type, contrasts, baseline,
                                         roi):
    
    rois = dict(
                Left_Cerebellum_6  =  4562,
                Right_Cerebellum_6 =  2705,
                left_thalamus    =  7186,
                right_thalamus   =  7188,
                Left_Putamen     =  7207,
                Right_Putamen    =  7213,
                left_insula      =  7850,
                right_insula     =  7859,
                left_SI          = 12083,
                right_SI         = 12067,
                Left_Inferior_Parietal_Cortex = 10148,
                right_IPC        = 10160,
                left_MCC         =  9580,
                right_MCC        =  9556,
                )
    
    vertex = rois[roi]
    
    org_name_input = input_file
    org_name_output = roi + '_' + output_file

    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for contrast in contrasts:
            print('Running contrast: ' + str(contrasts))
            paths = \
    get_load_and_save_paths_for_common_filter_hilbert(subject, date, spacing,
                                                      contrast, h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output, reg,
                                                      weight_norm,
                                                      channel_type,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=None,
                                                      tmin=tmin, tmax=tmax,
                                                      ISI_type=ISI_type,
                                                      proj=proj,
                                                      baseline=baseline,
                                                      figures=True)
            load_paths = paths[0]
            save_paths = paths[1]
            plt.close('all')
            fig_stc = plt.figure() ## no. 1
            fig_itc = plt.figure() ## no. 2
            fig_stc.set_size_inches(9, 6)
            fig_itc.set_size_inches(9, 6)
            for (load_path, save_path) in zip(load_paths, save_paths):
            
                if should_we_run(save_path, overwrite=True):
                    stc = mne.read_source_estimate(load_path)
                    if 'stcs' in load_path:
                        plt.figure(1)
                    elif 'itcs' in load_path:
                        plt.figure(2)
                    data_index = np.where(stc.vertices == vertex)[0]
                    plt.plot(stc.times * 1e3, stc.data[data_index, :].T)
            ## stc fig
            plt.figure(1)
            plt.title('Subject: ' + subject + ' Stc: ' + \
                      roi.replace('_', ' ' ))
            plt.xlim(tmin * 1e3, tmax * 1e3)
            plt.xlabel('Time (ms)')
            plt.ylabel('Activation (arbitrary unit)')
            plt.plot((0, 0), fig_stc.axes[0].get_ylim(), 'k--')
            plt.legend(contrast)
            print('Saving: ' + save_paths[4])
            fig_stc.savefig(save_paths[4])   
            
            ## itc fig
            plt.figure(2)
            plt.title('Subject: ' + subject + ' ITC: ' + \
                      roi.replace('_', ' ' ))
            plt.xlim(tmin * 1e3, tmax * 1e3)
            plt.xlabel('Time (ms)')
            plt.ylabel('Inter-trial coherence')
            plt.plot((0, 0), fig_itc.axes[0].get_ylim(), 'k--')
            plt.legend(contrast)
            print('Saving: ' + save_paths[5])
            fig_itc.savefig(save_paths[5])
                                           
                        
def plot_whitened_data_cov(subject, date, input_file_1, input_file_2,
                               output_file, overwrite,
                               collapse_cond, whiten, h_freqs, l_freqs,
                               spacing, reg, weight_norm, channel_type,
                               contrasts):

    
    org_name_input_1 = input_file_1
    org_name_input_2 = input_file_2                 
    org_name_output = output_file              
    bem_folder = join(subjects_dir, subject, 'bem')
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file_1 = hp_lp_string(h_freq, l_freq) + org_name_input_1
        input_file_2 = hp_lp_string(h_freq, l_freq) + org_name_input_2
        load_path_1 = join(get_save_path(subject, date), 'hilbert',
                           input_file_1)
        load_path_2 = join(get_save_path(subject, date), 'hilbert',
                           input_file_2)
    
        has_been_loaded = False
        for contrast_index, contrast in enumerate(contrasts):
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                      contrast[1])
            output_file = channel_type + '_' + \
                            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[0],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + org_name_output
   
            save_path = join(get_fig_path(subject, date), 'hilbert', 'cov',
                       output_file)
            
            if should_we_run(save_path, overwrite):
            
                if not has_been_loaded:
                    raw = mne.io.read_raw_fif(load_path_1, preload=False)
                    epochs_hilbert = mne.read_epochs(load_path_2,
                                                     preload=False,
                                                     proj=False)
                    if channel_type != 'both':
                        if whiten:
                            raise RuntimeError('Only whiten if you use ' + \
                                               'channel_type="both"!')     
                        picks = mne.pick_types(epochs_hilbert.info,
                                               meg=channel_type)
                    else:
                        picks = mne.pick_types(epochs_hilbert.info,
                                               meg=True)
                                         
                    if collapse_cond:
                        epochs_hilbert = collapse_conditions(epochs_hilbert)
                    
                events = epochs_hilbert.events
                tmin = epochs_hilbert.tmin
                tmax = epochs_hilbert.tmax
                baseline = epochs_hilbert.baseline
                reject = epochs_hilbert.reject
                if channel_type == 'mag' and contrast_index == 0:
                    del reject['grad']
                if channel_type == 'grad' and contrast_index == 0:
                    del reject['mag']
                event_ids = epochs_hilbert.event_id
                
                new_event_id = dict()
                for event in event_ids:
                    if event in contrast:
                        new_event_id[event] = event_ids[event]
                
                epochs_cov = mne.Epochs(raw, events, new_event_id, tmin, tmax,
                                        baseline, proj=False, reject=reject,
                                        preload=True, picks=picks)
                epochs_cov.del_proj()
                
                if not has_been_loaded:
                    ico_string = '5120'
                    trans = join(bem_folder, subject + '-trans.fif')
                    src = join(bem_folder,
                               'volume-' + str(spacing) + 'mm-src.fif')
                    bem = join(bem_folder,
                               subject + '-' + ico_string + '-bem-sol.fif')
                    
                    fwd = mne.make_forward_solution(epochs_cov.info, trans,
                                                    src, bem)
                    has_been_loaded = True
                
                if 'sss' in input_file_1:
                    rank = 'info'
                else:
                    rank = None
                data_cov = mne.compute_covariance(epochs_cov, tmin=0,
                                                  tmax=None, rank=rank)
                if whiten:
                    noise_cov = mne.compute_covariance(epochs_cov, tmin=None,
                                                       tmax=baseline[1],
                                                       rank=rank)
                else:
                    noise_cov = None
                        
                if weight_norm == 'None':
                    weight_norm = None
    
                (cov, whitener) = \
    mne.beamformer.make_lcmv_retrieve_whitened_data_cov(epochs_cov.info,
                                                    fwd,
                                                    data_cov=data_cov,
                                                    noise_cov=noise_cov,
                                                    pick_ori='max-power',
                                                    weight_norm=weight_norm,
                                                    reg=reg)
                s = linalg.svd(cov, compute_uv=False)
                fig = plt.figure()
                plt.plot(s)
                plt.xlabel = '# Component'
                plt.ylabel('Eigenvalue')
                print('Saving: ' + save_path)
                fig.savefig(save_path)
        
        
#%% FUNCTION CALLS


if function_to_run == 'high_and_lowpass':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freq = int(argv[8])
    l_freq = int(argv[9])

    high_and_lowpass(subject, date, split_recording, input_file, output_file,
                     overwrite, h_freq, l_freq)
if function_to_run == 'find_events':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))

    find_events(subject, date, input_file, output_file, overwrite)

if function_to_run == 'epoch_data':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freq = int(argv[8])
    l_freq = int(argv[9])
    events_file = argv[10]

    epoch_data(subject, date, split_recording, input_file, output_file,
               overwrite, h_freq, l_freq, events_file)

if function_to_run == 'get_evokeds':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])
    collapse_cond = bool(int(argv[9]))

    get_evokeds(subject, date, input_file, output_file, overwrite,
                     h_freq, l_freq, collapse_cond)

if function_to_run == 'grand_average_evokeds':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freq = int(argv[5])
    l_freq = int(argv[6])

    grand_average_evokeds(subjects, dates, input_file, output_file, overwrite,
                          h_freq, l_freq)

if function_to_run == 'watershed':
    subject = argv[2]
    overwrite = bool(int(argv[3]))

    watershed(subject, overwrite)

if function_to_run == 'source_space_and_bem_model_and_bem_solution':
    subject = argv[2]

    source_space_and_bem_model_and_bem_solution(subject)

if function_to_run == 'volumetric_source_space':
    subject = argv[2]
    spacing = float(argv[3])

    volumetric_source_space(subject, spacing)

if function_to_run == 'morph_to_fsaverage':
    subject = argv[2]
    date = argv[3]
    output_file  = argv[4]
    overwrite = bool(int(argv[5]))
    ico = int(argv[6])

    morph_to_fsaverage(subject, date, output_file, overwrite, ico)

if function_to_run == 'morph_to_fsaverage_volume':
    subject = argv[2]
    date = argv[3]
    output_file  = argv[4]
    overwrite = bool(int(argv[5]))
    spacing = float(argv[6])

    morph_to_fsaverage_volume(subject, date, output_file, overwrite, spacing)

if function_to_run == 'forward_solution':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    spacing = argv[8]
    ico = int(argv[9])

    forward_solution(subject, date,  split_recording, input_file, output_file,
                     overwrite, spacing, ico)

if function_to_run == 'forward_solution_volumetric':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    spacing = int(argv[8])
    ico = int(argv[9])

    forward_solution_volumetric(subject, date,  split_recording, input_file,
                                output_file, overwrite, spacing, ico)

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
    split_recording = bool(int(argv[4]))
    input_file_1 = argv[5]
    input_file_2 = argv[6]
    input_file_3 = argv[7]
    output_file = argv[8]
    overwrite = bool(int(argv[9]))
    h_freq = int(argv[10])
    l_freq = int(argv[11])

    inverse_operator(subject, date, split_recording, input_file_1,
                     input_file_2, input_file_3, output_file, overwrite,
                     h_freq, l_freq)

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

if function_to_run == 'morph_minimum_norm_estimate':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    input_file_morph = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freq = int(argv[8])
    l_freq = int(argv[9])
    collapse_cond = bool(int(argv[10]))
    ico = int(argv[11])

    morph_minimum_norm_estimate(subject, date, input_file, input_file_morph,
                                output_file, overwrite, h_freq, l_freq,
                                collapse_cond, ico)

if function_to_run == 'grand_average_stc':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freq = int(argv[5])
    l_freq = int(argv[6])
    collapse_cond = bool(int(argv[7]))
    ico =  int(argv[8])
        
    grand_average_stc(subjects, dates, input_file, output_file, overwrite,
                      h_freq, l_freq, collapse_cond, ico)    

if function_to_run == 'statistics_stc':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freq = int(argv[5])
    l_freq = int(argv[6])
    contrasts = list_of_list_parser(argv[7])
    p_threshold = float(argv[8])
    tmin = float(argv[9])
    tmax = float(argv[10])
    ico = int(argv[11])
    n_jobs = int(argv[12])

    statistics_stc(subjects, dates, input_file, output_file, overwrite,
                   h_freq, l_freq, contrasts, p_threshold, tmin, tmax, ico,
                   n_jobs)
   
if function_to_run == 'filter_for_hilbert':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freqs = list(map(int, list_parser(argv[8])))
    l_freqs = list(map(int, list_parser(argv[9])))
    n_jobs = int(argv[10])
  
    filter_for_hilbert(subject, date, split_recording, input_file, output_file,
                       overwrite, h_freqs, l_freqs, n_jobs)    
    
if function_to_run == 'hilbert_transform':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    baseline_tmin = argv[9]
    baseline_tmax = argv[10]
    tmin = float(argv[11])
    tmax = float(argv[12])
    events_file = argv[13]
    ISI_type = argv[14]
  
    hilbert_transform(subject, date, input_file, output_file,
                      overwrite, h_freqs, l_freqs, baseline_tmin,
                      baseline_tmax, tmin, tmax, events_file, ISI_type)
    
if function_to_run == 'evoked_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))
    proj = bool(int(argv[8]))
    h_freqs = list(map(int, list_parser(argv[9])))
    l_freqs = list(map(int, list_parser(argv[10])))
    tmin = float(argv[11])
    tmax = float(argv[12])
    ISI_type = argv[13]
    full_baseline = bool(int(argv[14]))
  
    evoked_hilbert(subject, date, input_file, output_file, overwrite,
                   collapse_cond, proj, h_freqs, l_freqs, tmin, tmax, ISI_type,
                   full_baseline)
        
if function_to_run == 'z_transform_wilcoxon_contrasts':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))
    proj = bool(int(argv[8]))
    h_freqs = list(map(int, list_parser(argv[9])))
    l_freqs = list(map(int, list_parser(argv[10])))
    tmin = float(argv[11])
    tmax = float(argv[12])
    contrasts = list_of_list_parser(argv[13])
    ISI_type = argv[14]
    baseline = argv[15]
    n_projs = int(argv[16])
    
    z_transform_wilcoxon_contrasts(subject, date, input_file,
                                   output_file, overwrite, collapse_cond, proj,
                                   h_freqs, l_freqs, tmin, tmax, contrasts,
                                   ISI_type, baseline, n_projs)    
    
if function_to_run == 'grand_average_evoked_hilbert':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    tmin = float(argv[7])
    tmax = float(argv[8])
    ISI_type = argv[9]
    baseline = argv[10]

    grand_average_evoked_hilbert(subjects, dates, input_file, output_file,
                                 overwrite, h_freqs, l_freqs, tmin, tmax,
                                 ISI_type, baseline)
    
if function_to_run == 'statistics_evoked_hilbert':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    channel_type = argv[7]
    contrasts = list_of_list_parser((argv[8]))
    p_threshold = float(argv[9])
    tmin = float(argv[10])
    tmax = float(argv[11])
    stat_tmin = float(argv[12])
    stat_tmax = float(argv[13])
    ISI_type = argv[14]
    baseline = argv[15]
    n_jobs = int(argv[16])
    
    statistics_evoked_hilbert(subjects, dates, input_file, output_file,
                              overwrite, h_freqs, l_freqs, channel_type,
                              contrasts, p_threshold, tmin, tmax, stat_tmin,
                              stat_tmax, ISI_type, baseline, n_jobs)

    
if function_to_run == 'beamformer_contrast_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file_1 = argv[4]
    input_file_2 = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    collapse_cond = bool(int(argv[8]))
    whiten = bool(int(argv[9]))
    proj = bool(int(argv[10]))
    h_freqs = list(map(int, list_parser(argv[11])))
    l_freqs = list(map(int, list_parser(argv[12])))
    tmin = float(argv[13])
    tmax = float(argv[14])
    spacing = float(argv[15])
    reg = float(argv[16])
    weight_norm = argv[17]
    channel_type = argv[18]
    ISI_type = argv[19]
    contrasts = list_of_list_parser(argv[20])
    baseline_type = argv[21]
    n_projs = int(argv[22])
    save_single_stcs = bool(int(argv[23]))
    
    beamformer_contrast_hilbert(subject, date, input_file_1, input_file_2,
                                output_file, overwrite, collapse_cond, whiten,
                                proj, h_freqs, l_freqs, tmin, tmax, spacing,
                                reg,
                                weight_norm, channel_type, ISI_type, contrasts,
                                baseline_type, n_projs, save_single_stcs)
    
if function_to_run == 'morph_beamformer_contrast_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    input_file_morph = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    proj = bool(int(argv[8]))
    h_freqs = list(map(int, list_parser(argv[9])))
    l_freqs = list(map(int, list_parser(argv[10])))
    tmin = float(argv[11])
    tmax = float(argv[12])
    spacing = float(argv[13])
    reg = float(argv[14])
    weight_norm = argv[15]
    channel_type = argv[16]
    ISI_type = argv[17]
    contrasts = list_of_list_parser(argv[18])
    baseline = argv[19]
    save_single_stcs = bool(int(argv[20]))

    morph_beamformer_contrast_hilbert(subject, date, input_file,
                                      input_file_morph, output_file, overwrite,
                                      proj, h_freqs, l_freqs, tmin, tmax,
                                      spacing,
                                      reg, weight_norm, channel_type, ISI_type,
                                      contrasts, baseline, save_single_stcs)   
    
if function_to_run == 'grand_average_beamformer_contrast_hilbert':
    these_subjects = subjects#subjects[:10] + subjects[11:]
    these_dates = dates# dates[:10] + dates[11:]
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    proj = bool(int(argv[5]))
    h_freqs = list(map(int, list_parser(argv[6])))
    l_freqs = list(map(int, list_parser(argv[7])))
    tmin = float(argv[8])
    tmax = float(argv[9])
    spacing = float(argv[10])
    reg = float(argv[11])
    weight_norm = argv[12]
    channel_type = argv[13]
    ISI_type = argv[14]
    contrasts = list_of_list_parser(argv[15])
    baseline = argv[16]
    save_std = bool(int(argv[17]))

    grand_average_beamformer_contrast_hilbert(these_subjects, these_dates,
                                              input_file, output_file,
                                              overwrite, proj,
                                              h_freqs, l_freqs,
                                              tmin, tmax, spacing, reg,
                                              weight_norm, channel_type,
                                              ISI_type, contrasts,
                                              baseline, save_std)
    
if function_to_run == 'statistics_beamformer_contrast_hilbert':
    these_subjects = subjects
    these_dates = dates
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    proj = bool(int(argv[5]))
    h_freqs = list(map(int, list_parser(argv[6])))
    l_freqs = list(map(int, list_parser(argv[7])))
    spacing = float(argv[8])
    reg = float(argv[9])
    weight_norm = argv[10]
    channel_type = argv[11]
    contrasts = list_of_list_parser(argv[12])
    p_threshold = float(argv[13])
    tmin = float(argv[14])
    tmax = float(argv[15])
    stat_tmin = float(argv[16])
    stat_tmax = float(argv[17])
    Type = argv[18]
    ISI_type = argv[19]
    stat_function = argv[20]
    baseline = argv[21]
    n_jobs = int(argv[22])

    statistics_beamformer_contrast_hilbert(these_subjects, these_dates,
                                           input_file, output_file, overwrite,
                                           proj, h_freqs, l_freqs, spacing,
                                           reg, weight_norm, channel_type,
                                           contrasts,
                                           p_threshold, tmin, tmax,
                                           stat_tmin, stat_tmax, Type,
                                           ISI_type, stat_function,
                                           baseline, n_jobs)
    
if function_to_run == 'save_nifti':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    tmin = float(argv[9])
    tmax = float(argv[10])
    spacing = float(argv[11])
    reg = float(argv[12])
    weight_norm = argv[13]
    channel_type = argv[14]
    ISI_type = argv[15]
    contrasts = list_of_list_parser(argv[16])
    baseline = argv[17]


    
    save_nifti(subject, date, input_file, output_file, overwrite, h_freqs,
               l_freqs, tmin, tmax,  spacing, reg, weight_norm, channel_type,
               ISI_type, contrasts, baseline)
    
    
if function_to_run == 'save_nifti_stat':
    input_file_1 = argv[2]
    input_file_2 = argv[3]
    output_file = argv[4]
    overwrite = bool(int(argv[5]))
    h_freq = int(argv[6])
    l_freq = int(argv[7])
    spacing = float(argv[8])
    reg = float(argv[9])
    weight_norm = argv[10]
    channel_type = argv[11]
    contrasts = list_of_list_parser(argv[12])
    p_threshold = float(argv[13])
    tmin = float(argv[14])
    tmax = float(argv[15])
    stat_tmin = float(argv[16])
    stat_tmax = float(argv[17])
    Type = argv[18]
    ISI_type = argv[19]
    stat_function = argv[20]
    
    save_nifti_stat(input_file_1, input_file_2, output_file, overwrite, h_freq,
                    l_freq, spacing, reg, weight_norm, channel_type, contrasts,
                    p_threshold, tmin, tmax, stat_tmin, stat_tmax, Type,
                    ISI_type, stat_function)
    
if function_to_run == 'read_and_resave_epochs':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    ISI_type = argv[9]
    tmin = float(argv[10])
    tmax = float(argv[11])
    full_baseline = bool(int(argv[12]))
    
    read_and_resave_epochs(subject, date, input_file, output_file, overwrite,
                           h_freqs, l_freqs, ISI_type, tmin, tmax,
                           full_baseline)
    
    
if function_to_run == 'plot_alignment':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[7]
    overwrite = bool(int(argv[7]))
    
    plot_alignment(subject, date, split_recording, input_file, output_file,
                   overwrite)    
    
if function_to_run == 'plot_evoked_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    
    plot_evoked_hilbert(subject, date, input_file, output_file, overwrite,
                        h_freqs, l_freqs)
    
if function_to_run == 'plot_projs_evoked_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    ISI_type = argv[9]
    
    plot_projs_evoked_hilbert(subject, date, input_file, output_file, 
                              overwrite, h_freqs, l_freqs, ISI_type)    
    
if function_to_run == 'plot_all_together_evoked_hilbert':
    these_subjects = subjects
    these_dates = dates
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    channel_type = argv[7]
    ISI_type = argv[8]
    
    plot_all_together_evoked_hilbert(these_subjects, these_dates, input_file,
                                     output_file, overwrite, h_freqs, l_freqs,
                                     channel_type, ISI_type)
    
if function_to_run == 'plot_beamformer_contrast_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[6])))
    l_freqs = list(map(int, list_parser(argv[7])))
    tmin = float(argv[8])
    tmax = float(argv[9])
    spacing = float(argv[10])
    reg = float(argv[11])
    weight_norm = argv[12]
    channel_type = argv[13]
    ISI_type = argv[14]
    contrasts = list_of_list_parser(argv[15])
    initial_time = float(argv[16])
    initial_pos = list(map(float, list_parser(argv[17])))
    
    plot_beamformer_contrast_hilbert(subject, date, input_file, output_file,
                                     overwrite, h_freqs, l_freqs, tmin, tmax,
                                     spacing, reg, weight_norm, channel_type,
                                     ISI_type, contrasts,
                                     initial_time, 
                                     initial_pos)
    
if function_to_run == 'plot_beamformer_time_courses_hilbert':
    subject = argv[2]    
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    proj = bool(int(argv[7]))
    h_freqs = list(map(int, list_parser(argv[8])))
    l_freqs = list(map(int, list_parser(argv[9])))
    tmin = float(argv[10])
    tmax = float(argv[11])
    spacing = float(argv[12])
    reg = float(argv[13])
    weight_norm = argv[14]
    channel_type = argv[15]
    ISI_type = argv[16]
    contrasts = list_of_list_parser(argv[17])
    baseline = argv[18]
    roi = argv[19]
    
    plot_beamformer_time_courses_hilbert(subject, date, input_file,
                                         output_file, overwrite, proj,
                                         h_freqs, l_freqs, tmin, tmax,
                                         spacing, reg, weight_norm, 
                                         channel_type, ISI_type, contrasts,
                                         baseline, roi)

    
if function_to_run == 'plot_whitened_data_cov':
    subject = argv[2]
    date = argv[3]
    input_file_1 = argv[4]
    input_file_2 = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    collapse_cond = bool(int(argv[8]))
    whiten = bool(int(argv[9]))
    h_freqs = list(map(int, list_parser(argv[10])))
    l_freqs = list(map(int, list_parser(argv[11])))
    spacing = float(argv[12])
    reg = float(argv[13])
    weight_norm = argv[14]
    channel_type = argv[15]
    contrasts = list_of_list_parser(argv[16])

    plot_whitened_data_cov(subject, date, input_file_1, input_file_2,
                           output_file, overwrite, collapse_cond, whiten,
                           h_freqs, l_freqs, spacing, reg, weight_norm,
                           channel_type, contrasts)    