#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:22:20 2020

@author: lau
"""

#%% IMPORTS

import mne
import numpy as np
from scipy import stats as stats
from scipy import linalg
from subprocess import check_output
from os.path import join, isfile
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
print(function_to_run)
print(argv[:])

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

def collapse_conditions(epochs, ISI_type, full_baseline=False):
    to_collapse = {}
    to_collapse['ns']    = ['n1_0',  'n2_0',  'n3_0',  'n4_0',  'n5_0',
                            'n1_5',  'n2_5',  'n3_5',  'n4_5',  'n5_5',
                            'n1_15', 'n2_15', 'n3_15', 'n4_15', 'n5_15'
                            ]
    if ISI_type != 'local' and ISI_type != 'mean' and not full_baseline:
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
    string = condition + '_' + str(spacing) + '_mm_common_filter_' + \
             contrast[0] + '_' + contrast[1] + '_'
    return string

def reg_and_weight_norm_string(reg, weight_norm):
    string = 'reg_' + str(reg) + '_weight_' + weight_norm + '_'
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
                                                      full_baseline=False):
    
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
    if full_baseline:
        baseline_text = 'full_baseline_'
    else:
        baseline_text = ''
    
    input_file_1 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[0],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_string(tmin, tmax) + proj + org_name_input
    input_file_2 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[1],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_string(tmin, tmax) + proj +  org_name_input
    input_file_3 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_contrast_string(spacing, contrast) + \
                        hp_lp_string(h_freq, l_freq) + \
                                time_string(tmin, tmax) + proj + org_name_input
    output_file_1 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[0],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_string(tmin, tmax) + proj + \
                                    org_name_output
    output_file_2 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_condition_common_filter_string(spacing, contrast[1],
                                               contrast) + \
                            hp_lp_string(h_freq, l_freq) + \
                                time_string(tmin, tmax) + proj + \
                                    org_name_output
    output_file_3 = \
        ISI_string + baseline_text + channel_type + '_' + \
            reg_and_weight_norm_string(reg, weight_norm) + \
            spacing_contrast_string(spacing, contrast) + \
                        hp_lp_string(h_freq, l_freq)  + \
                                time_string(tmin, tmax) + proj + \
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
        save_path_5 = join(get_fig_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', output_file_3)
        save_path_6 = join(get_fig_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', 'contrasts', output_file_3)
    else:
        save_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_1)
        save_path_2 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', output_file_2)
        save_path_3 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_1)
        save_path_4 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', output_file_2)
        save_path_5 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                           'common_filter', 'contrasts', output_file_3)
        save_path_6 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                           'common_filter', 'contrasts', output_file_3)
            
    
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
            raw = handle_split_recording(subject, date, input_file)

        raw.notch_filter(freqs=np.arange(50, 251, 50))
        raw.save(save_path, overwrite=overwrite)
        
def epoch_data_tfr(subject, date, split_recording, input_file, output_file,
                   overwrite, events_file):
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
        tmin = -1.500
        tmax = 1.500
        baseline = (None, None)
        decim = 1
        reject =  dict(grad=4000e-13,
               mag=4e-12)

        epochs = mne.Epochs(raw, events, event_id, tmin, tmax, baseline,
                            proj=False, decim=decim, reject=reject)
        epochs.save(save_path, overwrite=overwrite)
        
def run_ICA(subject, date, input_file, output_file, overwrite,
            h_freq, l_freq):
    input_file = hp_lp_string(h_freq, l_freq) + input_file
    output_file = hp_lp_string(h_freq, l_freq) + output_file
    load_path = join(get_save_path(subject, date), input_file)

    save_path = join(get_save_path(subject, date), output_file)
    
    if should_we_run(save_path, overwrite):
        ica = mne.preprocessing.ICA(n_components=0.95, n_pca_components=None,
                                    max_pca_components=None)
        raw = mne.io.read_raw_fif(load_path, preload=True)
        ica.fit(raw)
        nMaxVog, nMaxHog, nMaxEkg = [1, 1, 3]

        vogIndices, vogScores = ica.find_bads_eog(raw, ch_name='EOG001')
        hogIndices, hogScores = ica.find_bads_eog(raw, ch_name='EOG002')
        ekgIndices, ekgScores = ica.find_bads_ecg(raw, ch_name='ECG003',
                                                  method='ctps')
    
        vogIndices = vogIndices[:nMaxVog]
        hogIndices = hogIndices[:nMaxHog]
        ekgIndices = ekgIndices[:nMaxEkg]
    
        ica.exclude += vogIndices
        ica.exclude += hogIndices
        ica.exclude += ekgIndices  
    
        ica.save(save_path)
        
def get_tfr(subject, date, input_file, output_file, overwrite, collapse_cond):
    load_path = join(get_save_path(subject, date), input_file)
    save_path = join(get_save_path(subject, date), output_file)
    if should_we_run(save_path, overwrite):
        epochs = mne.read_epochs(load_path, proj=False)
        if collapse_cond:
            epochs = collapse_conditions(epochs)
        tfrs = []
        freqs = np.arange(1, 41, 1)
        n_cycles = freqs
        for event in epochs.event_id:
            print('TFR on event: ' + event)
            tfr, itc = mne.time_frequency.tfr_multitaper(epochs[event], freqs,
                                                         n_cycles, n_jobs=-1)
            tfr.comment = event
            tfrs.append(tfr)

        mne.time_frequency.write_tfrs(save_path, tfrs, overwrite)

def grand_average_tfr(subjects, dates, input_file, output_file, overwrite,
                      low_memory):
    save_path = join(get_save_path('grand_averages', ''), output_file)
    if should_we_run(save_path, overwrite):
        if low_memory:
            print('Running low memory version')
            tfrs_dict = dict()
            n_subjects = len(subjects)
            for subject_index, subject in enumerate(subjects):
                date = dates[subject_index]
                load_path = join(get_save_path(subject, date), input_file)
                tfrs = mne.time_frequency.read_tfrs(load_path)
                for tfr in tfrs:
                    event = tfr.comment
                    print('Loading event: ' + event + ' for subject: ' + \
                          subject)
                    if low_memory:
                        if subject_index == 0:
                            tfrs_dict[event] = tfr.copy()
                        else:
                            if np.all(np.isnan(tfr.data)): ## FIXME
                                n_subjects -= 1
                                continue
                            tfrs_dict[event]._data += tfr.data
                    else:
                        if subject_index == 0:
                            tfrs_dict[event] = []
                        tfrs_dict[event].append(tfr)
                tfrs = None ## free up memory

            ga_tfrs = []
            if low_memory:
                for event in tfrs_dict:
                    tfrs_dict[event]._data /= n_subjects
                    tfrs_dict[event].comment = \
                            event + ": Grand average (n = %d)" % n_subjects
                    ga_tfrs.append(tfrs_dict[event])
            else:
                for event in tfrs_dict:
                    temp = tfrs_dict[event]
                    ga = mne.grand_average(temp)
                    ga.comment = event
                    ga_tfrs.append(ga)

            mne.time_frequency.write_tfrs(save_path, ga_tfrs,
                                          overwrite=overwrite)      
            
def beamformer(subject, date, input_file, output_file, overwrite,
               collapse_cond, fmin, fmax, tmin, tmax, conditions):
    print('Running on condicions: ' + str(conditions))
    load_path = join(get_save_path(subject, date), input_file)
    org_name_output = output_file
    epochs = mne.read_epochs(load_path, proj=False, preload=False)
    epochs.del_proj()
    bem_folder = join(subjects_dir, subject, 'bem')

    ## BUILD FORWARD MODEL ON THE FLY

    ico_string = '5120'
    spacing = 5
    trans = join(bem_folder, subject + '-trans.fif')
    src = join(bem_folder, 'volume-' + str(spacing) + 'mm-src.fif')
    bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')

    fwd = None ## setting it later

    if collapse_cond:
        epochs = collapse_conditions(epochs)

    for condition in conditions:
        print('Running condition: ' + condition)
        output_file = freq_time_string_condition(condition, fmin, fmax,
                                                 tmin, tmax) + org_name_output
        save_path = join(get_save_path(subject, date), 'dics', output_file)
        if should_we_run(save_path, overwrite):
            these_epochs = epochs[condition]
            these_epochs.load_data()
            these_epochs.pick_types(meg='grad')

            ## CSD
            freqs = np.arange(fmin, fmax + 1, 1)
            csd = mne.time_frequency.csd_morlet(these_epochs, freqs,
                                                tmin=tmin, tmax=tmax)
            if fwd is None: ## only make once
                fwd = mne.make_forward_solution(these_epochs.info, trans, src,
                                                bem)

            ## DICS
            filters = mne.beamformer.make_dics(these_epochs.info, fwd,
                                               csd.mean(),
                                               pick_ori='max-power',
                                               inversion='single',
                                               weight_norm='unit-noise-gain',
                                               normalize_fwd=False)
            stc, freqs = mne.beamformer.apply_dics_csd(csd.mean(), filters)

            stc.save(save_path)

def beamformer_contrast(subject, date, input_file, output_file,
                        overwrite, collapse_cond, fmin, fmax, tmin, tmax,
                        contrasts):
    load_path = join(get_save_path(subject, date), input_file)
    org_name_output = output_file
    epochs = mne.read_epochs(load_path, proj=False, preload=False)
    epochs.del_proj()
    bem_folder = join(subjects_dir, subject, 'bem')

    ## BUILD FORWARD MODEL ON THE FLY
    ico_string = '5120'
    spacing = 5
    trans = join(bem_folder, subject + '-trans.fif')
    src = join(bem_folder, 'volume-' + str(spacing) + 'mm-src.fif')
    bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')

    fwd = None ## built later

    if collapse_cond:
        epochs = collapse_conditions(epochs)

    for contrast_index, contrast in enumerate(contrasts):
        print('Running contrast: ' + contrast[0] + ' versus ' + contrast[1])
        output_file = freq_time_string_contrast(contrast, fmin, fmax,
                                                tmin, tmax) + \
                                                org_name_output
        save_path = join(get_save_path(subject, date), 'dics_contrasts',
                         output_file)
        if should_we_run(save_path, overwrite):

            combined = epochs[contrast]
            first  = epochs[contrast[0]]
            second = epochs[contrast[1]]

            combined.load_data()
            first.load_data()
            second.load_data()

            combined.pick_types(meg='grad')
            first.pick_types(meg='grad')
            second.pick_types(meg='grad')

            ## CSDS

            freqs = np.arange(fmin, fmax + 1 , 1)
            csd_combined = mne.time_frequency.csd_morlet(combined, freqs,
                                                         tmin=tmin, tmax=tmax)
            csd_first    = mne.time_frequency.csd_morlet(first, freqs,
                                                         tmin=tmin, tmax=tmax)
            csd_second   = mne.time_frequency.csd_morlet(second, freqs,
                                                         tmin=tmin, tmax=tmax)
            ## DICS
            
            if fwd is None: ## only make it once
                fwd = mne.make_forward_solution(combined.info, trans, src,
                                                 bem)

            filters = mne.beamformer.make_dics(combined.info, fwd,
                                               csd_combined.mean(),
                                               pick_ori='max-power')
            stc_first, freqs = mne.beamformer.apply_dics_csd(csd_first.mean(),
                                                             filters)
            stc_second, freqs = \
                mne.beamformer.apply_dics_csd(csd_second.mean(),
                                                             filters)
            stc_contrast = (stc_first - stc_second) / (stc_first + stc_second)

            stc_contrast.save(save_path)

def morph_beamformer(subject, date, input_file, input_file_morph,
                     output_file, overwrite, fmin, fmax, tmin, tmax,
                     conditions):
    org_name_input = input_file
    org_name_output = output_file
    load_path_morph = join(get_save_path(subject, date), input_file_morph)
    morph = mne.read_source_morph(load_path_morph)

    for condition in conditions:
        print('Running condition: ' + condition)
        input_file = freq_time_string_condition(condition, fmin, fmax,
                                               tmin, tmax) + \
                                                org_name_input
        output_file = freq_time_string_condition(condition, fmin, fmax,
                                                tmin, tmax) + \
                                                org_name_output

        load_path = join(get_save_path(subject, date), 'dics',
                         input_file)
        save_path = join(get_save_path(subject, date), 'dics',
                         output_file)
        if should_we_run(save_path, overwrite):
            stc = mne.read_source_estimate(load_path)
            stc_morph = morph.apply(stc)
            stc_morph.save(save_path)

def morph_beamformer_contrast(subject, date, input_file, input_file_morph,
                              output_file, overwrite, fmin, fmax, tmin, tmax,
                              contrasts):
    org_name_input = input_file
    org_name_output = output_file
    load_path_morph = join(get_save_path(subject, date), input_file_morph)
    morph = mne.read_source_morph(load_path_morph)

    for contrast in contrasts:
        print('Running contrast: ' + contrast[0] + ' versus ' + contrast[1])
        input_file = freq_time_string_contrast(contrast, fmin, fmax,
                                               tmin, tmax) + \
                                                org_name_input
        output_file = freq_time_string_contrast(contrast, fmin, fmax,
                                                tmin, tmax) + \
                                                org_name_output

        load_path = join(get_save_path(subject, date), 'dics_contrasts',
                         input_file)
        save_path = join(get_save_path(subject, date), 'dics_contrasts',
                         output_file)
        if should_we_run(save_path, overwrite):
            stc = mne.read_source_estimate(load_path)
            stc_morph = morph.apply(stc)
            stc_morph.save(save_path)

def grand_average_beamformer(subjects, dates, input_file,
                             output_file, overwrite, fmin, fmax,
                             tmin, tmax, conditions):
    org_name_input = input_file
    org_name_output = output_file

    stcs_dict = dict()
    for subject_index, subject in enumerate(subjects):
        date = dates[subject_index]
        for condition in conditions:
            print('Running condition: ' + condition)
            input_file = freq_time_string_condition(condition, fmin, fmax,
                                                    tmin, tmax) + \
                                                    org_name_input
            load_path = join(get_save_path(subject, date), 'dics',
                             input_file)
            print('Loading: ' + str(condition) + ' for subject: ' + subject)
            stc = mne.read_source_estimate(load_path)
            if subject_index == 0:
                stcs_dict[condition] = []
            stcs_dict[condition].append(stc)

    ga_stcs = dict()
    for condition in stcs_dict:
        stcs = sum(stcs_dict[condition]) / len(stcs_dict[condition])
        ga_stcs[condition] = stcs

    for condition in conditions:
        output_file = freq_time_string_condition(condition, fmin, fmax,
                                                tmin, tmax) + \
                                                  org_name_output
        save_path = join(get_save_path('grand_averages', ''),
                         'dics', output_file)
        stc = ga_stcs[condition]
        stc.save(save_path)


def grand_average_beamformer_contrast(subjects, dates, input_file,
                                      output_file, overwrite, fmin, fmax,
                                      tmin, tmax, contrasts):
    org_name_input = input_file
    org_name_output = output_file

    stcs_dict = dict()
    for subject_index, subject in enumerate(subjects):
        date = dates[subject_index]
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])
            contrast_name = contrast[0] + '_vs_' + contrast[1]
            input_file = freq_time_string_contrast(contrast, fmin, fmax,
                                                   tmin, tmax) + \
                                                  org_name_input
            load_path = join(get_save_path(subject, date), 'dics_contrasts',
                             input_file)
            print('Loading: ' + str(contrast) + ' for subject: ' + subject)
            stc = mne.read_source_estimate(load_path)
            if subject_index == 0:
                stcs_dict[contrast_name] = []
            stcs_dict[contrast_name].append(stc)

    ga_stcs = dict()
    for contrast in stcs_dict:
        stcs = sum(stcs_dict[contrast]) / len(stcs_dict[contrast])
        ga_stcs[contrast] = stcs

    for contrast in contrasts:
        contrast_name = contrast[0] + '_vs_' + contrast[1]
        output_file = freq_time_string_contrast(contrast, fmin, fmax,
                                                tmin, tmax) + \
                                                  org_name_output
        save_path = join(get_save_path('grand_averages', ''),
                         'dics_contrasts', output_file)
        stc = ga_stcs[contrast_name]
        stc.save(save_path)
        
def statistics_beamformer_contrast(subjects, dates, input_file, output_file,
                                   overwrite, fmin, fmax, tmin, tmax,
                                   contrasts):
    org_name_input = input_file
    org_name_output = output_file
    
    bem_folder = join(subjects_dir, 'fsaverage', 'bem')
    spacing = 5
    src = \
        mne.source_space.read_source_spaces(join(bem_folder,
                                                 'volume-' + str(spacing) + \
                                                     'mm-src.fif'))
    connectivity = mne.spatial_src_connectivity(src)
    n_subjects = len(subjects)
    p_threshold = 0.01
    t_threshold = -stats.distributions.t.ppf(p_threshold / 2.0, n_subjects - 1)

    stcs_dict = dict()
    for subject_index, subject in enumerate(subjects):
        date = dates[subject_index]
        for contrast in contrasts:
            print('Running contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1])
            contrast_name = contrast[0] + '_vs_' + contrast[1]
            input_file = freq_time_string_contrast(contrast, fmin, fmax,
                                                   tmin, tmax) + \
                                                  org_name_input
            load_path = join(get_save_path(subject, date), 'dics_contrasts',
                             input_file)
            print('Loading: ' + str(contrast) + ' for subject: ' + subject)
            stc = mne.read_source_estimate(load_path)
            if subject_index == 0:
                stcs_dict[contrast_name] = []
            stcs_dict[contrast_name].append(stc)
        
    stat_dict = dict()
    for contrast in stcs_dict:
        stat_array = None
        stcs = stcs_dict[contrast]
        for stc in stcs:
            if stat_array is None:
                stat_array = stc.data.T
                stat_array = np.expand_dims(stat_array, 0)
            else:
                temp = stc.data.T
                temp = np.expand_dims(temp, 0)
                stat_array = np.concatenate((stat_array, temp), 0)
        print('Clustering for contrast: ' + contrast)
        T_obs, clusters, cluster_p_values, H0 = clu = \
            mne.stats.spatio_temporal_cluster_1samp_test(stat_array,
                connectivity=connectivity, n_jobs=-1,
                threshold=t_threshold, verbose=True)
        stat_dict[contrast] = clu
        
    ## save dicts
    for contrast in contrasts:
        contrast_name = contrast[0] + '_vs_' + contrast[1]
        output_file = freq_time_string_contrast(contrast, fmin, fmax,
                                            tmin, tmax) + \
                                              org_name_output
        save_path = join(get_save_path('grand_averages', ''),
                     'statistics', output_file)
        clu = stat_dict[contrast_name]
        cluster_dict = dict()
        names = ['T_obs', 'clusters', 'cluster_p_values', 'H0']
        for name_index, name in enumerate(names):
            cluster_dict[name] = clu[name_index]
        np.save(save_path, cluster_dict)
        
def z_transform_wilcoxon(subject, date, input_file,
                         output_file, overwrite, collapse_cond, proj, h_freqs,
                         l_freqs, ISI_type):
    
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
    if proj:
        proj_string = 'proj_'
    else:
        proj_string = ''
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + proj_string + \
                     org_name_output
        load_path = join(get_save_path(subject, date), 'hilbert',
                           input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        if should_we_run(save_path, overwrite):
            epochs_hilbert = mne.read_epochs(load_path, proj=proj)
            if not proj:
                epochs_hilbert.del_proj()
            picks = mne.pick_types(epochs_hilbert.info, meg=True, eog=False,
                                   ecg=False, emg=False, misc=False)
            epochs_hilbert.pick(picks)
            if collapse_cond:
                epochs_hilbert = collapse_conditions(epochs_hilbert, ISI_type)
            z_transformed = []
 
            for event in epochs_hilbert.event_id:
                print('Transforming event: ' + event)
                data = epochs_hilbert[event].get_data()
                n_channels = data.shape[1]
                n_samples = data.shape[2]
                zs = np.zeros((n_channels, n_samples))
            
                for channel_index in range(n_channels):
                    for sample_index in range(n_samples):
                        this_data = data[:, channel_index, sample_index]
                        z = wilcoxon(this_data)
                        zs[channel_index, sample_index] = z
                z_transform = mne.EvokedArray(zs, info=epochs_hilbert.info,
                                              tmin=epochs_hilbert.tmin,
                                              comment=event,
                                              nave=this_data.shape[0])
                z_transformed.append(z_transform)
                
            mne.write_evokeds(save_path, z_transformed)
        epochs_hilbert = None ## memory control        
        
def tfr_z_transform(subject, date, input_file, output_file, overwrite,
                    h_freqs, l_freqs, ISI_type):
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_output
        load_path = join(get_save_path(subject, date), 'hilbert',
                           input_file)
        save_path = join(get_save_path(subject, date), 'hilbert', output_file)
        freqs = np.arange(h_freq, l_freq + 1, 1)
        if should_we_run(save_path, overwrite):
            evokeds = mne.read_evokeds(load_path)
            tfrs = []
            for evoked in evokeds:
                tfr = mne.time_frequency.tfr_multitaper(evoked, freqs, 3,
                                                        return_itc=False)
                tfr.comment = evoked.comment
                tfrs.append(tfr)
            mne.time_frequency.write_tfrs(save_path, tfrs, overwrite)
            
def grand_average_tfr_z_transform(subjects, dates, input_file, output_file,
                                  overwrite, h_freqs, l_freqs, ISI_type):
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
    n_subjects = len(subjects)
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_output
        
        tfrs_dict = dict()
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            load_path = join(get_save_path(subject, date), 'hilbert',
                             input_file)

            print('Loading tfr for subject: ' + subject)
            tfrs = mne.time_frequency.read_tfrs(load_path)
            for tfr in tfrs:
                if subject_index == 0:
                    tfrs_dict[tfr.comment] = []
                tfrs_dict[tfr.comment].append(tfr)
                
        ga_tfrs = []
        for event in tfrs_dict:
            print('Getting grand average for: ' + event)
            ga = mne.grand_average(tfrs_dict[event])
            ga.comment = event + ": Grand average (n = %d)" % n_subjects
            ga_tfrs.append(ga)
        save_path = join(get_save_path('grand_averages', ''),
                                       'hilbert', output_file)
        mne.time_frequency.write_tfrs(save_path, ga_tfrs, overwrite)
        
def grand_average_evoked_hilbert_contrast(subjects, dates, input_file,
                                          output_file, overwrite, h_freqs,
                                          l_freqs, contrasts):
    org_name_input = input_file
    org_name_output = output_file
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = hp_lp_string(h_freq, l_freq) + org_name_output
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

        ga_contrasts = []
        for contrast in contrasts:
            contrast_name = contrast[0] + '_vs_' + contrast[1]
            contrast_1 = evokeds_dict[contrast[0]]
            contrast_2 = evokeds_dict[contrast[1]]
            temp_list = []
            for (temp_1, temp_2) in zip(contrast_1, contrast_2):
                temp = temp_1.copy()
                temp._data = temp_1.data - temp_2.data
                temp_list.append(temp)
            ga = mne.grand_average(temp_list)
            ga.comment = contrast_name
            ga_contrasts.append(ga)

        mne.evoked.write_evokeds(save_path, ga_contrasts)        
                
            
            
def statistics_tfr_z_transform(sybjects, dates, input_file, output_file,
                               overwrite, h_freqs, l_freqs, contrast,
                               ISI_type, n_jobs):
                
    org_name_input = input_file
    org_name_output = output_file
    prepend = ISI_type + '_'
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_input
        output_file = prepend + hp_lp_string(h_freq, l_freq) + org_name_output
        n_subjects = len(subjects)
        
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            load_path = join(get_save_path(subject, date), 'hilbert',
                             input_file)

            print('Loading contrast: ' + contrast[0] + ' versus ' + \
                  contrast[1] + ' for subject: ' + subject)
            tfrs = []
            all_tfrs = mne.time_frequency.read_tfrs(load_path)
            for tfr in all_tfrs:
                if tfr.comment in contrast:
                    tfr.pick_types('mag')
                    tfrs.append(tfr)
            all_tfrs = None
            tfr = None
            if subject_index == 0:
                n_channels, n_freqs, n_samples = tfrs[0].data.shape
                array = np.zeros(shape=(n_subjects, n_channels, n_freqs,
                                          n_samples))                                                   
            array[subject_index, :, :, :] = \
                (tfrs[0].data - tfrs[1].data) / (tfrs[0].data + tfrs[1].data)
        
        ## test
        threshold = \
                        -stats.distributions.norm.ppf(0.05 / 2.0)
        t_obs, clusters, cluster_p_values, H0 = clu = \
            mne.stats.permutation_cluster_1samp_test(array, threshold,
                                                     n_permutations=10,
                                                     out_type='indices',
                                                     seed=7, n_jobs=n_jobs)
        cluster_dict = dict()
        names = ['t_obs', 'clusters', 'cluster_p_values', 'H0']
        for name_index, name in enumerate(names):
            cluster_dict[name] = clu[name_index]
        save_path = join(get_save_path('grand_averages', ''), 'statistics',
                                       'hilbert', output_file)
        np.save(save_path, cluster_dict)
        
def beamformer_hilbert(subject, date,
                       input_file_1, input_file_2, output_file_1,
                       output_file_2, overwrite, collapse_cond,
                       h_freqs, l_freqs, spacing, conditions):
    org_name_input_1 = input_file_1
    org_name_input_2 = input_file_2
    org_name_output_1 = output_file_1
    org_name_output_2 = output_file_2
    bem_folder = join(subjects_dir, subject, 'bem')
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        input_file_1 = hp_lp_string(h_freq, l_freq) + org_name_input_1
        input_file_2 = hp_lp_string(h_freq, l_freq) + org_name_input_2
        load_path_1 = join(get_save_path(subject, date), 'hilbert',
                           input_file_1)
        load_path_2 = join(get_save_path(subject, date), 'hilbert',
                           input_file_2)
   
        raw = mne.io.read_raw_fif(load_path_1, preload=False)
        epochs_hilbert = mne.read_epochs(load_path_2, preload=False,
                                         proj=False) ## preload mus be True
                                                     ## due to bug with repre-
                                                     ## senting them correctly
                                                     ## as complex128 (0.20.0)
        if collapse_cond:
            epochs_hilbert = collapse_conditions(epochs_hilbert)
        
        
        ico_string = '5120'
        trans = join(bem_folder, subject + '-trans.fif')
        src = join(bem_folder, 'volume-' + str(spacing) + 'mm-src.fif')
        bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')

        fwd = mne.make_forward_solution(epochs_hilbert.info, trans, src, bem)
        for condition in conditions:
            print('Running condition: ' + condition)
            event_id = epochs_hilbert.event_id[condition]            
            output_file_1 = condition + '_' + str(spacing) + '_mm_' + \
                            hp_lp_string(h_freq, l_freq) + org_name_output_1
            output_file_2 = condition + '_' + str(spacing) + '_mm_' + \
                            hp_lp_string(h_freq, l_freq) + org_name_output_2
            save_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                       output_file_1)
            save_path_2 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                       output_file_2)
            if should_we_run(save_path_1, overwrite) or \
                should_we_run(save_path_2, overwrite):
                events = epochs_hilbert.events
                tmin = epochs_hilbert.tmin
                tmax = epochs_hilbert.tmax
                baseline = epochs_hilbert.baseline
                decim = 1
                reject =  epochs_hilbert.reject
        
                epochs_cov = mne.Epochs(raw, events, event_id, tmin, tmax,
                                        baseline, proj=False, decim=decim,
                                        reject=reject)
                epochs_cov.del_proj()
                if 'sss' in input_file_1:
                    rank = 'info'
                else:
                    rank = None
                data_cov = mne.compute_covariance(epochs_cov, tmin=0,
                                                  tmax=None, rank=rank)
                noise_cov = mne.compute_covariance(epochs_cov, tmin=None,
                                                   tmax=baseline[1],
                                                   rank=rank)
                
                filters = mne.beamformer.make_lcmv(epochs_cov.info, fwd,
                                                   data_cov=data_cov,
                                                   noise_cov=noise_cov,
                                                   pick_ori='max-power',
                                                   weight_norm='nai',
                                                   reg=0.05)
                epochs_cov = None
                stcs = \
                    mne.beamformer.apply_lcmv_epochs(epochs_hilbert[condition],
                                                     filters,
                                                     max_ori_out='signed')
                              
                ## itcs
                itc = compute_source_itc(stcs)
                itc_stc = stcs[0].copy()
                itc_stc._data = itc
                ## Hilbert stc
                for stc in stcs:
                    stc._data = np.array(np.abs(stc.data),
                                              dtype='float64')
                    stc.resample(250.0)
                    
                ## get means
                stc_mean = stcs[0].copy()
                mean_data = np.mean([stc.data for stc in stcs], axis=0)
                stc_mean._data = mean_data
                
                stc_mean.save(save_path_1)
                itc_stc.save(save_path_2)
                
                ## memory control
                stcs = None
                drop_indices = epochs_hilbert.events[:, 2] == event_id
                epochs_hilbert.drop(drop_indices, reason='memory_control')
                
def morph_beamformer_hilbert(subject, date, input_file, input_file_morph,
                             output_file, overwrite,
                             h_freqs, l_freqs, spacing, conditions):
    org_name_input = input_file
    org_name_output = output_file
    input_file_morph = str(spacing) + '_mm_' + input_file_morph
    load_path_morph = join(get_save_path(subject, date), input_file_morph)
    morph = mne.read_source_morph(load_path_morph)
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):

        for condition in conditions:
            print('Running condition: ' + condition)
            input_file = condition + '_' + str(spacing) + '_mm_' + \
                            hp_lp_string(h_freq, l_freq) + org_name_input
   
            output_file = condition + '_' + str(spacing) + '_mm_' + \
                            hp_lp_string(h_freq, l_freq) + org_name_output

            load_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                             input_file)
            load_path_2 = join(get_save_path(subject, date), 'hilbert', 'itcs',
                             input_file)
            save_path_1 = join(get_save_path(subject, date), 'hilbert', 'stcs',
                             output_file)
            save_path_2 = join(get_save_path(subject, date), 'hilbert', 'itcs',  
                               output_file)
     
            if should_we_run(save_path_1, overwrite):
                stc = mne.read_source_estimate(load_path_1)
                stc_morph = morph.apply(stc)
                stc_morph.save(save_path_1)
            if should_we_run(save_path_2, overwrite):
                itc = mne.read_source_estimate(load_path_2)
                itc_morph = morph.apply(itc)
                itc_morph.save(save_path_2)                

def grand_average_beamformer_hilbert(subjects, dates, input_file,
                                     output_file, overwrite, h_freqs, l_freqs,
                                     spacing, conditions):
    org_name_input = input_file
    org_name_output = output_file

    stcs_dict = dict()
    itcs_dict = dict()
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        print('Running frequency range: ' + str(h_freq) + '-' + str(l_freq) + \
              ' Hz')
        for subject_index, subject in enumerate(subjects):
            date = dates[subject_index]
            for condition in conditions:
                output_file = condition + '_' + str(spacing) + '_mm_' + \
                    hp_lp_string(h_freq, l_freq) + org_name_output
                save_path_1 = join(get_save_path('grand_averages', ''),
                               'hilbert', 'stcs', output_file)
                save_path_2 = join(get_save_path('grand_averages', ''),
                               'hilbert', 'itcs', output_file)
                if not(should_we_run(save_path_1, overwrite) or \
                    should_we_run(save_path_2, overwrite)):
                    print('Not overwriting')
                    continue
                
                print('Running condition: ' + condition)
                input_file = condition + '_' + str(spacing) + '_mm_' + \
                    hp_lp_string(h_freq, l_freq) + org_name_input
                load_path_1 = join(get_save_path(subject, date), 'hilbert',
                                   'stcs', input_file)
                print('Loading stc for : ' + str(condition) + \
                      ' for subject: ' + subject)
                stc = mne.read_source_estimate(load_path_1)
                if subject_index == 0:
                    stcs_dict[condition] = []
                stcs_dict[condition].append(stc)
                ## itc part
                input_file = condition + '_' + str(spacing) + '_mm_' + \
                    hp_lp_string(h_freq, l_freq) + org_name_input
                load_path_2 = join(get_save_path(subject, date), 'hilbert',
                                   'itcs', input_file)
                print('Loading itc for : ' + str(condition) + \
                      ' for subject: ' + subject)
                itc = mne.read_source_estimate(load_path_2)
                if subject_index == 0:
                    itcs_dict[condition] = []
                itcs_dict[condition].append(itc)
    
        ga_stcs = dict()
        ga_itcs = dict()
        for condition in stcs_dict:
            stcs = sum(stcs_dict[condition]) / len(stcs_dict[condition])
            itcs = sum(itcs_dict[condition]) / len(itcs_dict[condition])
            ga_stcs[condition] = stcs
            ga_itcs[condition] = itcs
            
    
        for condition in stcs_dict:
            output_file = condition + '_' + str(spacing) + '_mm_' + \
                hp_lp_string(h_freq, l_freq) + org_name_output
            save_path_1 = join(get_save_path('grand_averages', ''),
                               'hilbert', 'stcs', output_file)
            stc = ga_stcs[condition]
            stc.save(save_path_1)
            
            save_path_2 = join(get_save_path('grand_averages', ''),
                               'hilbert', 'itcs', output_file)
            itc = ga_itcs[condition]
            itc.save(save_path_2)
            
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

if function_to_run == 'epoch_data_tfr':
    subject = argv[2]
    date = argv[3]
    split_recording = bool(int(argv[4]))
    input_file = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    events_file = argv[8]

    epoch_data_tfr(subject, date, split_recording, input_file, output_file,
                   overwrite,events_file)
    
if function_to_run == 'run_ICA':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freq = int(argv[7])
    l_freq = int(argv[8])   
    
    run_ICA(subject, date, input_file, output_file, overwrite, h_freq, l_freq)
    
if function_to_run == 'get_tfr':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))

    get_tfr(subject, date, input_file, output_file, overwrite, collapse_cond)

if function_to_run == 'grand_average_tfr':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    low_memory = bool(int(argv[5]))
    these_subjects = subjects[:12] + subjects[14:]
    these_dates = dates[:12] + dates[14:]

    grand_average_tfr(these_subjects, these_dates, input_file, output_file,
                      overwrite, low_memory)    
    
if function_to_run == 'beamformer':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))
    fmin = int(argv[8])
    fmax = int(argv[9])
    tmin = float(argv[10])
    tmax = float(argv[11])
    conditions = list_parser(argv[12])

    beamformer(subject, date, input_file, output_file, overwrite,
               collapse_cond, fmin, fmax, tmin, tmax, conditions)

if function_to_run == 'morph_beamformer':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    input_file_morph = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    fmin = int(argv[8])
    fmax = int(argv[9])
    tmin = float(argv[10])
    tmax = float(argv[11])
    conditions = list_parser(argv[12])

    morph_beamformer(subject, date, input_file, input_file_morph,
                     output_file, overwrite, fmin, fmax, tmin, tmax,
                     conditions)

if function_to_run == 'grand_average_beamformer':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    fmin = int(argv[5])
    fmax = int(argv[6])
    tmin = float(argv[7])
    tmax = float(argv[8])
    these_subjects = subjects
    these_dates = dates
    conditions = list_parser(argv[9])

    grand_average_beamformer(these_subjects, these_dates, input_file,
                             output_file, overwrite, fmin, fmax,
                             tmin, tmax, conditions)

if function_to_run == 'beamformer_contrast':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))
    fmin = int(argv[8])
    fmax = int(argv[9])
    tmin = float(argv[10])
    tmax = float(argv[11])
    contrasts = list_of_list_parser(argv[12])

    beamformer_contrast(subject, date, input_file, output_file, overwrite,
                        collapse_cond, fmin, fmax, tmin, tmax, contrasts)

if function_to_run == 'morph_beamformer_contrast':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    input_file_morph = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    fmin = int(argv[8])
    fmax = int(argv[9])
    tmin = float(argv[10])
    tmax = float(argv[11])
    contrasts = list_of_list_parser(argv[12])

    morph_beamformer_contrast(subject, date, input_file, input_file_morph,
                              output_file, overwrite, fmin, fmax, tmin, tmax,
                              contrasts)

if function_to_run == 'grand_average_beamformer_contrast':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    fmin = int(argv[5])
    fmax = int(argv[6])
    tmin = float(argv[7])
    tmax = float(argv[8])
    these_subjects = subjects
    these_dates = dates
    contrasts = list_of_list_parser(argv[9])

    grand_average_beamformer_contrast(these_subjects, these_dates, input_file,
                                      output_file, overwrite, fmin, fmax,
                                      tmin, tmax, contrasts)
    
if function_to_run == 'statistics_beamformer_contrast':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    fmin = int(argv[5])
    fmax = int(argv[6])
    tmin = float(argv[7])
    tmax = float(argv[8])
    these_subjects = subjects
    these_dates = dates
    contrasts = list_of_list_parser(argv[9])

    statistics_beamformer_contrast(these_subjects, these_dates, input_file,
                                   output_file, overwrite, fmin, fmax,
                                   tmin, tmax, contrasts)    
    
if function_to_run == 'tfr_z_transform':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    h_freqs = list(map(int, list_parser(argv[7])))
    l_freqs = list(map(int, list_parser(argv[8])))
    ISI_type = argv[9]
    
    tfr_z_transform(subject, date, input_file, output_file, overwrite,
                    h_freqs, l_freqs, ISI_type)
    
if function_to_run == 'grand_average_tfr_z_transform':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    ISI_type = argv[7]
    
    grand_average_tfr_z_transform(subjects, dates, input_file, output_file,
                                  overwrite, h_freqs, l_freqs, ISI_type)
    
if function_to_run == 'statistics_tfr_z_transform':
    these_subjects = subjects
    these_dates = dates
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    contrast = list_parser(argv[7])
    ISI_type = argv[8]
    n_jobs = int(argv[9])
    
    statistics_tfr_z_transform(subjects, dates, input_file, output_file,
                               overwrite, h_freqs, l_freqs, contrast, ISI_type,
                               n_jobs)  
    
if function_to_run == 'z_transform_wilcoxon':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    output_file = argv[5]
    overwrite = bool(int(argv[6]))
    collapse_cond = bool(int(argv[7]))
    proj = bool(int(argv[8]))
    h_freqs = list(map(int, list_parser(argv[9])))
    l_freqs = list(map(int, list_parser(argv[10])))
    ISI_type = argv[11]
    
    z_transform_wilcoxon(subject, date, input_file,
                         output_file, overwrite, collapse_cond, proj, h_freqs,
                         l_freqs, ISI_type) 
    
if function_to_run == 'grand_average_evoked_hilbert_contrast':
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    contrasts = list_of_list_parser(argv[7])
    
    grand_average_evoked_hilbert_contrast(subjects, dates, input_file,
                                          output_file, overwrite, h_freqs,
                                          l_freqs, contrasts)    
    
if function_to_run == 'beamformer_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file_1 = argv[4]
    input_file_2 = argv[5]
    output_file_1 = argv[6]
    output_file_2 = argv[7]
    overwrite = bool(int(argv[8]))
    collapse_cond = bool(int(argv[9]))
    h_freqs = list(map(int, list_parser(argv[10])))
    l_freqs = list(map(int, list_parser(argv[11])))
    spacing = float(argv[12])
    conditions = list_parser(argv[13])

    beamformer_hilbert(subject, date, input_file_1,
                       input_file_2, output_file_1, output_file_2, overwrite,
                       collapse_cond, h_freqs, l_freqs, spacing, conditions)   
    
if function_to_run == 'morph_beamformer_hilbert':
    subject = argv[2]
    date = argv[3]
    input_file = argv[4]
    input_file_morph = argv[5]
    output_file = argv[6]
    overwrite = bool(int(argv[7]))
    h_freqs = list(map(int, list_parser(argv[8])))
    l_freqs = list(map(int, list_parser(argv[9])))
    spacing = float(argv[10])
    conditions = list_parser(argv[11])

    morph_beamformer_hilbert(subject, date, input_file, input_file_morph,
                              output_file,  overwrite,
                              h_freqs, l_freqs, spacing, conditions)    
    
if function_to_run == 'grand_average_beamformer_hilbert':
    these_subjects = subjects[:26]
    these_dates = dates[:26]
    input_file = argv[2]
    output_file = argv[3]
    overwrite = bool(int(argv[4]))
    h_freqs = list(map(int, list_parser(argv[5])))
    l_freqs = list(map(int, list_parser(argv[6])))
    spacing = float(argv[7])
    conditions = list_parser(argv[8])

    grand_average_beamformer_hilbert(these_subjects, these_dates, input_file,
                                     output_file, overwrite, h_freqs, l_freqs,
                                     spacing, conditions)    