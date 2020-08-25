#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 09:16:00 2020

@author: lau
"""

#%% PLAY WITH NILEARN

from nilearn import image, datasets, plotting
from os import chdir
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import mne

script_path = '/home/lau/projects/cerebellar_clock/scripts/tactile_jitter/' + \
              'python/'
chdir(script_path)              

from qsubs import sensor_functions

scratch_path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/'
subjects_dir = join(scratch_path, 'freesurfer')

chdir('/home/lau/Skrivebord')
img = image.load_img('test.nii')
data = np.asanyarray(img.dataobj)

time_index = 751
img_time = image.index_img(img, time_index)

## source space and an stc

src = mne.read_source_spaces(join('/', 'home', 'lau', 'projects',
                             'cerebellar_clock', 'scratch', 'tactile_jitter',
                             'freesurfer','fsaverage', 'bem',
                                  'volume-7.5mm-src.fif'))

contrast = ['o0', 'o15']
stat_tmin, stat_tmax = (-0.400, 0.400)
tmin, tmax = (-0.750, 0.750)
freqs = [14, 30]
ISI_type = 'normal'
proj_string = ''#'proj_'
baseline_string =  'no_baseline_'


freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])

# contrast_name = contrast[0] + '_vs_' + contrast[1]
contrast_string = ''
for part in contrast:
    contrast_string += part + '_'
contrast_string += 'contrast_'    
    
time_string = sensor_functions.time_string(tmin, tmax)
data_path = join(scratch_path, 'MEG', 'grand_averages')
hilbert_path = join('hilbert', 'stcs', 'common_filter', 'contrasts')
filename = ISI_type + '_ISI_' + baseline_string + \
            'mag_reg_0.0_weight_unit-noise-gain_' + \
            contrast_string + '7.5_mm_' + freq_string + \
            time_string + proj_string + 'tactile_jitter_morph-vl.stc'
stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))



## labels of interest

rois = [
        'Cerebelum_6_L',
        'Cerebelum_6_R',
        'Cerebelum_8_L',
        'Cerebelum_8_R',
        'Hippocampus_L',
        'Hippocampus_R',
        'Thalamus_L',
        'Thalamus_R',
        'Putamen_L',
        'Putamen_R',
        'Precuneus_L',
        'Precuneus_R',
        'Insula_L',
        'Insula_R',
        'Postcentral_L',
        'Postcentral_R',
        'Parietal_Inf_L',
        'Parietal_Inf_R',
        'Cingulum_Mid_L',
        'Cingulum_Mid_R',
        
        ]

## aal atlas
atlas_dataset = datasets.fetch_atlas_aal()
atlas_filepath = atlas_dataset['maps']
atlas_img = image.load_img(atlas_filepath)
labels = atlas_dataset['labels']
indices = list(map(int, atlas_dataset['indices']))

## harvard oxford
# atlas_filename = 'sub-maxprob-thr0-1mm'
# atlas_dataset = datasets.fetch_atlas_harvard_oxford(atlas_filename)
# atlas_filepath = atlas_dataset.maps
# atlas_img = image.load_img(atlas_filepath)
# atlas_data = np.asanyarray(atlas_img.dataobj)
# labels = atlas_dataset['labels']
# n_labels = len(atlas_dataset['labels'])

## cerebellar atlas
# atlas_type = 'flirt'
# atlas_path = '/usr/local/fsl/data/atlases/'
# atlas_img = image.load_img(join(atlas_path, 'Cerebellum',
#                                   'Cerebellum-MNI' + atlas_type + \
#                                   '-maxprob-thr0-1mm.nii.gz'))
# atlas_data = np.asanyarray(atlas_img.dataobj)
# label_path = join(atlas_path, 'Cerebellum_MNI' + atlas_type + '.xml')
# label_file = open(label_path, 'r')
# labels = []
# for line in label_file:
#     if line[:6] != '<label':
#         continue
#     first_index = None
#     last_index = None
#     for character_index, character in enumerate(line):
#         if character == '>' and first_index is None:
#             first_index = character_index + 1
#         if character == '<' and first_index is not None:
#             last_index = character_index
#     if first_index is not None and last_index is not None:
#         labels.append(line[first_index:last_index])
# indices = range(1, len(labels))


## resampling
atlas_interpolated = image.resample_to_img(atlas_img, img, 'nearest')
atlas_interpolated_data = np.asanyarray(atlas_interpolated.dataobj)

## summarize time courses (mean)


times = np.arange(-0.750, 0.751, 0.001)
n_times = len(times)
n_labels = len(labels)

label_time_courses = np.zeros((n_labels, n_times))
for label_index, label in enumerate(indices):
    for time_index in range(n_times):
        this_data = data[:, :, :, time_index]
        this_mask = atlas_interpolated_data == label
        label_data = this_data[this_mask]
        label_time_courses[label_index, time_index] = np.mean(label_data)
        
## identify max voxel and its voxel coordinates and extract time course
max_coordinates = np.zeros((n_labels, 3))
max_time_courses = np.zeros((n_labels, n_times))
for label_index, label in enumerate(indices):
    this_mask = atlas_interpolated_data == label
    label_data = data[this_mask, :]
    argmax = np.argmax(np.max(label_data, axis=1))
    max_time_courses[label_index, :] = label_data[argmax, :]
    xs, ys, zs = np.where(this_mask)
    max_coordinates[label_index] = xs[argmax], ys[argmax], zs[argmax]
    
## find vertices
# FIXME
stc_voxels = np.array(
    np.unravel_index(stc.vertices, img.shape[:3], order='F')).T

vert_nos = np.zeros(n_labels)
n_xs, n_ys, n_zs = src[0]['shape']

for label_index, label in enumerate(labels):
    max_coordinate = max_coordinates[label_index]
    distances = np.zeros(len(stc_voxels))
    for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
        distances[stc_voxel_index] = np.linalg.norm(max_coordinate - stc_voxel)
    
    vert_nos[label_index] = stc.vertices[np.argmin(distances)]
        
## and map them to rois
roi_vertices = dict()
for label_index, label in enumerate(labels):
    if label in rois:
        roi_vertices[label] = int(vert_nos[label_index])

## plot label time courses (mean)
plt.close('all')
       
for label_index, label in enumerate(labels):
    if label in rois:
        plt.figure()
        plt.plot(times, label_time_courses[label_index, :].T)
        plt.title(labels[label_index])
        plt.ylim(-0.015, 0.015)
        plt.plot((times[0], times[-1]), (0, 0), 'k--')
        plt.xlim((times[0]), times[-1])
        
## plot maximally responding time courses
plt.close('all')
       
for label_index, label in enumerate(labels):
    if label in rois:
        plt.figure()
        plt.plot(times, max_time_courses[label_index, :].T)
        plt.title(labels[label_index])
        plt.ylim(-0.025, 0.025)
        plt.plot((times[0], times[-1]), (0, 0), 'k--')
        plt.xlim((times[0]), times[-1])

#%% get max vertices

scratch_path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/'
subjects_dir = join(scratch_path, 'freesurfer')
src = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                      'volume-7.5mm-src.fif'))

def get_max_vertices_per_roi(atlas, contrast):

    from nilearn import image, datasets
    from os import chdir
    from os.path import join
    import numpy as np
    import mne

    script_path = '/home/lau/projects/cerebellar_clock/scripts/tactile_jitter/' + \
                  'python/'
    chdir(script_path)              
    
    from qsubs import sensor_functions

    scratch_path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/'
    # subjects_dir = join(scratch_path, 'freesurfer')

    path = join(scratch_path, 'MEG', 'grand_averages', 'hilbert', 'stcs',
                'common_filter', 'contrasts')
    if contrast == 's1_vs_s2':
        filename = 'normal_ISI_no_baseline_Ss_mag_reg_0.0_' + \
            'weight_unit-noise-gain_s1_s2_contrast_7.5_mm_' + \
            'hp_4_Hz_lp_7_Hz_-400_ms_to_1100_ms_tactile_jitter_morph-vl.nii'
    elif contrast == 'o0_vs_o15':
        filename = 'normal_ISI_no_baseline_mag_reg_0.0_' + \
            'weight_unit-noise-gain_o0_o15_contrast_7.5_mm_' + \
            'hp_14_Hz_lp_30_Hz_-750_ms_to_750_ms_tactile_jitter_morph-vl.nii'
    img = image.load_img(join(path, filename))
    data = np.asanyarray(img.dataobj)


    ## source space and an stc
    

    
    ## any fsaverage stc will do - just need the vertices
    contrast = ['o0', 'o15']
    tmin, tmax = (-0.750, 0.750)
    freqs = [14, 30]
    ISI_type = 'normal'
    proj_string = ''#'proj_'
    baseline_string =  'no_baseline_'
    
    
    freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])
    
    # contrast_name = contrast[0] + '_vs_' + contrast[1]
    contrast_string = ''
    for part in contrast:
        contrast_string += part + '_'
    contrast_string += 'contrast_'    
        
    time_string = sensor_functions.time_string(tmin, tmax)
    data_path = join(scratch_path, 'MEG', 'grand_averages')
    hilbert_path = join('hilbert', 'stcs', 'common_filter', 'contrasts')
    filename = ISI_type + '_ISI_' + baseline_string + \
                'mag_reg_0.0_weight_unit-noise-gain_' + \
                contrast_string + '7.5_mm_' + freq_string + \
                time_string + proj_string + 'tactile_jitter_morph-vl.stc'
    stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))



## labels of interest
    if atlas == 'aal':
        rois = [
                'Cerebelum_6_L',
                'Cerebelum_6_R',
                'Cerebelum_8_L',
                'Cerebelum_8_R',
                'Thalamus_L',
                'Thalamus_R',
                'Putamen_L',
                'Putamen_R',
                'Precuneus_L',
                'Precuneus_R',
                'Insula_L',
                'Insula_R',
                'Postcentral_L',
                'Postcentral_R',
                'Parietal_Inf_L',
                'Parietal_Inf_R',
                'Cingulum_Mid_L',
                'Cingulum_Mid_R',
                ]
        
        ## aal atlas
        atlas_dataset = datasets.fetch_atlas_aal()
        atlas_filepath = atlas_dataset['maps']
        atlas_img = image.load_img(atlas_filepath)
        labels = atlas_dataset['labels']
        indices = list(map(int, atlas_dataset['indices']))
    elif atlas == 'harvard_oxford':
        # harvard oxford
        rois = [
          'Parietal Operculum Cortex',
          ]
        atlas_filename = 'cort-maxprob-thr0-1mm'
        atlas_dataset = datasets.fetch_atlas_harvard_oxford(atlas_filename)
        atlas_filepath = atlas_dataset.maps
        atlas_img = image.load_img(atlas_filepath)
        labels = atlas_dataset['labels']
        indices = list(map(int, np.unique(atlas_img.get_fdata())))
    elif atlas == 'cerebellar':

        # cerebellar atlas
        atlas_type = 'flirt'
        atlas_path = '/usr/local/fsl/data/atlases/'
        atlas_img = image.load_img(join(atlas_path, 'Cerebellum',
                                          'Cerebellum-MNI' + atlas_type + \
                                          '-maxprob-thr0-1mm.nii.gz'))
        label_path = join(atlas_path, 'Cerebellum_MNI' + atlas_type + '.xml')
        label_file = open(label_path, 'r')
        labels = []
        for line in label_file:
            if line[:6] != '<label':
                continue
            first_index = None
            last_index = None
            for character_index, character in enumerate(line):
                if character == '>' and first_index is None:
                    first_index = character_index + 1
                if character == '<' and first_index is not None:
                    last_index = character_index
            if first_index is not None and last_index is not None:
                labels.append(line[first_index:last_index])
        indices = range(1, len(labels))     
        

    ## resampling
    atlas_interpolated = image.resample_to_img(atlas_img, img, 'nearest')
    atlas_interpolated_data = np.asanyarray(atlas_interpolated.dataobj)
    
    ## find peak voxels
    
    n_labels = len(labels)

   
    if atlas == 'aal':     
        ## identify max voxel and its voxel coordinates and extract time course
        max_coordinates = np.zeros((n_labels, 3))
        for label_index, label in enumerate(indices):
            this_mask = atlas_interpolated_data == label
            label_data = data[this_mask, :]
            argmax = np.argmax(np.max(label_data, axis=1))
            xs, ys, zs = np.where(this_mask)
            max_coordinates[label_index] = xs[argmax], ys[argmax], zs[argmax]
        
        ## find vertices
        # 
        stc_voxels = np.array(
            np.unravel_index(stc.vertices, img.shape[:3], order='F')).T
        
        vert_nos = [None] * n_labels
        n_xs, n_ys, n_zs = src[0]['shape']
        
        for label_index, label in enumerate(labels):
            max_coordinate = max_coordinates[label_index]
            distances = np.zeros(len(stc_voxels))
            for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
                distances[stc_voxel_index] = \
                    np.linalg.norm(max_coordinate - stc_voxel)
            vert_nos[label_index] = stc.vertices[np.argmin(distances)]
        
    if atlas == 'harvard_oxford':
        vert_nos = [None] * n_labels


        stc_voxels = np.array(
            np.unravel_index(stc.vertices, img.shape[:3], order='F')).T
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
                closest_vertex = stc.vertices[np.argmin(distances)]
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
            vert_nos[label_index].append(stc.vertices[np.argmin(distances)])
            
            ## right
            argmax = np.argmax(np.max(right_data, axis=1))
            max_coordinate = xs[rights][argmax], ys[rights][argmax], \
                            zs[rights][argmax]
            distances = np.zeros(len(stc_voxels))
            for stc_voxel_index, stc_voxel in enumerate(stc_voxels):
                distances[stc_voxel_index] = \
                    np.linalg.norm(max_coordinate - stc_voxel)
            vert_nos[label_index].append(stc.vertices[np.argmin(distances)])
                    
                
    ## and map them to rois
    roi_vertices = dict()
    for label_index, label in enumerate(labels):
        if label in rois:
            if atlas == 'aal':
                roi_vertices[label] = vert_nos[label_index]
            if atlas == 'harvard_oxford':
                roi_vertices[label + '_L'] = vert_nos[label_index][0]
                roi_vertices[label + '_R'] = vert_nos[label_index][1]
            
    ## fix names of rois
    for roi in roi_vertices:
        if 'Cerebelum' in roi:
            new_roi = roi.replace('Cerebelum', 'Cerebellum')
            roi_vertices[new_roi] = roi_vertices.pop(roi)
        if 'Postcentral' in roi:
            new_roi = roi.replace('Postcentral', 'SI')
            roi_vertices[new_roi] = roi_vertices.pop(roi)
        if 'Parietal Operculum Cortex' in roi:
            new_roi = roi.replace('Parietal Operculum Cortex', 'SII')
            roi_vertices[new_roi] = roi_vertices.pop(roi)
            
            
    return roi_vertices


#%% PLAY WITH ENVELOPE CORRELATION

import mne
from os.path import join
import numpy as np
import time
import matplotlib.pyplot as plt


project_path = '/home/lau/projects/cerebellar_clock'
raw_path = join(project_path, 'raw')
scratch_path = join(project_path, 'scratch', 'tactile_jitter', 'MEG')
subjects_dir = join(project_path, 'scratch', 'tactile_jitter', 'freesurfer')
figures_path = join(project_path, 'scratch', 'tactile_jitter', 'figures')

subject = '0013'
date = '20200124_000000'
subject_path = join(scratch_path, subject, date, 'hilbert')
filename = 'normal_no_baseline_hp_14_Hz_lp_30_Hz_' + \
            '-750_ms_to_750_ms_tactile_jitter_hilbert-epo.fif'
            
epochs_hilbert = mne.read_epochs(join(subject_path, filename), proj=False,
                         preload=False)            
omissions = mne.epochs.concatenate_epochs([epochs_hilbert['o0'],
                                           epochs_hilbert['o15']])
omissions.pick_types('mag')
omissions.del_proj()


bem_folder = join(subjects_dir, subject, 'bem')

src = mne.read_source_spaces(join(bem_folder, 'volume-7.5mm-src.fif'))
trans = mne.read_trans(join(bem_folder, subject + '-trans.fif'))
bem = mne.read_bem_solution(join(bem_folder, subject + '-5120-bem-sol.fif'))
fwd = mne.make_forward_solution(omissions.info, trans, src, bem, n_jobs=-1)

data_cov = mne.compute_covariance(omissions)
filters = mne.beamformer.make_lcmv(omissions.info, fwd, data_cov, 0.00,
                                   pick_ori='max-power',
                                   weight_norm='unit-noise-gain',
                                   noise_cov=None)

## keep 10 epochs from each for speed
drop_indices = np.concatenate((np.arange(10, 150), np.arange(160, 300)))
omissions.drop(drop_indices)

stcs_list = []
events = ['o0', 'o15']
for event in events:
    stcs = mne.beamformer.apply_lcmv_epochs(omissions[event],
                                                      filters,
                                                      return_generator=False)
    stcs_list.append(stcs)

    
    
morph = mne.read_source_morph(join(subject_path, '..',
                                    '7.5_mm_tactile_jitter_vol-morph.h5'))

stcs_morph_list = []
counter = 0
for stcs in stcs_list:
    stcs_morph = []
    for stc in stcs:
        stc.crop(-0.400, 0.400)
        print('Morphing stc ' + str(counter))
        begin_time = time.time()
        stcs_morph.append(morph.apply(stc))
        print(str(int(time.time() - begin_time)) + ' s elapsed')
        counter += 1
    stcs_morph_list.append(stcs_morph)

# aal_vertices = get_max_vertices_per_roi('aal', 'o0_vs_o15')    
# SII_vertices = get_max_vertices_per_roi('harvard_oxford', 'o0_vs_o15')



vertices = {
     # 'Insula_L': 7229,
     # 'Insula_R': 6595,
     # 'Cingulum_Mid_L': 10222,
     # 'Cingulum_Mid_R': 9557,
     'Parietal_Inf_L': 9527,
     'Parietal_Inf_R': 9538,
     # 'Precuneus_L': 7048,
     # 'Precuneus_R': 7671,
     'Putamen_L': 7230,
     'Putamen_R': 5972,
     # 'Thalamus_L': 6564,
     # 'Thalamus_R': 7809,
     'Cerebelum_8_L': 2654,
     'Cerebelum_8_R': 2706,
     'SI_L': 10285,
     'SI_R': 10229,
     # 'Cerebellum_6_L': 3365,
     # 'Cerebellum_6_R': 3947,
     'SII_L': 8952,
     'SII_R': 8367
     }

n_epochs = len(stcs_morph_list[0])
n_vertices = len(vertices)


stc_morph = stcs_morph_list[0][0]
data_indices = np.zeros(n_vertices, dtype=int)
for vertex_index, vertex in enumerate(vertices):
    data_indices[vertex_index] = \
        np.where(stc_morph.vertices[0] == vertices[vertex])[0]


## reduce data
corr_data_list = []
for stcs_morph in stcs_morph_list:
    corr_data = np.zeros((n_epochs, n_vertices, len(stcs_morph[0].times)),
                         dtype=complex)
    for stc_morph_index, stc_morph in enumerate(stcs_morph):
        corr_data[stc_morph_index, :, :] = stc_morph.data[data_indices, :]
    corr_data_list.append(corr_data)
        
corr_list = []
for event_index in range(len(corr_data_list)):  
    print('Running event: ' + str(event_index))
    begin_time = time.time()
    corr = mne.connectivity.envelope_correlation(corr_data_list[event_index])
    print(str(int(time.time() - begin_time)) + ' s elapsed')
    corr_list.append(corr)



src_fsaverage = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                            'volume-7.5mm-src.fif'))    

n_sources = stcs_morph_list[0][0].shape[0]
## degrees
plt.close('all')
stc_degrees = []
for corr in corr_list:
    degree = mne.connectivity.degree(corr, 0.5)
    dipoles = np.zeros(n_sources)
    dipoles[data_indices] = degree
    stc = mne.VolSourceEstimate(dipoles,
                                stcs_morph_list[0][0].vertices, 0, 1,
                                'fsaverage')
    stc.plot(src_fsaverage, clim=dict(kind='value', lims=[1,
                                                          np.median(degree),
                                                          max(degree)]),
             mode='glass_brain')
    
mne.viz.plot_connectivity_circle    
#%%# load evokeds_dict

import numpy as np
from scipy import stats as stats

n_subjects = 30
p_threshold = 0.05
t_threshold = -stats.distributions.t.ppf(p_threshold / 2.0,
                                         n_subjects - 1)

array = np.zeros((30, 306, 1501))
hilb_array = np.zeros((30, 312, 1501))
for subject_index in range(n_subjects):
    array[subject_index, :, :] = evokeds_dict['o0_vs_o15'][subject_index].data
                                  # evokeds_dict['o15'][subject_index].data
    # array[subject_index, :, :] = evokeds_dict['o0'][subject_index].data - \
    #                               evokeds_dict['o15'][subject_index].data
                                 
    # hilb_array[subject_index, :, :] = np.abs(evokeds_dict['o0'][subject_index].apply_hilbert().data) - \
    #                               np.abs(evokeds_dict['o15'][subject_index].apply_hilbert().data)

temp = evokeds_dict['o0_vs_o15'][0].copy()                                 
temp = evokeds_dict['o0'][0].copy()
hilb_temp = evokeds_dict['o0'][0].copy()
temp._data = array
hilb_temp._data = hilb_array
picks = mne.pick_types(temp.info, 'mag')

stat_array = temp.data[:, picks, 201:601]


# picks = mne.channels.layout._pair_grad_sensors(temp.info, topomap_coords=False)
# data = np.zeros((30, 102, 1501))
# for subject_index in range(n_subjects):
#     data[subject_index, :, :] = \
#         mne.channels.layout._merge_grad_data(array[subject_index, picks, :])
        
# temp.pick_types('mag')
# temp._data = data
# for subject_index in range(n_subjects):
#     for channel_index in range(102):
#         temp._data[subject_index, channel_index, :] -= \
#             np.mean(temp.data[subject_index, channel_index, :])

t = np.mean(temp.data, 0) / (np.std(temp.data, 0) / np.sqrt(n_subjects))
hilb_t = np.mean(hilb_temp.data, 0) / (np.std(hilb_temp.data, 0) / np.sqrt(n_subjects))
mean = np.mean(temp.data, 0)
hilb_mean = np.mean(hilb_temp.data, 0)
t[np.logical_and(t < t_threshold, t > -t_threshold)] = 0
hilb_t[np.logical_and(hilb_t < t_threshold, hilb_t > -t_threshold)] = 0
temp._data = t
hilb_temp._data = hilb_t
temp._data = mean
hilb_temp._data = hilb_mean

temp.plot_topo(scalings=dict(mag=1, grad=1))
hilb_temp.plot_topo(scalings=dict(mag=1, grad=1))



t_obs, clusters, cluster_p_values, H0 = clu = \
    mne.stats.permutation_cluster_1samp_test(stat_array,
                                               n_permutations=1024, n_jobs=8)
temp = evokeds_dict['o0_vs_o15'][0].copy()                                 
cl = temp.copy()
cl._data = cl.data[:, 201:601]
temp._data = array
temp.pick_types('mag')
temp._data = temp.data[:, :, 201:601]


cl.pick_types('mag')

mean = np.mean(temp.data, 0)
# cl.crop(tmin, tmax)
old_data = mean
cl._data = np.zeros(cl.data.shape)
min_cluster_index = np.argmin(cluster_p_values)
lala, cluster_indices = np.unique(cluster_p_values, return_index=True)
min_cluster_index = cluster_indices[2]
cluster = clusters[min_cluster_index]
cl._data[cluster] = old_data[cluster]
cl._times = cl.times[201:601]

cl.plot_topo(scalings=dict(mag=1))





## plot hilbert ratios

import mne
from os.path import join
import matplotlib.pyplot as plt


plt.close('all')

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

path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/MEG'
n_subjects = len(subjects)


comparisons = []

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index] + '_000000'

    subject_path = join(path, subject, date, 'hilbert')
    filename = 'hp_14_Hz_lp_30_Hz_tactile_jitter_hilbert-ave.fif'
    
    evokeds = mne.read_evokeds(join(subject_path, filename), ['o0', 'o15'])
    comparison = evokeds[0].copy()
    comparison._data = (evokeds[0].data - evokeds[1].data) / \
                       (evokeds[0].data + evokeds[1].data)
    comparison.plot(picks='mag')
    comparisons.append(comparison)
    
ga = mne.grand_average(comparisons)

## evetns

import mne

path = '/home/lau/projects/cerebellar_clock/' + \
    'raw/0001/20200124_000000/MEG/001.tactile_jitter_raw/files/' + \
        'tactile_jitter_raw.fif'
        
raw = mne.io.read_raw_fif(path)        
events = mne.find_events(raw, min_duration=0.002)
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
        mean_ISI = int(np.mean([events[s4_index, 0] - events[s3_index, 0],
                                events[s5_index, 0] - events[s4_index, 0],
                                events[s6_index, 0] - events[s5_index, 0]]))
        # print('Local ISI: ' + str(local_ISI) + ' ms')
        # print('Mean ISI: ' + str(mean_ISI) + ' ms')
        local_sample = events[s6_index, 0] + local_ISI
        mean_sample  = events[s6_index, 0] + mean_ISI
        
        local_events[index, 0] = local_sample
        mean_events[index, 0] = mean_sample
        
## for plotting
events_for_plotting = np.concatenate((events[events[:, 2] == 28],
            local_events[local_events[:, 2] == 28] + np.array([0, 0, 1]),
            mean_events[mean_events[:, 2] == 28] + np.array([0, 0, 2]),
            events[events[:, 2] == 38],
            local_events[local_events[:, 2] == 38] + np.array([0, 0, 1]),
            mean_events[mean_events[:, 2] == 38] + np.array([0, 0, 2])))

mne.viz.plot_events(events_for_plotting)
        
        


def compare_with_and_without_nai(subject, date, spacing, contrast, h_freq,
                                 l_freq, org_name_input, org_name_output,
                                 reg, weight_norm):

    from os.path import join
    from os import listdir, chdir
    import numpy as np
    import matplotlib.pyplot as plt
    import mne
    from sensor_functions import get_load_and_save_paths_for_common_filter_hilbert, get_save_path 
    
    paths = get_load_and_save_paths_for_common_filter_hilbert(subject, date,
                                                      spacing, contrast,
                                                      h_freq, l_freq,
                                                      org_name_input,
                                                      org_name_output,
                                                      reg, weight_norm,
                                                      grand_average=False,
                                                      statistics=False,
                                                      p_threshold=False)
    
    path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter'

    subjects_dir = join(path, 'freesurfer')
    

    
    src = mne.source_space.read_source_spaces(join(subjects_dir, subject,
                                                       'bem',
                                                       'volume-7.5mm-src.fif'))
        
    print(paths[0][4])
    stc = mne.read_source_estimate(paths[0][4])
    stc.plot(src, src[0]['subject_his_id'])        





## plot t-values and cluster

          
# look at clusters after getting t_obs

# stc, tmin, tmax and cluster, p_threshold should be present

import numpy as np
from scipy import stats as stats
mne.viz.set_3d_backend('pyvista')
n_subjects = 30

t_threshold = -stats.distributions.t.ppf(p_threshold / 2.0,
                                         n_subjects - 1)

t_obs = clu['t_obs']
clusters = clu['clusters']
cluster_p_values = clu['cluster_p_values']


## tfr

# tfr = tfrs[0].copy()
# tfr._data = np.mean(array, axis=0)
# tfr.pick_types('mag')
# tfr._data = t_obs

# tfr._data[tfr.data < t_threshold] = 0

# tfr.plot_topo()

## stc
# t = stc.copy()
# t.crop(tmin, tmax)
# t._data = t_obs.T

# t.plot(src, src[0]['subject_his_id'])
# t.plot('fsaverage', hemi='both')




cl = stc.copy()
cl.crop(tmin, tmax)
old_data = cl.data
cl._data = np.zeros(cl.data.shape)
min_cluster_index = np.argmin(cluster_p_values)
cluster = clusters[min_cluster_index]
cl._data[cluster[1], cluster[0]] = old_data[cluster[1], cluster[0]]
# cl._data[cluster[1], cluster[0]] = t_obs[cluster[0], cluster[1]]

cl_size = np.sum(t_obs[cluster[0], cluster[1]])
# clim = dict(kind='value', lims=(0, t_threshold, np.quantile(cl.data, 0.99)))
# cl.plot(src, src[0]['subject_his_id'], clim=clim, initial_time=0,
        # mode='stat_map')
# clim = dict(kind='value', lims=(0, np.quantile(cl.data, 0.95),
#                                 np.quantile(cl.data, 0.99)))
# clim = dict(kind='value', lims=(0.0,
#                                 np.quantile(cl.data, 0.975),
#                                 np.quantile(cl.data, 0.99)))
# cl.plot(src, src[0]['subject_his_id'], clim=clim, initial_time=0)
cl.plot(src, src[0]['subject_his_id'], initial_time=0)

# cl.plot('fsaverage', hemi='both', clim=clim)


mne.source_space.setup_volume_source_space
mne.stats.summarize_clusters_stc
## 
import mne
import numpy as np
from os import chdir
from os.path import join

mne.viz.set_3d_backend('pyvista')

subject_path = '/home/lau/mounts/hyades/scratch5' + \
                '/MINDLAB2019_MEG-CerebellarClock/tactile_jitter/MEG' + \
                '/grand_averages/'
               
chdir(subject_path)

stat = np.load(join(subject_path, 'statistics', 'hilbert', 'itcs',
                    'common_filter', 'contrasts',
        'p_0.05_o5_ns_contrast_7.5_mm_hp_8_Hz_lp_12_Hz_tactile_jitter_morph-vl.npy'),
                allow_pickle=True).item()
stc = mne.read_source_estimate(join(subject_path, 'hilbert', 'itcs',
                                    'common_filter', 'contrasts',
    'o5_ns_contrast_7.5_mm_hp_8_Hz_lp_12_Hz_tactile_jitter_morph-vl.stc'))
src = \
  mne.source_space.read_source_spaces('/home/lau/projects/cerebellar_clock/' + \
          'scratch/tactile_jitter/freesurfer/fsaverage/bem/volume-7.5mm-src.fif')


#%% TEST LCMV

import mne
from os import chdir
from os.path import join
import numpy as np

def collapse_conditions(epochs):
    to_collapse = {}
    to_collapse['ns']    = ['n1_0',  'n2_0',  'n3_0',  'n4_0',  'n5_0',
                            'n1_5',  'n2_5',  'n3_5',  'n4_5',  'n5_5',
                            'n1_15', 'n2_15', 'n3_15', 'n4_15', 'n5_15'
                            ]
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


path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/MEG/' + \
       '0001/20200124_000000/hilbert'
chdir(path)

raw = mne.io.read_raw_fif('hp_14_Hz_lp_30_Hz_tactile_jitter_raw.fif')
epochs_hilbert = \
    mne.read_epochs('hp_14_Hz_lp_30_Hz_tactile_jitter_hilbert-epo.fif',
                    proj=False, preload=False)
epochs_hilbert = collapse_conditions(epochs_hilbert)

subject = '0001'
subjects_dir = path[:59] + 'freesurfer'
bem_folder = join(subjects_dir, subject, 'bem')

ico_string = '5120'
spacing = 7.5
trans = join(bem_folder, subject + '-trans.fif')
src = join(bem_folder, 'volume-' + str(spacing) + 'mm-src.fif')
bem = join(bem_folder, subject + '-' + ico_string + '-bem-sol.fif')

fwd = mne.make_forward_solution(epochs_hilbert.info, trans, src, bem)
src = mne.read_source_spaces(src) 
events = epochs_hilbert.events
tmin = epochs_hilbert.tmin
tmax = epochs_hilbert.tmax
baseline = epochs_hilbert.baseline
decim = 1
reject =  epochs_hilbert.reject
event_id = 18

epochs_cov = mne.Epochs(raw, events, event_id, tmin, tmax,
                        baseline, proj=False, decim=decim,
                        reject=reject)
epochs_cov.del_proj()
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
epochs_hilbert['o0'].load_data()
stcs = \
    mne.beamformer.apply_lcmv_epochs(epochs_hilbert['o0'],
                                                     filters,
                                                     max_ori_out='signed')
    
for stc in stcs:
    stc._data = np.array(np.abs(stc.data),
                 dtype='float64')    
    
stc_mean = stcs[0].copy()
mean_data = np.mean([stc.data for stc in stcs], axis=0)
stc_mean._data = mean_data
stc_mean.plot(src, src[0]['subject_his_id'])


#%% MAKE VOLUME SPACE WITH LABELS

from os.path import join
import mne
path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter'
meg_path = join(path, 'MEG')
subjects_dir = join(path, 'freesurfer')
subject = '0001'
bem_dir = join(subjects_dir, subject, 'bem')
fname_aseg = join(subjects_dir, subject, 'mri', 'aparc+aseg.mgz')
fname_model = join(bem_dir, '%s-5120-bem.fif' % subject)

labels_vol = [
            'Left-Cerebral-White-Matter',
                #     'Left-Amygdala',
              # 'Left-Thalamus-Proper',
              # 'Left-Cerebellum-Cortex',
              # 'Brain-Stem',
              # 'Right-Amygdala',
              # 'Right-Thalamus-Proper',
              # 'Right-Cerebellum-Cortex',
              'ctx-rh-precuneus'
              ]

src = mne.setup_volume_source_space(subject, 7.5, bem=fname_model)

labels_aseg = mne.get_volume_labels_from_src(src, subject, subjects_dir)

#%% TEST FOR ZERO VOXELS

import mne
from os.path import join
from os import chdir


path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter'
meg_path = join(path, 'MEG')
subjects = [
    '0001', '0002', '0004', '0005', '0006', '0007',
    '0008', '0009', '0010', '0011', '0012', '0013',
    '0014', '0015', '0016', '0017', '0018', '0019',
    '0020', '0021', '0022', '0023', '0024', '0025',
    '0026', #'0027',
    '0028', '0029', '0030', '0031',
    ]
dates = [
    '20200124', '20200121', '20200121', '20200128', '20200120', '20200120',
    '20200120', '20200121', '20200122', '20200122', '20200122', '20200124',
    '20200127', '20200128', '20200128', '20200129', '20200129', '20200129',
    '20200131', '20200131', '20200131', '20200203', '20200203', '20200203',
    '20200204', #'20200204',
    '20200204', '20200205', '20200205', '20200205'
         ]


for subject_index, subject in enumerate(subjects):
    date = dates[subject_index] + '_000000'
    chdir(join(meg_path, subject, date, 'hilbert', 'stcs'))
    print('Loading subject: ' + subject)
    stc_morph = mne.read_source_estimate('o0_hp_8_Hz_lp_12_Hz_tactile_jitter_morph-vl.stc')
    print('Number of voxels equal to 0: ' + str(np.sum(stc_morph.data == 0) // 1501))

    

#%% Plot hilbert beamformer subject


def plot_hilbert(subject, freq='all'):

    from os.path import join
    from os import listdir, chdir
    import numpy as np
    import matplotlib.pyplot as plt
    import mne
    
    path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter'
    meg_path = join(path, 'MEG')
    subjects_dir = join(path, 'freesurfer')
    
    if freq == 'all':
        h_freqs = ['4', '8', '14', '60']
        l_freqs = ['7', '12', '30', '85']
    
    if subject != 'fsaverage':
        dates = listdir(join(meg_path, subject))
        for date in dates:
            if '2020' in date:
                break
            
    if subject == 'fsaverage':
        morph_text = '_morph'
        chdir(join(meg_path, 'grand_averages', 'hilbert', 'stcs'))
    else:
        chdir(join(meg_path, subject, date, 'hilbert', 'stcs'))
        morph_text = ''

    conditions = ['o0', 'o15']
    src = mne.source_space.read_source_spaces(join(subjects_dir, subject,
                                                       'bem',
                                                       'volume-7.5mm-src.fif'))
    plt.close('all')
    # conditions = ['o0', 'o5', 'o15']
    # conditions = ['s1', 's2']
    for (h_freq, l_freq) in zip(h_freqs, l_freqs):
        stcs = []
        ratios = []
        diffs =  []
        for condition in conditions:
            stc = mne.read_source_estimate(condition + '_7.5_mm_' +  \
                                   'hp_' + h_freq + '_Hz_lp_' + l_freq + \
                                '_Hz_tactile_jitter' + \
                                   morph_text + '-vl.stc')
            stcs.append(stc)
        

        # contrasts = [[0, 1], [0, 2], [1, 2]]
        contrasts = [[0, 1]]
        for contrast in contrasts:
            ratio = stcs[contrast[0]].copy()
            copy_1 = stcs[contrast[0]].copy()
            copy_2 = stcs[contrast[1]].copy()
            copy_1._data = (copy_1.data - np.mean(copy_1.data)) / \
                                np.std(copy_1.data)
            copy_2._data = (copy_2.data - np.mean(copy_2.data)) / \
                                np.std(copy_2.data)                                
            
            # ratio._data = (stcs[contrast[0]].data - stcs[contrast[1]].data) / \
            #               (stcs[contrast[0]].data + stcs[contrast[1]].data)
            ratio._data = (copy_1.data - copy_2.data) / \
                          (copy_1.data + copy_2.data)
            ratios.append(ratio)
            
            diff = stcs[contrast[0]].copy()
            diff._data = stcs[contrast[0]].data - stcs[contrast[1]].data
            diffs.append(diff)
            
        
        for stc in stcs:
            max_voxel, max_time = np.unravel_index(stc.data.argmax(),
                                                   stc.data.shape)
            max_voxels = stc.data[:, max_time]
            lims = (np.quantile(max_voxels, 0.95), np.quantile(max_voxels, 
                                                             0.975),
                    np.quantile(max_voxels, 0.99))
            stc.plot(src, src[0]['subject_his_id'], clim=dict(kind='value',
                                                              lims=lims))
            
        for ratio in ratios:
            lims = (np.quantile(ratio.data, 0.95), np.quantile(ratio.data, 
                                                             0.975),
                    np.quantile(ratio.data, 0.99))
            ratio.plot(src, src[0]['subject_his_id'],# clim=dict(kind='value',
                                                         #       lims=lims),
                       initial_time=-0.028)
        
    # for diff in diffs:
    #     diff.plot(src, src[0]['subject_his_id'])
        

#%% Test beamformer contrast on subject 1

import mne
import numpy as np
from os import chdir

subject_path = '/home/lau/mounts/hyades/scratch5' + \
               '/MINDLAB2019_MEG-CerebellarClock/tactile_jitter/MEG' + \
               '/0001/20200124_000000'
               
chdir(subject_path)               

## epochs

epochs = mne.read_epochs('tactile_jitter_tfr_sss-epo.fif',
                         preload=False, proj=False)

comb = epochs[['n1_0', 'n1_5', 'n1_15',
             'n2_0', 'n2_5', 'n2_15',
             'n3_0', 'n3_5', 'n3_15',
             'n4_0', 'n4_5', 'n4_15',
             'n5_0', 'n5_5', 'n5_15',
             'o0']]
o0 = epochs['o0']
ns = epochs[['n1_0', 'n1_5', 'n1_15',
             'n2_0', 'n2_5', 'n2_15',
             'n3_0', 'n3_5', 'n3_15',
             'n4_0', 'n4_5', 'n4_15',
             'n5_0', 'n5_5', 'n5_15']]

comb.load_data()
o0.load_data()
ns.load_data()
comb.pick_types(meg='grad')
o0.pick_types(meg='grad')
ns.pick_types(meg='grad')

## csds

freqs = np.arange(10, 15, 1)

csd = mne.time_frequency.csd_morlet(comb, freqs, tmin=-0.400, tmax=1.000,
                                    decim=1, n_jobs=1)
csd_o0 = mne.time_frequency.csd_morlet(o0, freqs, tmin=-0.400, tmax=1.000,
                                       decim=1)
csd_ns = mne.time_frequency.csd_morlet(ns, freqs, tmin=-0.400, tmax=1.000,
                                       decim=1)

## forward

subject = '0001'
bem = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-5120-bem.fif'
        
src = mne.source_space.setup_volume_source_space(subject, bem=bem)  
filename = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/volume-5mm-src.fif'
# mne.source_space.write_source_spaces(filename, src)        

trans = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-trans.fif'
 
bem_sol = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-5120-bem-sol.fif'    
fwd = mne.make_forward_solution(o0.info, trans, src, bem_sol)        

## DICS

filters = mne.beamformer.make_dics(epochs.info, fwd, csd.mean(),
                                   pick_ori='max-power')

o0_source_power, freqs = mne.beamformer.apply_dics_csd(csd_o0.mean(), filters)
ns_source_power, freqs = mne.beamformer.apply_dics_csd(csd_ns.mean(), filters)

stc = o0_source_power / ns_source_power

lims = [0.025, 0.05, 0.10]
stc.plot(src=fwd['src'], subject=subject, clim=dict(kind='value',
                                                    pos_lims=lims))

stc.plot(src=fwd['src'], subject=subject, clim=dict(kind='value',
                                                    lims=lims),
                     mode='glass_brain')

#%% test beamformer  on subject 1

import mne
import numpy as np
from os import chdir

subject_path = '/home/lau/mounts/hyades/scratch5' + \
               '/MINDLAB2019_MEG-CerebellarClock/tactile_jitter/MEG' + \
               '/0001/20200124_000000'
               
chdir(subject_path)

epochs = mne.read_epochs('tactile_jitter_tfr_sss-epo.fif', preload=False)
s1 = epochs['s1']
s1.load_data()
s1.pick_types(meg='grad')

freqs = np.arange(14, 19, 1)

csd = mne.time_frequency.csd_morlet(s1, freqs, tmin=0.500, tmax=1.000,
                                    decim=1, n_jobs=1)

## forward

subject = '0001'
bem = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-5120-bem.fif'
        
# src = mne.source_space.setup_volume_source_space(subject, bem=bem)  
filename = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/volume-5mm-src.fif'
src = mne.source_space.read_source_spaces(filename)        

trans = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-trans.fif'
 
bem_sol = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/' + \
        'freesurfer/0001/bem/0001-5120-bem-sol.fif'    
fwd = mne.make_forward_solution(s1.info, trans, src, bem_sol, n_jobs=-1)

filters = mne.beamformer.make_dics(s1.info, fwd, csd.mean(),
                                   pick_ori=None, real_filter=False,
                                   weight_norm='unit-noise-gain',
                                   normalize_fwd=False, reg=0.20)

stc, freqs = mne.beamformer.apply_dics_csd(csd.mean(), filters)

lims = [np.quantile(stc.data, 0.90), np.quantile(stc.data, 0.95),
        np.quantile(stc.data, 0.975)]

stc.plot(src, src[0]['subject_his_id'], clim=dict(kind='value', lims=lims,
                                                  mode='glass_brain'))


#%% plot stat

import mne
import numpy as np
from os import chdir
from os.path import join

mne.viz.set_3d_backend('pyvista')

subject_path = '/home/lau/mounts/hyades/scratch5' + \
               '/MINDLAB2019_MEG-CerebellarClock/tactile_jitter/MEG' + \
               '/grand_averages/'
               
chdir(subject_path)

stat = np.load(join(subject_path, 'statistics', 'hilbert', 'itcs',
                    'common_filter', 'contrasts',
        'p_0.05_o5_ns_contrast_7.5_mm_hp_8_Hz_lp_12_Hz_tactile_jitter_morph-vl.npy'),
               allow_pickle=True).item()
stc = mne.read_source_estimate(join(subject_path, 'hilbert', 'itcs',
                                    'common_filter', 'contrasts',
    'o5_ns_contrast_7.5_mm_hp_8_Hz_lp_12_Hz_tactile_jitter_morph-vl.stc'))
src = \
 mne.source_space.read_source_spaces('/home/lau/projects/cerebellar_clock/' + \
          'scratch/tactile_jitter/freesurfer/fsaverage/bem/volume-7.5mm-src.fif')
     
cluster_indices = np.where(stat['cluster_p_values'] < 0.24)[0]

all_clusters = None

for cluster_index in cluster_indices[:1]:
    cluster = stat['clusters'][cluster_index]
    if all_clusters is None:
        all_clusters = [cluster]
    else:
        all_clusters = np.append(all_clusters, [cluster])


stat_stc = stc.copy()
stat_stc._data = np.zeros(stat_stc.shape)
for cluster in all_clusters:
    cluster_times = cluster[0]
    cluster_vertices = cluster[1]
    
    for (time, vertex) in zip(cluster_times, cluster_vertices):
        stat_stc._data[vertex, time] = stc.data[vertex, time]
          
lims = [0.0001, 0.001, 0.01]
stat_stc.plot(src, src[0]['subject_his_id'], clim=dict(kind='value',
               pos_lims=lims))
    
#%% play with hilberts

from os import chdir
from os.path import join
import mne
from mne.datasets import somato

data_path = somato.data_path()

load_path = join(data_path, 'sub-01', 'meg')
save_path = join(data_path, 'derivatives', 'sub-01')
chdir(load_path)

raw = mne.io.read_raw_fif('sub-01_task-somato_meg.fif', preload=True)
events = mne.find_events(raw)
raw_hilbert = raw.copy()
raw_hilbert.filter(8, 12)
raw_hilbert.apply_hilbert()

event_id, tmin, tmax = 1, -.2, 2     
baseline = (-.2, -0.05)

epochs_hilbert = mne.Epochs(raw_hilbert, events, event_id, tmin, tmax,
                            baseline)

print('Original has right dtype')
print(epochs_hilbert.get_data().dtype)

## now checking when loading them again

single_name = join(save_path, 'hilbert_single-epo.fif')
double_name = join(save_path, 'hilbert_double-epo.fif')

epochs_hilbert.save(single_name, fmt='single')
epochs_hilbert.save(double_name, fmt='double')

epochs_hilbert_read_single = mne.read_epochs(single_name, preload=True)
epochs_hilbert_read_double = mne.read_epochs(double_name, preload=True)

print('With preload=True')
print(epochs_hilbert_read_single.get_data().dtype)
print(epochs_hilbert_read_double.get_data().dtype)

epochs_hilbert_read_single = mne.read_epochs(single_name, preload=False)
epochs_hilbert_read_double = mne.read_epochs(double_name, preload=False)

print('With preload=False')
print(epochs_hilbert_read_single.get_data().dtype)
print(epochs_hilbert_read_double.get_data().dtype)