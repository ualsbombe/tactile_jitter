#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:46:34 2020

@author: lau
"""
#%% IMPORTS

import mne
from os.path import join
from os import chdir
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as stats
import warnings

script_path = '/home/lau/projects/cerebellar_clock/scripts/tactile_jitter/' + \
              'python/'
chdir(script_path)              

from qsubs import sensor_functions

scratch_path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter/'
subjects_dir = join(scratch_path, 'freesurfer')


#%% FUNCTIONS
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
            print(roi_vertices.keys())
        if 'Postcentral' in roi:
            new_roi = roi.replace('Postcentral', 'SI')
            roi_vertices[new_roi] = roi_vertices.pop(roi)
        if 'Parietal Operculum Cortex' in roi:
            new_roi = roi.replace('Parietal Operculum Cortex', 'SII')
            roi_vertices[new_roi] = roi_vertices.pop(roi)
            
            
    return roi_vertices




#%% EVOKEDS FIGURE

plt.close('all')

save_fig = True
data_path = join(scratch_path, 'MEG', 'grand_averages')
filename = 'lp_40_Hz_tactile_jitter_sss-ave.fif'
gas = mne.read_evokeds(join(data_path, filename))

mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'

ts_args = dict(ylim=dict(grad=[-30, 30]), spatial_colors=False, time_unit='ms')
topomap_args = dict(vmin=0, vmax=10, time_unit='ms')

events = dict(s1='First Stimulation', s2='Repeated Stimulation (2)',
              # s3='Repeated Stimulation (3)',
              o0='Omission (0)', ns='Non-Stimulation')

figure_path = join(scratch_path, 'figures', 'grand_averages', 'evokeds')

for ga in gas:
    if ga.comment in events.keys():
        title = events[ga.comment]
        fig = mne.viz.plot_evoked_joint(ga, times=(0.045, 0.120), picks='grad',
                                        ts_args=ts_args,
                                        topomap_args=topomap_args,
                                        title=title)
        fig.axes[3].text(-3, 11, 'fT/cm', {'size':10})
        fig_name = title + '.png'
        if save_fig:
            fig.savefig(join(figure_path, fig_name))

#%# STC GA

SI_vertex  = 144752
SII_vertex = 131610
data_path = join(scratch_path, 'MEG', 'grand_averages', 'stcs')
src = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                  'fsaverage-ico-5-src.fif'))
# first get Repstim2 time courses

event = 's2'
filename = event + '_lp_40_Hz_tactile_jitter_sss_morph_ico-5'

stc_2 = mne.read_source_estimate(join(data_path, filename))
SI_stc_vertex = src[0]['nearest'][SI_vertex]
SII_stc_vertex = src[0]['nearest'][SII_vertex]

# then get first stim

event = 's1'
filename = event + '_lp_40_Hz_tactile_jitter_sss_morph_ico-5'

stc = mne.read_source_estimate(join(data_path, filename))


clim = dict(kind='value', lims=(4.88, 5.46, 7.11))
mne.viz.set_3d_backend('pyvista')
# mne.viz.set_3d_backend('mayavi')

mpl.rcParams['font.size'] = 10
mpl.rcParams['font.weight'] = 'bold'

times = [45.0, 120.0]
names = ('First Stimulation', 'Repeated Stimulation (2)')
colours = ('b', 'r')

for time in times:
    
    brain = stc.plot('fsaverage', hemi='both', subjects_dir=subjects_dir,
                     initial_time=time, time_unit='ms', clim=clim,
                     show_traces=True, background='white', colorbar=True,
                     foreground='black')
    
    brain.add_text(0.35, 0.9, 'First Stimulation (' + str(int(time)) + ' ms)',
                   color=(0, 0, 0))

    axis = brain.time_viewer.mpl_canvas.axes
    brain.time_viewer.clear_points() ## remove default points
    brain.add_foci(SI_vertex, True, hemi='lh', color='green')
    brain.add_foci(SII_vertex, True, hemi='lh', color='green')
    brain.time_viewer.plot_time_course('lh', SI_vertex, 'blue')
    axis.plot(stc_2.times * 1e3, stc_2.data[SI_stc_vertex, :].T, 'r')
    brain.time_viewer.plot_time_course('lh', SII_vertex, 'blue')
    axis.get_lines()[-1].set_linestyle('--')
    brain.time_viewer.toggle_interface()
    ylim = axis.get_ylim()
    axis.plot(stc_2.times * 1e3, stc_2.data[SII_stc_vertex, :].T, 'r--')
    axis.plot((times[0], times[0]), ylim, 'k--')
    axis.plot((times[1], times[1]), ylim, 'k--')
    axis.legend_.remove()
    axis.legend(axis.get_lines()[1:3], names)

        
    org_line = axis.get_lines()[0]
    org_line.set_color('white')
    axis.set_xlabel('Time (ms)')
    axis.set_ylabel('Activation (dSPM)')
    
    if save_fig:
        figure_path = join(scratch_path, 'figures', 'grand_averages', 'stcs')
        filename = event + '_' + str(int(time)) + '_ms.png'
        brain.save_image(join(figure_path, filename))
        
        filename = event + '_time_course.png'
        brain.time_viewer.mpl_canvas.fig.savefig(join(figure_path, filename))
    
#%% BETA SENSOR SPACE CLUSTER CORRECTED

plt.close('all')

# contrast = ['s1', 's2']
# freqs = [4, 7]
# ISI_type = 'normal'
# times = 0.230
# stat_tmin, stat_tmax = (-0.400, 1.100)
# tmin, tmax = (-0.400, 1.100)
# baseline_string =  'no_baseline_Ss_'
# v_extremum = 1.0
# save_fig = False

contrast = ['o0', 'o15']
freqs = [14, 30]
ISI_type = 'normal'
times = -0.200 #np.arange(-0.300, 0.100, 0.025) #[-0.100, -0.050, 0.090]
stat_tmin, stat_tmax = (-0.400, 0.400)
tmin, tmax = (-0.750, 0.750)
baseline_string =  'no_baseline_'
v_extremum = 0.5
save_fig = False

time_string = sensor_functions.time_string(tmin, tmax)
stat_time_string = sensor_functions.time_string(stat_tmin, stat_tmax)

contrast_name = contrast[0] + '_vs_' + contrast[1]
contrast_string = contrast[0] + '_' + contrast[1] + '_contrast_'
p_threshold = 0.05
freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])



data_path = join(scratch_path, 'MEG', 'grand_averages')
hilbert_path = join('hilbert')

filename = ISI_type + '_' + baseline_string + freq_string  + time_string + \
            'proj_tactile_jitter_hilbert_z_transform_contrast-ave.fif'
evokeds = mne.read_evokeds(join(data_path, hilbert_path, filename))
filename  = 't-test_p_0.05_' + ISI_type + '_' + baseline_string + 'mag_' + \
            contrast_string + freq_string + 'stat_' + stat_time_string + \
            'range_' + time_string + \
            'proj_tactile_jitter_hilbert_z_transform_contrast.npy'
clu = np.load(join(data_path, 'statistics', hilbert_path, filename),
                   allow_pickle=True).item()

figure_path = join(scratch_path, 'figures', 'grand_averages', 'hilbert',
                   'evokeds')        

## find relevant contrast
for evoked in evokeds:
    if evoked.comment == contrast_name:
        evoked.pick_types('mag') ## could be channel_type
        break
    
t_obs = clu['t_obs']
clusters = clu['clusters']
cluster_p_values = clu['cluster_p_values']
H0 = clu['H0']


thresholded_evoked = evoked.copy()
old_data = thresholded_evoked.data
thresholded_evoked._data = np.zeros(thresholded_evoked.data.shape)
sig_indices = np.where(cluster_p_values < p_threshold)[0]
for sig_index in sig_indices[:]:
    cluster = clusters[sig_index]
    thresholded_evoked._data[cluster[1], cluster[0]] = \
        old_data[cluster[1], cluster[0]]

## topo plot
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['axes.titlesize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'x-large'
    
fig = thresholded_evoked.plot_topomap(times=times, scalings=dict(mag=1,
                                                                  grad=1),
                                      units=dict(mag='z-score'),
                                      vmin=-v_extremum, vmax=v_extremum,
                                      colorbar=True,
                                      time_unit='ms')

# lala = evoked.plot_topomap(times=times, scalings=dict(mag=1,
#                                                                   grad=1),
#                                       units=dict(mag='z-score'),
#                                       vmin=-1.0, vmax=1.0, colorbar=True,
#                                       time_unit='ms')
filename = contrast_name + '_' + freq_string + baseline_string + \
        sensor_functions.time_string(stat_tmin, stat_tmax) + \
            'cluster_magnetometers_topolot.png'
if save_fig:
    fig.savefig(join(figure_path, filename))

## try 3d plot ## exporting to FieldTrip

# field_map = mne.make_field_map(thresholded_evoked, None)
# fig = thresholded_evoked.plot_field(field_map, time=0.200)
## still need to figure out how to scale it...
mne.write_evokeds(join(figure_path, 'thresholded-ave.fif'), thresholded_evoked)
## butterfly plot
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['axes.titlesize'] = 'x-large'
mpl.rcParams['axes.labelsize'] = 'large'

n_channels = len(thresholded_evoked.info['ch_names'])
sig_channel_indices = []
for channel_index in range(n_channels):
    if np.sum(thresholded_evoked.data[channel_index, :]) != 0:
        sig_channel_indices.append(channel_index)

sig_evoked = evoked.copy()
sig_evoked.pick(sig_channel_indices)
mag_ylim = (-0.5, 1.8) # z

fig = sig_evoked.plot(spatial_colors=True, scalings=dict(grad=1, mag=1),
                      ylim=dict(mag=mag_ylim),
                      units=dict(mag='z-score',
                                 grad='z-score'),
                      titles=dict(mag='Magnetometers in clusters'),
                      time_unit='ms', xlim=(tmin * 1e3, tmax * 1e3))
axis = fig.axes[0]

for line in axis.lines:
    line.set_linewidth(mpl.rcParams['lines.linewidth'])

axis.plot((sig_evoked.times[0] * 1e3, sig_evoked.times[-1] * 1e3), [0, 0],
          'k--')    
axis.plot([0, 0], mag_ylim, 'k--')
fig.set_size_inches(12, 6)

filename = contrast_name + '_' + freq_string + baseline_string + \
        sensor_functions.time_string(stat_tmin, stat_tmax) + \
        'cluster_magnetometers.png'
if save_fig:
    fig.savefig(join(figure_path, filename))

    
#%% STC CLUSTER CORRECTED

plt.close('all')
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.weight'] = 'bold'

# contrast = ['s1', 's2']
# stat_tmin, stat_tmax = (-0.200, 0.600)
# tmin, tmax = (-0.400, 1.100)
# freqs = [4, 7]
# ISI_type = 'normal'
# proj_string = '' ##'proj_' # 'proj_' or ''
# baseline_string =  'no_baseline_Ss_'
# vertices = get_max_vertices_per_roi('aal', 's1_vs_s2')
# # vertices = get_max_vertices_per_roi('harvard_oxford', 's1_vs_s2')
# save_fig = True
# save_time_course = False
# # roi = 'None'
# # roi = 'Cerebellum_6_L'
# roi = 'Parietal_Inf_L'

contrast = ['o0', 'o15']
stat_tmin, stat_tmax = (-0.400, 0.400)
tmin, tmax = (-0.750, 0.750)
freqs = [14, 30]
ISI_type = 'local'
proj_string = ''#'proj_'
baseline_string =  'no_baseline_'
vertices = get_max_vertices_per_roi('aal', 'o0_vs_o15')
# vertices = get_max_vertices_per_roi('harvard_oxford', 'o0_vs_o15')
save_fig = True
save_time_course = True

roi = 'Putamen_L'
# roi = 'Cerebelum_8_R'
# roi = 'SII_L'

if save_time_course:
    time_course_string = ''
else:
    time_course_string = 'no_time_course'

p_threshold = 0.05
initial_time = None
# initial_time = 300

src = mne.read_source_spaces(join(subjects_dir, 'fsaverage', 'bem',
                                  'volume-7.5mm-src.fif'))

if roi == 'None':
    initial_pos =  None # max value (-22, 0, 7) # Putamen
else:
    initial_pos = src[0]['rr'][vertices[roi], :]

freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])

# contrast_name = contrast[0] + '_vs_' + contrast[1]
contrast_string = ''
for part in contrast:
    contrast_string += part + '_'
contrast_string += 'contrast_'    
    
time_string = sensor_functions.time_string(tmin, tmax)
stat_time_string = sensor_functions.time_string(stat_tmin, stat_tmax)


fig_filename = ISI_type + '_' + contrast_string + baseline_string + \
                freq_string + 'stat_' + \
                stat_time_string + 'range_' + time_string + proj_string + \
                roi + '_' +  time_course_string + '.png'

data_path = join(scratch_path, 'MEG', 'grand_averages')
hilbert_path = join('hilbert', 'stcs', 'common_filter', 'contrasts')
filename = ISI_type + '_ISI_' + baseline_string + \
            'mag_reg_0.0_weight_unit-noise-gain_' + \
            contrast_string + '7.5_mm_' + freq_string + \
            time_string + proj_string + 'tactile_jitter_morph-vl.stc'
stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))
filename = 't-test_p_0.05_' + ISI_type + \
            '_ISI_' + baseline_string + \
            'mag_reg_0.0_weight_unit-noise-gain_' + \
            contrast_string + '7.5_mm_' + freq_string + 'stat_' + \
            stat_time_string +  'range_' + time_string + proj_string  + \
            'tactile_jitter_morph-vl.npy'
clu = np.load(join(data_path, 'statistics', hilbert_path, filename),
                   allow_pickle=True).item()



figure_path = join(scratch_path, 'figures', 'grand_averages', 'hilbert',
                   'stcs')        


t_obs = clu['t_obs']
clusters = clu['clusters']
cluster_p_values = clu['cluster_p_values']
H0 = clu['H0']


thresholded_stc = stc.copy()
thresholded_stc.crop(stat_tmin, stat_tmax)
old_data = thresholded_stc.data
thresholded_stc._data = np.zeros(thresholded_stc.data.shape)
sig_indices = np.where(cluster_p_values < p_threshold)[0]

thresholded_stc._times = thresholded_stc.times * 1e3

if sig_indices.size != 0:

    for sig_index in sig_indices[:]:
        cluster = clusters[sig_index]
        thresholded_stc._data[cluster[1], cluster[0]] = \
            old_data[cluster[1], cluster[0]]
    thresholded_stc._data *= 1e2 # change to percentage                                        
    ## smallest positive value that is not zero
    min_sig = np.min(thresholded_stc.data[thresholded_stc.data > 0])
    # min_sig = np.quantile(thresholded_stc.data[thresholded_stc.data > 0], 0.70)
    mid_sig = np.median(thresholded_stc.data[thresholded_stc.data > 0])
    # mid_sig = np.quantile(thresholded_stc.data[thresholded_stc.data > 0], 0.80)
    max_sig = np.max(thresholded_stc.data[thresholded_stc.data > 0])
    clim = dict(kind='value', lims=(min_sig, mid_sig, max_sig))     
    fig = thresholded_stc.plot(src, src[0]['subject_his_id'],
                                       initial_time=initial_time,
                                       clim=clim, initial_pos=initial_pos
                                       )
    time = str(int(fig.axes[1].get_lines()[1].get_xdata()[0]))
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    fig.show()

    axis_1 = fig.axes[0]
    extent = \
        axis_1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    cbar_text = axis_1.get_children()[-3]
    cbar_text.set_text('Power Change (%): (x1-x2) / (x1+x2) at: ' + time + \
                       ' ms')
    cbar_text.set_fontsize('small')
    
    yticks = np.round(fig.axes[1].get_yticks(), 1)
    labels = [str(yticks[0]), str(yticks[1]), str(yticks[2])]
    fig.axes[1].set_yticklabels(labels)
    fig.axes[1].set_xlabel('Time (ms)')
    fig.axes[1].set_ylabel('Power Change (%)\n(x1-x2) / (x1+x2)')
    print('p = ' + str(np.unique(cluster_p_values)[0]))

    if save_fig:
        if save_time_course:
            fig.savefig(join(figure_path, fig_filename))
        else:
            fig.savefig(join(figure_path, fig_filename),
                        bbox_inches=extent.expanded(1.0, 1.20))
    
else:
    print('Null hypothesis cannot be rejected')
    print('p = ' + str(np.unique(cluster_p_values)[0]))
    
    # mne.viz.plot_volume_source_estimates
  
#%% PLOT CONTRASTS WITH 95 % CI

plt.close('all')
n_subjects = 30

figure_path = join(scratch_path, 'figures', 'grand_averages', 'hilbert',
                   'stcs')    

mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['lines.linewidth'] = 3

contrast = ['o0', 'o15']
freqs = [14, 30]
ISI_type = 'normal'
baseline_string =   'no_baseline_'#'full_baseline_'
times = (-0.750, 0.750)
stat_tmin, stat_tmax = (-0.400, 0.400)
proj_string = ''#'proj_'# ''
ylim = (-3.00, 3.00)
vertices = get_max_vertices_per_roi('aal', 'o0_vs_o15')
# vertices = get_max_vertices_per_roi('harvard_oxford', 'o0_vs_o15')
save_fig = True

# contrast = ['s1', 's2']
# freqs = [4, 7]
# ISI_type = 'normal'
# baseline_string =  'no_baseline_Ss_'#'full_baseline_'
# times = (-0.400, 1.100)
# stat_tmin, stat_tmax = (-0.200, 0.600)
# proj_string = ''
# ylim = (-3.00, 3.00)
# vertices = get_max_vertices_per_roi('aal', 's1_vs_s2')
# # vertices = get_max_vertices_per_roi('harvard_oxford', 's1_vs_s2')
# save_fig = True

p_threshold = 0.05
time_string = sensor_functions.time_string(times[0], times[1])
stat_time_string = sensor_functions.time_string(stat_tmin, stat_tmax)
## from mpl.rcParans['axes-prop_cycle']
colours = ['#1f77b4', '#ff7f0e', '#2ca02c']


# vertices = dict(
#                 Left_SII           = 7732,
#                 Right_SII          = 7746,
                # )
freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])
contrast_string = ''
for part in contrast:
    contrast_string += part + '_'
# contrast_string = contrast[0] + '_' + contrast[1] + '_'

data_path = join(scratch_path, 'MEG', 'grand_averages')
# data_path = join(scratch_path, 'MEG', '0001', '20200124_000000')
hilbert_path = join('hilbert', 'stcs', 'common_filter', 'contrasts')


filename = ISI_type + '_ISI_' + baseline_string + \
            'mag_reg_0.0_weight_unit-noise-gain_' +  contrast_string + \
            'contrast' + '_7.5_mm_' + \
            freq_string + time_string + proj_string + \
            'tactile_jitter_morph-vl.stc'
std_filename = filename[:-7] + '_std' + filename[-7:]
clu_filename = 't-test_p_0.05_' + ISI_type + \
            '_ISI_' + baseline_string + \
            'mag_reg_0.0_weight_unit-noise-gain_' + \
            contrast_string + 'contrast_7.5_mm_' + freq_string + 'stat_' + \
            stat_time_string +  'range_' + time_string + proj_string  + \
            'tactile_jitter_morph-vl.npy'
clu = np.load(join(data_path, 'statistics', hilbert_path, clu_filename),
                   allow_pickle=True).item()
 
stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))
std_stc = mne.read_source_estimate(join(data_path, hilbert_path,
                                        std_filename))
sem_stc = std_stc.copy()
sem_stc._data /= np.sqrt(n_subjects)
stc_plus = stc.copy()
stc_minus = stc.copy()
stc_plus._data  += stats.norm.ppf(0.975) * sem_stc.data
stc_minus._data -= stats.norm.ppf(0.975) * sem_stc.data

t_obs = clu['t_obs']
clusters = clu['clusters']
cluster_p_values = clu['cluster_p_values']
H0 = clu['H0']


thresholded_stc = stc.copy()
thresholded_stc.crop(stat_tmin, stat_tmax)
thresholded_stc._data = np.zeros(thresholded_stc.data.shape)
sig_indices = np.where(cluster_p_values < p_threshold)[0]

for sig_index in sig_indices[:]:
    cluster = clusters[sig_index]
    thresholded_stc._data[cluster[1], cluster[0]] = 1


##
for vertex in vertices:
    fig = plt.figure()
    fig.set_size_inches(9, 6)
    fig_filename = contrast_string + 'contrast_time_courses_' + \
        baseline_string + \
        sensor_functions.time_string(times[0], times[1]) + vertex + '.png'

    data_index = np.where(stc.vertices == vertices[vertex])[0]
    time_indices = np.where(thresholded_stc.data[data_index, :] == 1)[1]
    cluster_times = thresholded_stc.times[time_indices] * 1e3
    colour = colours[0]
    plt.plot(stc.times * 1e3, stc.data[data_index, :].T * 1e2, colour)
    plt.fill_between(stc.times * 1e3,
                     stc_plus.data[data_index, :].T.flatten() * 1e2,
                     stc.data[data_index, :].T.flatten() * 1e2,
                     facecolor=colour,
                     color=colour,
                     alpha=0.2)
    plt.fill_between(stc.times * 1e3,
                     stc_minus.data[data_index, :].T.flatten() * 1e2,
                     stc.data[data_index, :].T.flatten() * 1e2,
                     facecolor=colour,
                     color=colour,
                     alpha=0.2)
    

    plt.title(vertex.replace('_', ' '))
# plt.legend(legend_names, loc='upper right')
    plt.xlim((times[0] * 1e3, times[1] * 1e3))
    plt.ylim(ylim)
    # ylim = fig.axes[0].get_ylim()
    # new_ylim = (0.985 * ylim[0], 1.015 * ylim[1])
    # plt.ylim(new_ylim)
    y_cluster_line = ylim[0] * 0.50
    
             
    # text and cluster_line
    # x_text = cluster_times[len(cluster_times) // 2] # // integer division   
    if cluster_times.shape[0] != 0: 
        plt.plot(cluster_times, np.repeat(y_cluster_line, len(cluster_times)),
                 color=colour, linestyle='-')
        x_text = cluster_times[0]  
        y_text = ylim[0] * 0.65
        plt.text(x_text, y_text, 'Cluster Extent', fontdict=dict(fontsize=14))
    
    plt.xlabel('Time (ms)')
    plt.ylabel('Power Change (%)\n((x1-x2) / (x1+x2))')
    plt.plot(fig.axes[0].get_xlim(), (0, 0), 'k--')
    # plt.plot((0, 0), new_ylim, 'k--')
    if save_fig:
        fig.savefig(join(figure_path, fig_filename))  

#%% PLOT TIME COURSES FOR VERTICES FROM NON-CONTRASTS

plt.close('all')

figure_path = join(scratch_path, 'figures', 'grand_averages', 'hilbert',
                   'stcs')    

mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['lines.linewidth'] = 3

contrast = ['o0', 'o5', 'o15']
freqs = [14, 30]
ISI_type = 'normal'
baseline_string =   'no_baseline_'#'full_baseline_'
times = (-0.750, 0.750)
proj_string = ''#'proj_'# ''
legend_names = ['Omission 0', 'Omission 5', 'Omission 15']
vertices = get_max_vertices_per_roi('aal', 'o0_vs_o15')
# # vertices = get_max_vertices_per_roi('harvard_oxford', 'o0_vs_o15')
save_fig = False
roi = 'Thalamus_L'
# roi = 'Putamen_L'
# roi = 'SI_L'


# contrast = ['s1', 's2', 's3']
# freqs = [4, 7]
# ISI_type = 'normal'
# baseline_string =  'no_baseline_Ss_'#'full_baseline_'
# times = (-0.400, 1.100)
# proj_string = ''
# legend_names = ['First Stimulation', 'Repeated Stimulation (2)',
#                 'Repeated Stimulation (3)']
# vertices = get_max_vertices_per_roi('aal', 's1_vs_s2')
# # vertices = get_max_vertices_per_roi('harvard_oxford', 's1_vs_s2')
# save_fig = True
# # roi = 'Cerebellum_6_L'
# roi = 'Parietal_Inf_L'


time_string = sensor_functions.time_string(times[0], times[1])
freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])
contrast_string = ''
for part in contrast:
    contrast_string += part + '_'
# contrast_string = contrast[0] + '_' + contrast[1] + '_'

data_path = join(scratch_path, 'MEG', 'grand_averages')
# data_path = join(scratch_path, 'MEG', '0001', '20200124_000000')
hilbert_path = join('hilbert', 'stcs', 'common_filter')
fig = plt.figure()
fig.set_size_inches(9, 6)
for condition_index, condition in enumerate(contrast):
    filename = ISI_type + '_ISI_' + baseline_string + \
                'mag_reg_0.0_weight_unit-noise-gain_' + \
                condition + '_7.5_mm_common_filter_' + contrast_string + \
                freq_string + time_string + proj_string + \
                'tactile_jitter_morph-vl.stc'
 
    stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))
    fig_filename = contrast_string + 'time_courses_' + baseline_string + \
        sensor_functions.time_string(times[0], times[1]) + roi + '.png'

    data_index = np.where(stc.vertices == vertices[roi])[0]
    plt.plot(stc.times * 1e3, stc.data[data_index, :].T)
    if condition_index == 0:
        plt.title(roi.replace('_', ' '))
plt.legend(legend_names, loc='upper right')
plt.xlim((times[0] * 1e3, times[1] * 1e3))
ylim = fig.axes[0].get_ylim()
new_ylim = (0.985 * ylim[0], 1.015 * ylim[1])
plt.ylim(new_ylim)
plt.xlabel('Time (ms)')
plt.ylabel('Activation (arbitrary unit)')

plt.plot((0, 0), new_ylim, 'k--')
if save_fig:
    fig.savefig(join(figure_path, fig_filename))
      
    
#%% PLOT TIME COURSES FOR VERTICES FROM NON-CONTRASTS WITH STD!

# plt.close('all')
n_subjects = 30

figure_path = join(scratch_path, 'figures', 'grand_averages', 'hilbert',
                   'stcs')    

mpl.rcParams['font.size'] = 14
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['lines.linewidth'] = 3

# contrast = ['o0', 'o15']
# freqs = [14, 30]
# ISI_type = 'normal'
# baseline_string =   'no_baseline_'#'full_baseline_'
# times = (-0.750, 0.750)
# proj_string = ''#'proj_'# ''
# legend_names = ['Omission 0', 'Omission 5', 'Omission 15']
# save_fig = False

contrast = ['s1', 's2', 's3']
freqs = [4, 7]
ISI_type = 'normal'
baseline_string =  'no_baseline_Ss_'#'full_baseline_'
times = (-0.400, 1.100)
proj_string = ''
legend_names = ['First Stimulation', 'Repeated Stimulation (2)',
                'Repeated Stimulation (3)']
save_fig = False

time_string = sensor_functions.time_string(times[0], times[1])
## from mpl.rcParans['axes-prop_cycle']
colours = ['#1f77b4', '#ff7f0e', '#2ca02c'] 
vertices = dict(
                # Left_Cerebellum_6  =  4562,
                # Right_Cerebellum_6 =  2705,
                # left_thalamus    =  7186,
                # right_thalamus   =  7188,
                # Left_Putamen     =  7207,
                # Right_Putamen    =  7213,
                # Right_Precuneus = 9534
                # left_insula      =  7850,
                # right_insula     =  7859,
                # left_SI          = 12083,
                # right_SI         = 12067,
                Left_Inferior_Parietal_Cortex = 10148,
                # right_IPC        = 10160,
                # left_MCC         =  9580,
                # right_MCC        =  9556,
                # lala               = 4909 # frontal sup orb
                )
freq_string = sensor_functions.hp_lp_string(freqs[0], freqs[1])
contrast_string = ''
for part in contrast:
    contrast_string += part + '_'
# contrast_string = contrast[0] + '_' + contrast[1] + '_'

data_path = join(scratch_path, 'MEG', 'grand_averages')
# data_path = join(scratch_path, 'MEG', '0001', '20200124_000000')
hilbert_path = join('hilbert', 'stcs', 'common_filter')
fig = plt.figure()
fig.set_size_inches(9, 6)
for condition_index, condition in enumerate(contrast):

    filename = ISI_type + '_ISI_' + baseline_string + \
                'mag_reg_0.0_weight_unit-noise-gain_' + \
                condition + '_7.5_mm_common_filter_' + contrast_string + \
                freq_string + time_string + proj_string + \
                'tactile_jitter_morph-vl.stc'
    std_filename = filename[:-7] + '_std' + filename[-7:]
 
    stc = mne.read_source_estimate(join(data_path, hilbert_path, filename))
    std_stc = mne.read_source_estimate(join(data_path, hilbert_path,
                                            std_filename))
    sem_stc = std_stc.copy()
    sem_stc._data /= np.sqrt(n_subjects)
    stc_plus = stc.copy()
    stc_minus = stc.copy()
    stc_plus._data  += sem_stc.data
    stc_minus._data -= sem_stc.data
    for vertex in vertices:
        fig_filename = contrast_string + 'time_courses_' + baseline_string + \
            sensor_functions.time_string(times[0], times[1]) + vertex + '.png'

        data_index = np.where(stc.vertices == vertices[vertex])[0]
        colour = colours[condition_index]
        plt.plot(stc.times * 1e3, stc.data[data_index, :].T, colour)
        plt.fill_between(stc.times * 1e3,
                         stc_plus.data[data_index, :].T.flatten(),
                         stc.data[data_index, :].T.flatten(),
                         facecolor=colour,
                         color=colour,
                         alpha=0.2)
        plt.fill_between(stc.times * 1e3,
                         stc_minus.data[data_index, :].T.flatten(),
                         stc.data[data_index, :].T.flatten(),
                         facecolor=colour,
                         color=colour,
                         alpha=0.2)
        if condition_index == 0 and len(vertices) == 1:
            plt.title(vertex.replace('_', ' '))
# plt.legend(legend_names, loc='upper right')
plt.xlim((times[0] * 1e3, times[1] * 1e3))
ylim = fig.axes[0].get_ylim()
new_ylim = (0.985 * ylim[0], 1.015 * ylim[1])
plt.ylim(new_ylim)
plt.xlabel('Time (ms)')
plt.ylabel('Activation (arbitrary unit)')

plt.plot((0, 0), new_ylim, 'k--')
if save_fig:
    fig.savefig(join(figure_path, fig_filename))  
