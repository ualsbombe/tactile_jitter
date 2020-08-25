#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:03:19 2020

@author: lau
"""
#%% IMPORTS

from stormdb.cluster import ClusterJob
from subprocess import check_output
from os.path import join
from qsubs import sensor_functions

#%% PATHS

user = check_output('uname -n', shell=True)
if user == b'lau\n':
    project_path = '/home/lau/projects/cerebellar_clock'
elif user[:6] == b'hyades':
    project_path = '/projects/MINDLAB2019_MEG-CerebellarClock'
qsub_path = join(project_path, 'scripts', 'tactile_jitter', 'python', 'qsubs')

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

#%% SUBMIT TO CLUSTER


 #%% FIND EVENTS

for subject_index in range(1, n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'find_events' " + subject + ' ' + \
          date + " 'notch_tactile_jitter_raw.fif' " + \
          "'tactile_jitter-eve.fif' 0"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='event_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% HIGH AND LOWPASS FILTER

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'high_and_lowpass' " + subject + ' ' + \
          date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_raw_sss.fif' 0 0 40"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=4,
                    job_name='hlow_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% EPOCH DATA (for evokeds)
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]


    cmd = "python sensor_functions.py 'epoch_data' " + subject + ' ' + \
          date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_sss-epo.fif' 0 0 40 'tactile_jitter-eve.fif'"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=4,
                    job_name='epo_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit() 

#%% GET EVOKEDS
for subject_index in range(1, n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'get_evokeds' " + subject + \
          ' ' + date + " 'tactile_jitter_sss-epo.fif' " + \
          "'tactile_jitter_sss-ave.fif' 1 0 40 1"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=3,
                    job_name='ave_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% GRAND AVERAGE EVOKEDS

cmd = "python sensor_functions.py 'grand_average_evokeds' " + \
      " 'tactile_jitter_sss-ave.fif' " + \
      "'tactile_jitter_sss-ave.fif' 0 0 40"

cj = ClusterJob(cmd=cmd,
                queue='all.q',
                n_threads=4,
                job_name= 'ga_evokeds',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% WATERSHED
for subject_index in range(1):
    subject = subjects[subject_index]
    # subject = 'fsaverage'

    cmd = "python sensor_functions.py 'watershed' " + subject + " 1"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='wshed_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% SOURCE SPACE, BEM MODEL AND BEM SOLUTION
for subject_index in range(1):
    subject = subjects[subject_index]
    subject = 'fsaverage'

    cmd = "python sensor_functions.py " + \
        "'source_space_and_bem_model_and_bem_solution' " + subject
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='bem_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% VOLUMETRIC SOURCE SPACE
for subject_index in range(1):
    subject = subjects[subject_index]
    subject = 'fsaverage'

    cmd = "python sensor_functions.py " + \
        "'volumetric_source_space' " + subject + " 7.5"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='vsrc_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% MORPH TO FSAVERAGE
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py " + \
        "'morph_to_fsaverage' " + subject + ' ' + date + \
        " 'tactile_jitter_surf-morph.h5' 0 4"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='morph_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% MORPH TO FSAVERAGE VOLUME
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py " + \
        "'morph_to_fsaverage_volume' " + subject + ' ' + date + \
        " 'tactile_jitter_vol-morph.h5' 0 7.5"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='vmph_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% FORWARD SOLUTION
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'forward_solution' " + subject + ' ' + \
          date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_sss-fwd.fif' 0 oct6 4"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='fwd_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% FORWARD SOLUTION VOLUMETRIC
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'forward_solution_volumetric' " + \
          subject + ' ' + date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_sss_vol-fwd.fif' 0 5 4"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='vfwd_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% ESTIMATE COVARIANCE
for subject_index in range(1, 2):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'estimate_covariance' " + subject + \
          ' ' + date + " 'tactile_jitter-epo.fif' " + \
          "'tactile_jitter-cov.fif' 0 0 0"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=3,
                    job_name='cov_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% INVERSE OPERATOR
for subject_index in range(1):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'inverse_operator' " + subject + \
          ' ' + date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_sss-fwd.fif' 'tactile_jitter_sss-cov.fif' " + \
          " 'tactile_jitter_sss-inv.fif' 0 0 40"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='inv_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%%# MINIMUM NORM ESTIMATE
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'minimum_norm_estimate' " + subject + \
          ' ' + date + " 'tactile_jitter_sss-ave.fif' " + \
          "'tactile_jitter_sss-inv.fif' 'tactile_jitter_sss' " + " 1 0 40"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='stc_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%%# MORPH MINIMUM NORM ESTIMATES
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'morph_minimum_norm_estimate' " + \
          subject + ' ' + date + " 'tactile_jitter_sss' " + \
          "'tactile_jitter_surf-morph.h5' 'tactile_jitter_sss_morph' " + \
          " 0 0 40 1 4"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='mstc_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% GRAND AVERAGE STC

cmd = "python sensor_functions.py 'grand_average_stc' " + \
      " 'tactile_jitter_sss_morph' " + \
      "'tactile_jitter_sss_morph' 0 0 40 1 4"

cj = ClusterJob(cmd=cmd,
                queue='short.q',
                n_threads=4,
                job_name= 'ga_stc',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% STATISTICS STC

n_jobs = 12

cmd = "python sensor_functions.py 'statistics_stc' " + \
      " 'tactile_jitter_sss_morph' " + \
      "'tactile_jitter_sss_morph' 0 0 40 [['s1','s2']] 0.05 0.020 0.220 4 " + \
      str(n_jobs)

cj = ClusterJob(cmd=cmd,
                queue='all.q',
                n_threads=n_jobs,
                job_name= 'stat_stc',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% FILTER FOR HILBERT

n_jobs = 1
for subject_index in range(10, 11):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'filter_for_hilbert' " + subject + \
          ' ' + date + ' ' + str(split_recording) + \
          " 'tactile_jitter_raw.fif' " + \
          "'tactile_jitter_raw.fif' 1 [25] [35] " + str(n_jobs)
          #"[4,8,14] [7,12,30]"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=2, #n_jobs, ## too high on memory
                    job_name='fhil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% HILBERT TRANSFORM

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'hilbert_transform' " + subject + \
          ' ' + date + \
          " 'tactile_jitter_raw.fif' " + \
          "'tactile_jitter_hilbert-epo.fif' 0 " + \
          "[4,8,14,60] [7,12,30,85] no_baseline_Ss no_baseline_Ss " + \
          " -0.400 1.100 " + \
          " 'tactile_jitter-eve.fif' normal"
         
    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=4,
                    job_name='hilb_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% EVOKED HILBERT 

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'evoked_hilbert' " + subject + \
          ' ' + date + \
          " 'tactile_jitter_hilbert-epo.fif' " + \
          "'tactile_jitter_hilbert-ave.fif' 0 1 1 " + \
          "[14] [30] normal"
         
    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=4,
                    job_name='ehil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()    
    
#%% Z TRANSFORM HILBERT CONTRASTS

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'z_transform_wilcoxon_contrasts' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_hilbert-epo.fif' " + \
          "'tactile_jitter_hilbert_z_transform_contrast-ave.fif' 1 1 1 " + \
          "[14] [30] -0.750 0.750 " + \
          "[['o0','o5'],['o0','o15'],['o5','o15']," + \
          "['o0','ns'],['o5','ns'],['o15','ns']] " + \
          "normal no_baseline 4"

                   # "[['s1','s2'],['s2','s3']] " + \

    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=2,
                    job_name='zhic_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()     
          
    
#%% GRAND AVERAGE EVOKED HILBERT

cmd = "python sensor_functions.py 'grand_average_evoked_hilbert' " + \
      " 'proj_tactile_jitter_hilbert_z_transform_contrast-ave.fif' " + \
      "'proj_tactile_jitter_hilbert_z_transform_contrast-ave.fif' " + \
      "1 [14] [30] -0.750 0.750 normal no_baseline"

cj = ClusterJob(cmd=cmd,
                queue='short.q',
                n_threads=4,
                job_name= 'ga_ev_hilb',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% STATISTICS EVOKED HILBERT

h_freqs = [[4], [8], [14], [60]]
l_freqs = [[7], [12], [30], [85]]
contrasts = ["[o0,o5]", "[o0,o15]", "[o5,o15]",
             "[o0,ns]", "[o5,ns]", "[o15,ns]",
             "[s1,s2]", "[s2,s3]"]
n_jobs = 8

for (h_freq, l_freq) in zip(h_freqs[2:3], l_freqs[2:3]):
    for contrast in contrasts[3:6]:

        cmd = "python sensor_functions.py " + \
          " 'statistics_evoked_hilbert' " + \
          " 'proj_tactile_jitter_hilbert_z_transform_contrast-ave.fif' " + \
          "'proj_tactile_jitter_hilbert_z_transform_contrast.npy' 1 " + \
              str(h_freq) + " " + \
              str(l_freq) + " mag " + \
              str([str(contrast)]) + " " + \
              " 0.05 -0.750 0.750 -0.750 0.750 normal no_baseline " + \
              str(n_jobs)
                  
        cj = ClusterJob(cmd=cmd,
                        queue='short.q',
                        n_threads=n_jobs,
                        job_name= 'stat_evh',
                        proj_name=proj_name,
                        working_dir=qsub_path)
        cj.submit()  

    
#%% BEAMFORMER CONTRAST HILBERT

# time.sleep(3600)
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'beamformer_contrast_hilbert' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_raw.fif' " + \
          "'tactile_jitter_hilbert-epo.fif' 'tactile_jitter-vl-stc.h5' " + \
          "0 1 0 0 [14] [30] -0.750 0.750 7.5 0.00 " + \
          "unit-noise-gain mag normal " + \
          "[['o0','o15']] no_baseline 4 1"    
          # "[['s1','s2'],['s2','s3']] no_baseline_Ss 4"    
          # "[['s1','s2','s3']] no_baseline_Ss"

          # "[['o0','o5','o15']] no_baseline"

          # "[['o0','o5'],['o0','o15'],['o5','o15']]"
           # "[['o0','ns'],['o5','ns'],['o15','ns']]"# + \
          # ",['s4_0','s5_0'],['s5_0','s6_0']" + \
          # ",['s4_5','s5_5'],['s5_5','s6_5']" + \
          # ",['s4_15','s5_15'],['s5_15','s6_15']" + \
          # ",['s4_0','s4_5'],['s4_5','s4_15']" + \
          # ",['s5_0','s5_5'],['s5_5','s5_15']" + \
          # ",['s6_0','s6_5'],['s6_5','s6_15']" + \
          # ",['o0','o5'],['o0','o15'],['o5','o15']]" #+ \
                        # "[['o0','o15'],['o5','o15']]"
                                  # "[['s1','s2'],['s2','s3']]"


    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=2,
                    job_name='bchil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

    
#%% MORPH BEAMFORMER CONTRAST HILBERT

# time.sleep(24 * 3600)
for subject_index in range(18, n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'morph_beamformer_contrast_hilbert' " + \
          subject + ' ' + date + " 'tactile_jitter-vl-stc.h5' " + \
          "'tactile_jitter_vol-morph.h5' " + \
          "'tactile_jitter_morph-vl-stc.h5' " + \
          "0 0 [14] [30] -0.750 0.750 7.5 0.00 " + \
          "unit-noise-gain mag normal " + \
          "[['o0','o15']] no_baseline 1"    
          # "[['o0','o5'],['o0','o15'],['o5','o15']] no_baseline"

          # "[['s1','s2'],['s2','s3']] no_baseline_Ss"


          # "[['o0','o5'],['o0','o15'],['o5','o15']] no_baseline"
          # "[['o0','o5','o15']] no_baseline"

          


          # "[['o0','o5'],['o0','o15'],['o5','o15']] 0"



          # ",['o0','ns'],['o5','ns'],['o15','ns']" + \
          # ",['s4_0','s5_0'],['s5_0','s6_0']" + \
          # ",['s4_5','s5_5'],['s5_5','s6_5']" + \
          # ",['s4_15','s5_15'],['s5_15','s6_15']" + \
          # ",['s4_0','s4_5'],['s4_5','s4_15']" + \
          # ",['s5_0','s5_5'],['s5_5','s5_15']" + \
          # ",['s6_0','s6_5'],['s6_5','s6_15']]"

    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=1,
                    job_name='mbchi_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()    
    


#%% GRAND AVERAGE BEAMFORMER CONTRAST HILBERT

cmd = "python sensor_functions.py " + \
      "'grand_average_beamformer_contrast_hilbert' " + \
      " 'tactile_jitter_morph-vl.stc' " + \
      "'tactile_jitter_morph-vl.stc' 1 1 [14] [30] " + \
      "-0.750 0.750 7.5 0.00 unit-noise-gain mag normal " + \
      "[['o0','o15']] no_baseline 1"      
      # "[['s1','s2'],['s2','s3']] no_baseline_Ss 1"   
      # "[['o0','o5'],['o0','o15'],['o5','o15']] no_baseline 1"

    # "[['o0','o5'],['o0','o15'],['o5','o15']] no_baseline 1"
       # "[['o0','o5','o15']] no_baseline 1"

      # "[['o0','o5'],['o0','o15'],['o5','o15']]" #+ \
      # "[['o0','ns'],['o5','ns'],['o15','ns']" + \

          
cj = ClusterJob(cmd=cmd,
                queue='short.q',    
                n_threads=2,
                job_name= 'ga_bchil',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% STATISTICS BEAMFORMER CONTRAST HILBERT

h_freqs = [[4], [8], [14], [60]]
l_freqs = [[7], [12], [30], [85]]
contrasts = ["[o0,o5]", "[o0,o15]", "[o5,o15]",
             "[o0,ns]", "[o5,ns]", "[o15,ns]",
             "[s1,s2]", "[s2,s3]"]
n_jobs = 12

for (h_freq, l_freq) in zip(h_freqs[2:3], l_freqs[2:3]):
    for contrast in contrasts[1:2]:

        cmd = "python sensor_functions.py " + \
              " 'statistics_beamformer_contrast_hilbert' " + \
              " 'tactile_jitter_morph-vl.stc' " + \
              "'tactile_jitter_morph-vl.npy' 1 1 " + str(h_freq) + " " + \
              str(l_freq) + " 7.5 0.00 unit-noise-gain mag " + \
              str([str(contrast)]) + " " + \
              "0.05 -0.750 0.750 -0.400 0.400 stc normal t-test " + \
              "no_baseline " + str(n_jobs)
                  
        cj = ClusterJob(cmd=cmd,
                        queue='highmem.q',
                        n_threads=n_jobs,
                        job_name= 'stat_bc',
                        proj_name=proj_name,
                        working_dir=qsub_path)
        cj.submit()
        
#%% SAVE NIFTI
for subject_index in range(1, n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'save_nifti' " + \
          subject + ' ' + date + " 'tactile_jitter_morph-vl.stc' " + \
          " 'tactile_jitter_morph-vl.nii' " + \
          "0 [14] [30] -0.750 0.750 7.5 0.00 " + \
          "unit-noise-gain mag normal " + \
          "[['o0','o5'],['o0','o15'],['o5','o15']] no_baseline "

    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='nifti_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()           
        
#%% SAVE NIFTI STAT

h_freqs = [4, 8, 14, 25, 60]
l_freqs = [7, 12, 30, 35, 85]
contrasts = ["[o0,o5]", "[o0,o15]", "[o5,o15]",
             "[o0,ns]", "[o5,ns]", "[o15,ns]",
             "[s1,s2]"]

for (h_freq, l_freq) in zip(h_freqs[0:1], l_freqs[0:1]):
    for contrast in contrasts[6:7]:

        cmd = "python sensor_functions.py " + \
              " 'save_nifti_stat' " + \
              " 'tactile_jitter_morph-vl.npy' " + \
              "'tactile_jitter_morph-vl.stc' " + \
              "'tactile_jitter_morph-vl.nii' " + \
              "0 " + str(h_freq) + " " + \
              str(l_freq) + " 7.5 0.00 unit-noise-gain mag " + \
              str([str(contrast)]) + " " + \
              "0.05 -0.400 1.100 -0.200 0.600 stc normal t-test "
                  
        cj = ClusterJob(cmd=cmd,
                        queue='short.q',
                        n_threads=1,
                        job_name= 'nifti_stat',
                        proj_name=proj_name,
                        working_dir=qsub_path)
        cj.submit()       

#%% READ AND RESAVE EPOCHS
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'read_and_resave_epochs' " + \
          subject + ' ' + date + " 'tactile_jitter_hilbert-epo.fif' " + \
          " 'tactile_jitter_hilbert-epo.fif' " + \
          "0 [4,8,14,60] [7,12,30,85] normal -0.400 1.100 0"

    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=4,
                    job_name='resav_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()  
        
#%% PLOT ALIGNMENT

for subject_index in range(1):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'plot_alignment' " + subject + \
          ' ' + date +  ' ' + str(split_recording) + \
          " 'tactile_jitter_raw.fif' " + \
          "'alignment.png' 0"
         
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='align_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()    

#%% PLOT ALIGNMENT (locally)

for subject_index in range(9, 10):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    sensor_functions.plot_alignment(subject, date, bool(split_recording),
                                    'tactile_jitter_raw.fif', 'second_try_alignment.svg',
                                    False)    
        
#%% PLOT EVOKED HILBERT

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'plot_evoked_hilbert' " + subject + \
          ' ' + date + \
          " 'tactile_jitter_hilbert-ave.fif' " + \
          "'evokeds.png' 1 [4,8,14,25,60] [7,12,30,35,85]"
         
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='phil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% PLOT PROJS EVOKED HILBERT

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'plot_projs_evoked_hilbert' " + \
          subject + ' ' + date + \
          " 'proj_tactile_jitter_hilbert-ave.fif' " + \
          "'evokeds_proj.png' 0 [14] [30] normal"
         
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='pphil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

    
#%% PLOT ALL TOGETHERS EVOKED HILBERT

cmd = "python sensor_functions.py 'plot_all_together_evoked_hilbert' " + \
      "'tactile_jitter_hilbert-ave.fif' " + \
      "'evokeds.png' 0 [14] [30] mag mean"
     
cj = ClusterJob(cmd=cmd,
                queue='all.q',
                n_threads=4,
                job_name='pathil',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()


#%% PLOT BEAMFORMER CONTRAST HILBERT

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'plot_beamformer_contrast_hilbert' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_morph-vl.stc' " + \
          "'beamfomer_contrast.png' 0 [14] [30] -0.750 0.750 7.5 0.00 ' "+ \
          "unit-noise-gain mag normal " + \
          "[['o0','o5'],['o0','o15'],['o5','o15']] 0.000 " + \
          "['0.025','-0.055','-0.048']"
         
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='pbhilc_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% PLOT BEAMFORMER TIME COURSES HILBERT

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py " + \
          "'plot_beamformer_time_courses_hilbert' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_morph-vl.stc' " + \
          "'beamfomer_time_courses.png' 0 0 [14] [30] -0.750 0.750 7.5 " + \
          "0.00 unit-noise-gain mag normal " + \
          "[['o0','o15']] no_baseline Right_Cerebellum_6"
         
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='ptchi_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()    


#%% PLOT WHITENED DATA COV

for subject_index in range(1):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'plot_whitened_data_cov' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_raw.fif' " + \
          "'tactile_jitter_hilbert-epo.fif' 'whitened_cov.png' " + \
          "0 1 1 [4,8,14,25,60] [7,12,30,35,85] 7.5 0.00 " + \
          "unit-noise-gain both " + \
           "[['o0','o5'],['o0','o15'],['o5','o15']]"# + \
          # ",['s4_0','s5_0'],['s5_0','s6_0']" + \
          # ",['s4_5','s5_5'],['s5_5','s6_5']" + \
          # ",['s4_15','s5_15'],['s5_15','s6_15']" + \
          # ",['s4_0','s4_5'],['s4_5','s4_15']" + \
          # ",['s5_0','s5_5'],['s5_5','s5_15']" + \
          # ",['s6_0','s6_5'],['s6_5','s6_15']" + \
          # ",['o0','o5'],['o0','o15'],['o5','o15']]" #+ \
                        # "[['o0','o15'],['o5','o15']]"
                                  # "[['s1','s2'],['s2','s3']]"



    
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='pwcov_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()   