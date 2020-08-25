#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 13:40:23 2020

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

#%% NOTCH FILTER

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'notch_filter' " + subject + ' ' + \
          date + ' ' + \
          str(split_recording) + " 'tactile_jitter_raw.fif' " + \
          "'notch_tactile_jitter_raw.fif' 0"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=4,
                    job_name='notch_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% EPOCH DATA (for tfr)
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]


    cmd = "python sensor_functions.py 'epoch_data_tfr' " + subject + ' ' + \
          date + ' ' + str(split_recording) + \
        " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_tfr_sss-epo.fif' 1 'tactile_jitter-eve.fif'"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=4,
                    job_name='eptfr_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% RUN ICA
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]


    cmd = "python sensor_functions.py 'run_ICA' " + subject + ' ' + \
          date + \
          " 'tactile_jitter_raw_sss.fif' " + \
          "'tactile_jitter_raw_sss-ica.fif' 1 0 0"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=12,
                    job_name='ica_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()      

#%% GET TFR

for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'get_tfr' " + subject + \
          ' ' + date + " 'tactile_jitter_tfr_sss-epo.fif' " + \
          "'tactile_jitter_sss-tfr.h5' 1 1"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=4,
                    job_name='tfr_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% GET TFR (locally)

for subject_index in range(29, 30):
    subject = subjects[subject_index]
    date = dates[subject_index]

    sensor_functions.get_tfr(subject, date,
                             input_file='tactile_jitter_tfr_sss-epo.fif',
                             output_file='tactile_jitter_sss-tfr.h5',
                             overwrite=False,
                             collapse_cond=True)

#%% GRAND AVERAGE TFR

cmd = "python sensor_functions.py 'grand_average_tfr' " + \
      " 'tactile_jitter_sss-tfr.h5' " + \
      "'tactile_jitter_sss-tfr.h5' 1 1"

cj = ClusterJob(cmd=cmd,
                queue='highmem.q',
                n_threads=4,
                job_name= 'ga_tfr',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%%# BEAMFORMER
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'beamformer' " + subject + \
          ' ' + date + " 'tactile_jitter_tfr_sss-epo.fif' " + \
          "'tactile_jitter_sss-vl.stc' " + " 1 1 14 18 0.500 1.000 " + \
          "['s1','s2','s3','s4_0','s5_0','s6_0'" + \
          ",'s4_5','s5_5','s6_5','s4_15','s5_15','s6_15']"
          #['o0','o5','o15']"
    cj = ClusterJob(cmd=cmd,
                    queue='all.q',
                    n_threads=1,
                    job_name='beam_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%%# BEAMFORMER CONTRAST
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'beamformer_contrast' " + subject + \
          ' ' + date + " 'tactile_jitter_tfr_sss-epo.fif' " + \
          "'tactile_jitter_sss-vl.stc' " + " 0 1 8 12 -0.100 0.300 " + \
              "[['o0','o15']]"
        # "[['o0','ns'],['o5','ns'],['o15','ns']" + \
        # ",['o0','o5'],['o0','o15'],['o5','o15']]"
        
                 # "[['s1','s2'],['s1','s3'],['s2','s3'],['s1','ns']]" 

    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=2,
                    job_name='beamc_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%%# MORPH BEAMFORMER
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'morph_beamformer' " + \
          subject + ' ' + date + " 'tactile_jitter_sss-vl.stc' " + \
          "'tactile_jitter_vol-morph.h5' " + \
          "'tactile_jitter_sss_morph-vl.stc' " + \
          " 0 8 12 0.100 0.800 " + \
          "['s1','s2','s3','s4_0','s5_0','s6_0'" + \
          ",'s4_5','s5_5','s6_5','s4_15','s5_15','s6_15']"
          #['o0','o5','o15']"
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='mbeam_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%%# MORPH BEAMFORMER CONTRASTS
for subject_index in range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'morph_beamformer_contrast' " + \
          subject + ' ' + date + " 'tactile_jitter_sss-vl.stc' " + \
          "'tactile_jitter_vol-morph.h5' " + \
          "'tactile_jitter_sss_morph-vl.stc' " + \
          " 0 8 12 -0.100 0.300 " + \
              "[['o0','o15']]"
          # "[['o0','ns'],['o5','ns'],['o15','ns']" + \
          # ",['o0','o5'],['o0','o15'],['o5','o15']]"
                 # "[['s1','s2'],['s1','s3'],['s2','s3'],['s1','ns']]"

        
    cj = ClusterJob(cmd=cmd,
                    queue='short.q',
                    n_threads=1,
                    job_name='mbeac_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()


#%% GRAND AVERAGE BEAMFORMER

cmd = "python sensor_functions.py 'grand_average_beamformer' " + \
      " 'tactile_jitter_sss_morph-vl.stc' " + \
      "'tactile_jitter_sss_morph-vl.stc' 0 14 18 0.500 1.000 " + \
      "['s1','s2','s3','s4_0','s5_0','s6_0'" + \
          ",'s4_5','s5_5','s6_5','s4_15','s5_15','s6_15']"
          #['o0','o5','o15']"
cj = ClusterJob(cmd=cmd,
                queue='highmem_short.q',
                n_threads=1,
                job_name= 'ga_beam',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% GRAND AVERAGE BEAMFORMER CONTRAST

cmd = "python sensor_functions.py 'grand_average_beamformer_contrast' " + \
      " 'tactile_jitter_sss_morph-vl.stc' " + \
      "'tactile_jitter_sss_morph-vl.stc' 0 8 12 -0.100 0.300 " + \
      "[['o0','o15']]"
        # "[['o0','ns'],['o5','ns'],['o15','ns']" + \
        # ",['o0','o5'],['o0','o15'],['o5','o15']]"
#"[['s1','s2'],['s1','s3'],['s2','s3'],['s1','ns']]"
cj = ClusterJob(cmd=cmd,
                queue='highmem_short.q',
                n_threads=1,
                job_name= 'ga_beamc',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% STATISTICS BEAMFORMER CONTRAST

cmd = "python sensor_functions.py 'statistics_beamformer_contrast' " + \
      " 'tactile_jitter_sss_morph-vl.stc' " + \
      "'tactile_jitter_sss_morph-vl.npy' 0 8 12 -0.100 0.300 " + \
       "[['o0','o15']]"
        #"[['o0','ns'],['o5','ns'],['o15','ns']" + \
        #",['o0','o5'],['o0','o15'],['o5','o15']]"
#"[['s1','s2'],['s1','s3'],['s2','s3'],['s1','ns']]"
cj = ClusterJob(cmd=cmd,
                queue='highmem_short.q',
                n_threads=1,
                job_name= 'stat_bc',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()

#%% Z TRANSFORM HILBERT 

for subject_index in range(3, n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'z_transform_wilcoxon' " + subject + \
          ' ' + date + \
          " 'tactile_jitter_hilbert-epo.fif' " + \
          "'tactile_jitter_hilbert_z_transform-ave.fif' 1 1 1 " + \
          "[14] [30] normal "
         
    cj = ClusterJob(cmd=cmd,
                    queue='highmem.q',
                    n_threads=4,
                    job_name='zhil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit() 


#%% TFR Z TRANSFORM

for subject_index in range(n_subjects): 
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'tfr_z_transform' " + subject + \
          ' ' + date + \
          " 'tactile_jitter_hilbert_z_transform-ave.fif' " + \
          "'tactile_jitter_hilbert_z_transform-tfr.h5' 1 " + \
          "[14] [30] normal "
         
    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=1,
                    job_name='tfr_z_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()
    
#%% GRAND AVERAGE TFR Z TRANSFORM

cmd = "python sensor_functions.py 'grand_average_tfr_z_transform' " + \
      " 'tactile_jitter_hilbert_z_transform-tfr.h5' " + \
      "'tactile_jitter_hilbert_z_transform-tfr.h5' 0 [14] [30] normal"

cj = ClusterJob(cmd=cmd,
                queue='highmem.q',
                n_threads=6,
                job_name= 'ga_tfr_z',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()    

#%% STATISTICS TFR Z TRANSFORM

h_freqs = [[4], [8], [14], [60]]
l_freqs = [[7], [12], [30], [85]]
contrasts = ["[o0,o5]", "[o0,o15]", "[o5,o15]",
             "[o0,ns]", "[o5,ns]", "[o15,ns]",
             "[s1,s2]"]
n_jobs = 12

for (h_freq, l_freq) in zip(h_freqs[2:3], l_freqs[2:3]):
    for contrast in contrasts[1:2]:

        cmd = "python sensor_functions.py " + \
              " 'statistics_tfr_z_transform' " + \
              " 'tactile_jitter_hilbert_z_transform-tfr.h5' " + \
              "'tactile_jitter_hilbert_z_transform.npy' 0 " + \
              str(h_freq) + " " + \
              str(l_freq) + " " + contrast + " normal " + \
              str(n_jobs)
                  
        cj = ClusterJob(cmd=cmd,
                        queue='highmem.q',
                        n_threads=n_jobs,
                        job_name= 'stat_z',
                        proj_name=proj_name,
                        working_dir=qsub_path)
        cj.submit()
        
#%% GRAND AVERAGE EVOKED HILBERT CONTRAST

cmd = "python sensor_functions.py 'grand_average_evoked_hilbert_contrast' " + \
      " 'tactile_jitter_hilbert-ave.fif' " + \
      "'tactile_jitter_hilbert_contrast-ave.fif' 0 [14] [30] " + \
      "[['o0','o5'],['o0','o15'],['o5','o15']]"

cj = ClusterJob(cmd=cmd,
                queue='highmem_short.q',
                n_threads=2,
                job_name= 'ga_ev_hilc',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit() 

#%% BEAMFORMER HILBERT

for subject_index in [14]:#range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]
    split_recording = split_recordings[subject_index]

    cmd = "python sensor_functions.py 'beamformer_hilbert' " + \
          subject + ' ' + date + \
          " 'tactile_jitter_raw.fif' " + \
          "'tactile_jitter_hilbert-epo.fif' 'tactile_jitter-vl.stc' " + \
          "'tactile_jitter-vl.stc' " +\
          "0 1 [4,8,14,60] [7,12,30,85] 7.5 " + \
          "['ns']"
           # "['s4_0','s4_5','s4_15','s5_0','s5_5','s5_15'" + \
           # ",'s6_0','s6_5','s6_15','o0','o5','o15','ns']"
          #"['o0','o5','o15','ns']"
          # "['s1','s2','s3']"
          # "['s4_0','s4_5','s4_15','s5_0','s5_5','s5_15'" + \
          # ",'s6_0','s6_5','s6_15']"

          # "['s1','s2','s3']"   

         
          #"['o0','o5','o15']
          
    cj = ClusterJob(cmd=cmd,
                    queue='long.q',
                    n_threads=4,
                    job_name='bhil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% BEAMFORMER CONTRAST HILBERT (locally)
 
for subject_index in range(1):
    subject = subjects[subject_index]
    date = dates[subject_index]
    
    sensor_functions.beamformer_contrast_hilbert(subject, date,
                                                 input_file_1='tactile_jitter_raw.fif',
                                                 input_file_2='tactile_jitter_hilbert-epo.fif',
                                                 output_file='tactile_jitter-vl.stc',
                                                 overwrite=False,
                                                 collapse_cond=True,
                                                 h_freqs=[4],
                                                 l_freqs=[7],
                                                 spacing=7.5,
                                                 reg=0.00,
                                                 weight_norm='None',
                                                 contrasts=[['o0', 'o15']])
   

#%%# MORPH BEAMFORMER HILBERT
for subject_index in [24]:#range(n_subjects):
    subject = subjects[subject_index]
    date = dates[subject_index]

    cmd = "python sensor_functions.py 'morph_beamformer_hilbert' " + \
          subject + ' ' + date + " 'tactile_jitter-vl.stc' " + \
          "'tactile_jitter_vol-morph.h5' " + \
          "'tactile_jitter_morph-vl.stc' " + \
          " 0 [4,8,14,60] [7,12,30,85] 7.5 " + \
          "['o0','o5','o15','s1','s2','s3','s4_0','s5_0','s6_0'" + \
          ",'s4_5','s5_5','s6_5','s4_15','s5_15','s6_15','ns']"
    cj = ClusterJob(cmd=cmd,
                    queue='highmem_short.q',
                    n_threads=2,
                    job_name='mbhil_' + subject,
                    proj_name=proj_name,
                    working_dir=qsub_path)
    cj.submit()

#%% GRAND AVERAGE BEAMFORMER HILBERT

cmd = "python sensor_functions.py 'grand_average_beamformer_hilbert' " + \
      " 'tactile_jitter_morph-vl.stc' " + \
      "'tactile_jitter_morph-vl.stc' 0 [4,8,14,60] [7,12,30,85] 7.5  " + \
      "['o0','o5','o15','s1','s2','s3','s4_0','s5_0','s6_0'" + \
          ",'s4_5','s5_5','s6_5','s4_15','s5_15','s6_15','ns']"
cj = ClusterJob(cmd=cmd,
                queue='highmem_short.q',
                n_threads=1,
                job_name= 'ga_bhil',
                proj_name=proj_name,
                working_dir=qsub_path)
cj.submit()        