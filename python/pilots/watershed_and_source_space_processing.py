#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:18:43 2019

@author: lau
"""

#%% setup

import mne

subjects = ['0001',
            '0002',
            '0003',
            '0004',
            '0005']

subjects_dir = '/projects/' + \
                'MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/' + \
                'freesurfer'
                
#%% watershed

for subject in subjects:
              
    mne.bem.make_watershed_bem(subject, subjects_dir)                