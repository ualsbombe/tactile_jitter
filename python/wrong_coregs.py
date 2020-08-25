#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 09:56:51 2020

@author: lau
"""

from os.path import join
import mne

raw_path = '/home/lau/projects/cerebellar_clock/raw'
save_path = '/home/lau/mounts/hyades/public/Lubell/laus_subjects/'
subjects = [
            '0001',
            '0002',
            '0007',
            '0011'
           ]
dates = [
         '20200124_000000',
         '20200121_000000',
         '20200120_000000',
         '20200122_000000'
        ]

tags = [
        'good',
        'good',
        'bad',
        'bad'
        ]
    
for subject_index, subject in enumerate(subjects):
    date = dates[subject_index]
    tag = tags[subject_index]
    
    full_path = join(raw_path, subject, date, 'MEG', '001.tactile_jitter_raw',
                     'files', 'tactile_jitter_raw.fif')
    raw = mne.io.read_raw_fif(full_path)
    raw.crop(0.0, 5.0)
    new_path = join(save_path, subject + '-' + tag, subject + '-raw.fif')
    raw.save(new_path)
    