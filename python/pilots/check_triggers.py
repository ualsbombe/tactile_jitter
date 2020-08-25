#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:48:23 2019

@author: lau
"""

## check events

import mne
import numpy as np
from os import chdir

sds = ['0001/20190923_000000/',
      '0002/20190930_000000/',
      '0003/20190930_000000/',
      '0004/20191008_000000/',
      '0005/20191008_000000/']

subject_index = 4
sd = sds[subject_index]

path = '/home/lau/mounts/raw/sorted/MINDLAB2019_MEG-CerebellarClock/' + sd + \
       'MEG/001.tactile_jitter_raw/files'
       
       
chdir(path)

raw = mne.io.read_raw_fif('tactile_jitter_raw.fif')
events = mne.find_events(raw, 'STI101', min_duration=0.002)

#%% 

values = np.unique(events[:, 2])
total_values = 0

for value in values:
    n_values = sum(events[:, 2] == value)
    print('There are ' + str(n_values) + \
          ' instances of value: ' + str(value))
    total_values += n_values
    
print('There were ' + str(total_values) + ' triggers')

#%% 

combined_values = [[1], [3], [5], [18], [23], [25], [27],
                   [71, 87],
                   [73, 89],
                   [75, 91],
                   [81, 97],
                   [83, 99],
                   [85, 101],
                   [92, 154],
                   [102, 164],
                   [112, 114, 116, 118, 120],
                   [122, 124, 126, 128, 130],
                   [132, 134, 136, 138, 140],
                   [174, 176, 178, 180, 182],
                   [184, 186, 188, 190, 192],
                   [194, 196, 198, 200, 202]]

for values in combined_values:
    n_values = 0
    for value in values:
        n_values += sum(events[:, 2] == value)
    print('There are ' + str(n_values) + \
          ' instances of value: ' + str(value))