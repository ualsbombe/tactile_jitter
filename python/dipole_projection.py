#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 09:23:35 2020

@author: lau
"""


import mne
from os.path import join
import numpy as np

path = '/home/lau/projects/cerebellar_clock/scratch/tactile_jitter'
meg_path = join(path, 'MEG')
subjects_dir = join(path, 'freesurfer')
subject = '0002'
date = '20200121_000000'

evoked = mne.read_evokeds(join(meg_path, subject, date,
                               'tactile_jitter-ave.fif'),
                          's1', proj=False)
evoked.del_proj()
evoked.crop(0.040, 0.040)

cov = mne.read_cov(join(meg_path, subject, date, 'tactile_jitter-cov.fif'))
bem = mne.read_bem_solution(join(subjects_dir, subject, 'bem',
                                 subject + '-5120-bem-sol.fif'))
trans = mne.read_trans(join(subjects_dir, subject, 'bem',
                            subject + '-trans.fif'))

dipole, residuals = mne.fit_dipole(evoked, cov, bem, trans)

dipole.plot_locations(trans, subject)
cerebellar_dipole = dipole.copy()
cerebellar_dipole.pos = np.array([[0.025, -0.030, 0.000]])
cerebellar_dipole.ori = np.array([[0, 0, 1]])
cerebellar_dipole.plot_locations(trans, subject)

fwd, stc = mne.make_forward_dipole(cerebellar_dipole, bem, evoked.info, trans)

simulated_evoked = mne.simulation.simulate_evoked(fwd, stc, evoked.info,
                                                  cov, evoked.nave)
simulated_evoked.plot_topomap(times=0.040)