%% INITIALIZE

clear variables
restoredefaultpath

fieldtrip_path = '/home/lau/matlab/fieldtrip/';

addpath(fieldtrip_path)
ft_defaults

figures_path = fullfile('/', 'home', 'lau', 'projects', ...
                        'cerebellar_clock', 'scratch', ...
                        'tactile_jitter', 'figures', 'grand_averages', ...
                        'hilbert', 'evokeds', 'thresholded-ave.fif');
                    
%% LOAD DATASET

cfg = [];
cfg.dataset = figures_path;

data = ft_preprocessing(cfg);

%% TOPO PLOT 3D

close all

time_index = data.time{1} == 0.000;
val = data.trial{1}(:, time_index);

figure
ft_plot_topo3d(data.grad.chanpos, val, 'refine', 1)
colorbar