%% INITIALIZE

clear variables
restoredefaultpath

fieldtrip_path = '/home/lau/matlab/fieldtrip/';

addpath(fieldtrip_path)
ft_defaults

%% LOAD MRI AND ATLAS

mri = ft_read_mri(fullfile(fieldtrip_path, ...
                            'template', 'anatomy', 'single_subj_T1.nii'));
                        
mri.coordsys = 'mni';                        
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                            'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));

                        
%% PLOT IT

% close all

cfg = [];
cfg.atlas = atlas;
cfg.location = [-31 -32 -33];

ft_sourceplot(cfg, mri);

fig = gcf;
set(fig, 'units', 'normalized', 'outerposition', [0 0 0.5 1]);