%%  CLEAR AND ADDPATHS

clear variables
restoredefaultpath

addpath /home/lau/matlab/fieldtrip/
ft_defaults

% sd = 'NatMEG_0245/190809/';
% sd = 'NatMEG_0177/190619/';
sd = 'NatMEG_0521/190808/';

raw_path = ['/archive/20068_tactile_prediction/MEG/' sd];
data_path = ['/home/lau/analyses/tactile_prediction_jitter/data/' ...
             sd];
mr_path = ['/home/lau/analyses/tactile_prediction_jitter/data/MRI/' ...
            sd(1:11) '/mri/T1-neuromag/sets'];
        
        
%% read mr

% mr_file =  ['/home/lau/analyses/tactile_prediction_jitter/data/MRI/dicoms/' ...
%             sd(1:11) 'T1/00000001.dcm'];
mr_file =  ['/home/lau/analyses/tactile_prediction_jitter/data/MRI/dicoms/' ...
            sd(1:11) '/T1/IM-0001-392433672-0001.dcm'];

mri = ft_read_mri(fullfile(mr_file));

%% plot

cfg = [];

ft_sourceplot(cfg, mri);

%% align with fiducals

cfg          = [];
cfg.method   = 'interactive';
cfg.coordsys = 'neuromag';

mri_aligned_fiducials = ft_volumerealign(cfg, mri);

%% alignment

headshape = ft_read_headshape(fullfile(data_path, 'tactile_jitter_avg.fif'));

cfg                     = [];
cfg.method              = 'headshape';
cfg.headshape.headshape = headshape;
cfg.headshape.icp       = 'yes'; % use iterative closest point procedure
cfg.coordsys            = 'neuromag';

mri_aligned_headshape = ft_volumerealign(cfg, mri_aligned_fiducials);

%% check

ft_volumerealign(cfg, mri_aligned_headshape);

%% reslice

cfg            = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri_aligned_headshape);

save(fullfile(data_path, 'mri_resliced.mat'), 'mri_resliced')
%% segment

cfg        = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, mri_resliced);

%% mesh

cfg             = [];
cfg.method      = 'projectmesh';
cfg.tissue      = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, mri_segmented);
mesh_brain = ft_convert_units(mesh_brain, 'm'); % Use SI Units


cfg             = [];
cfg.method      = 'projectmesh';
cfg.tissue      = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
mesh_scalp = ft_convert_units(mesh_scalp, 'm'); % Use SI Units

%% headmodel

cfg = [];
cfg.method = 'singleshell';

headmodel = ft_prepare_headmodel(cfg, mesh_brain);

save(fullfile(data_path, 'headmodel.mat'), 'headmodel')

%% sensors 

sens = ft_read_sens(fullfile(data_path, 'tactile_jitter_avg.fif'), ...
                    'senstype', 'meg');
                
sens = ft_convert_units(sens, 'm');                
                
%% plot

figure
hold on
% ft_plot_mesh(mesh_scalp)
ft_plot_sens(sens)
ft_plot_mesh(mesh_brain)

%% sourcemodel

cfg             = [];
cfg.headmodel   = headmodel; % used to estimate extent of grid
cfg.resolution  = 0.01; % a source per 0.01 m -> 1 cm

sourcemodel = ft_prepare_sourcemodel(cfg);

save(fullfile(data_path, 'sourcemodel'), 'sourcemodel')

%% leadfield

load(fullfile(data_path, 'headmodel.mat'))
load(fullfile(data_path, 'sourcemodel.mat'))

cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.grad        = sens; %% where are the sensors?
cfg.channel = 'MEGGRAD';
cfg.normalize = 'yes';

% how do sources and sensors connect?
sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);

save(fullfile(data_path, 'sourcemodel_and_leadfield'), 'sourcemodel_and_leadfield')

%% warped sourcemodel 

load(['/home/lau/matlab/fieldtrip/template/' ...
    'sourcemodel/standard_sourcemodel3d10mm']);
load(fullfile(data_path, 'mri_resliced.mat'))
template_grid = sourcemodel;
clear sourcemodel;

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_resliced;
cfg.unit      = 'm';

sourcemodel_warp = ft_prepare_sourcemodel(cfg);

save(fullfile(data_path, 'sourcemodel_warp.mat'), 'sourcemodel_warp');

%% add atlas to it

atlas = ft_read_atlas('~/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas, 'm');


% and call ft_sourceinterpolate:
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel_warp_atlas = ft_sourceinterpolate(cfg, atlas, sourcemodel_warp);
save(fullfile(data_path, 'sourcemodel_warp_atlas'), 'sourcemodel_warp_atlas')

%% warped leadfield

cfg = [];
cfg.sourcemodel = sourcemodel_warp;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.grad        = sens; %% where are the sensors?
cfg.channel = 'MEGGRAD';
cfg.normalize = 'yes';

% how do sources and sensors connect?
warped_leadfield = ft_prepare_leadfield(cfg);

save(fullfile(data_path, 'warped_leadfield'), 'warped_leadfield')