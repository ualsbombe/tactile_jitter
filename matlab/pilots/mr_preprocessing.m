%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

subject_index = 4;

[~, user] = system('uname -n');
if strcmp(user(1:(end-1)), 'lau-LIFEBOOK-E744')
    home_dir = '/home/lau/';
elseif strcmp(user, 'hyades02')
    home_dir = '/users/lau/';
else
    disp('"User" has not been specified yet')
end

subjects = {'0001' '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' '20190930_000000' '20190930_000000' ... 
         '20191008_000000' '20191008_000000'};
dates_MR = {'20180515_132117' '20191015_104445' '20191015_093141' ...
            '20191015_112257' '20191015_121553'};
sequences_MR = {'006' '005' '006' '007' '005'};
subject = subjects{subject_index};
date = dates{subject_index};
mr_date = dates_MR{subject_index};
mr_sequence = [sequences_MR{subject_index} '.t1_mprage_3D_sag_fatsat'];

project_name = 'MINDLAB2019_MEG-CerebellarClock';
raw_path = fullfile(home_dir, 'mounts',  'hyades', 'projects', ...
                    project_name, 'raw', subject, date, 'MEG', ...
                    '001.tactile_jitter_raw', 'files');
data_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                      project_name, 'scratch', 'tactile_jitter', ...
                      'MEG', subject, date);                
fieldtrip_path = '/home/lau/matlab/fieldtrip';
               
mr_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                    project_name,  'raw', subject, mr_date, 'MR', ...
                   mr_sequence, 'files');

addpath(fieldtrip_path)
ft_defaults

%% READ MRI

mr_file_tag = ['PROJ0427_SUBJ' subject '_SER' sequences_MR{subject_index}];
files = dir(fullfile(mr_path, [mr_file_tag '*']));
first_file = files(1).name;

mri = ft_read_mri(fullfile(mr_path, first_file));

%% PLOT MRI

cfg = [];

ft_sourceplot(cfg, mri);

%% ALIGN WITH FIDUCIALS

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';

mri_aligned_fiducials = ft_volumerealign(cfg, mri);

%% ALIGN WITH HEAD POINTS

headshape = ft_read_headshape(fullfile(raw_path, ...
                                    'tactile_jitter_raw.fif'));
                                
cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';

mri_aligned_headshape = ft_volumerealign(cfg, mri_aligned_fiducials);

%% CHECK ALIGNMENT

ft_volumerealign(cfg, mri_aligned_headshape);

%% RESLICE MRI

cfg = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri_aligned_headshape);

save(fullfile(data_path, 'mri_resliced'), 'mri_resliced')

%% SEGMENT BRAIN

cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, mri_resliced);

%% CREATE MESHES

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, mri_segmented);
mesh_brain = ft_convert_units(mesh_brain, 'm');

cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
mesh_scalp = ft_convert_units(mesh_scalp, 'm');

save(fullfile(data_path, 'mesh_scalp'), 'mesh_scalp')

%% CREATE HEAD MODEL

cfg = [];
cfg.method = 'singleshell';

headmodel = ft_prepare_headmodel(cfg, mesh_brain);

save(fullfile(data_path, 'headmodel'), 'headmodel')

%% CHECK CO-REGISTRATION

sens = ft_read_sens(fullfile(raw_path, 'tactile_jitter_raw.fif'), ...
                    'senstype', 'meg');
sens = ft_convert_units(sens, 'm');                

figure
hold on
ft_plot_mesh(mesh_scalp)
ft_plot_sens(sens)

%% WARPED SOURCE MODEL

load(fullfile(fieldtrip_path, 'template', 'sourcemodel', ...
              'standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel

cfg = [];
cfg.warpmni = 'yes';
cfg.template = template_grid;
cfg.nonlinear = 'yes';
cfg.mri = mri_resliced;
cfg.unit = 'm';

sourcemodel_warp = ft_prepare_sourcemodel(cfg);

save(fullfile(data_path, 'sourcemodel_warp'), 'sourcemodel_warp');

%% CHECK SOURCE MODEL

figure
hold on
ft_plot_mesh(mesh_scalp)
ft_plot_sens(sens)
ft_plot_mesh(sourcemodel_warp)