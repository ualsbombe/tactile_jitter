%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

[~, user] = system('uname -n');
if strcmp(user(1:(end-1)), 'lau-LIFEBOOK-E744')
    home_dir = '/home/lau/';
elseif strcmp(user, 'hyades02')
    home_dir = '/users/lau/';
else
    disp('"User" has not been specified yet')
end

subject_index = 1;

subjects = {'0001' '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' '20190930_000000' '20190930_000000' ...
         '20191008_000000' '20191008_000000'};
subject = subjects{subject_index};
date = dates{subject_index};

project_name = 'MINDLAB2019_MEG-CerebellarClock';
data_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                      project_name, 'scratch', 'tactile_jitter', ...
                      'MEG', subject, date);
figures_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                      project_name, 'scratch', 'tactile_jitter', 'figures', ...
                      subject, date);
fieldtrip_path = '/home/lau/matlab/fieldtrip';     

addpath(fieldtrip_path)
ft_defaults

%% SET TOILIM AND FREQUENCY

% cerebellum error?
toilim = [0.600 0.800]; 
frequency = 15;
tapsmofrq = 4;

% cerebellum and insula preparation?
% toilim = [-0.700 -0.300];
% frequency = 10;
% tapsmofrq = 4;

% insula (L) and somatosensory cortices for regular only
% toilim = [-0.600 -0.400];
% frequency = 60;
% tapsmofrq = 20;

% Insula (R) and cerebellum for irregular
% toilim = [0.600 0.900];
% frequency = 80;
% tapsmofrq = 20;

band_string = [num2str(frequency) '_Hz'];
time_string = [num2str(toilim(1)*1e3) '_to_' num2str(toilim(2)*1e3) '_ms'];

%% CHOOSE DATA

load(fullfile(data_path, 'cleaned_data.mat'));

cfg = [];
cfg.toilim = toilim;

data_t_win = ft_redefinetrial(cfg, cleaned_data);

%% POW AND CSD

events = {18 92 102 256};
n_events = length(events);

powcsds = cell(1, n_events);

for event_index = 1:n_events
    
    event = events{event_index};

    cfg = [];
    cfg.trials = find(data_t_win.trialinfo == event);
    cfg.channel = 'MEGGRAD';
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'powandcsd';
    cfg.keeptrials = 'no';
    cfg.foi = frequency;
    cfg.tapsmofrq = 4;
    cfg.pad = 'nextpow2';

    powcsds{event_index} = ft_freqanalysis(cfg, data_t_win);
    
end
%% WARPED LEADFIELD

load(fullfile(data_path, 'sourcemodel_warp.mat'))
load(fullfile(data_path, 'headmodel.mat'))
grad = cleaned_data.grad;
grad = ft_convert_units(grad, 'm');

cfg = [];
cfg.sourcemodel = sourcemodel_warp;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.grad        = grad; %% where are the sensors?
cfg.channel = 'MEGGRAD';
cfg.normalize = 'yes';

% how do sources and sensors connect?
warped_leadfield = ft_prepare_leadfield(cfg, cleaned_data);

save(fullfile(data_path, 'warped_leadfield'), 'warped_leadfield')
save(fullfile(data_path, 'grad'), 'grad');

%% QUALITY PLOT

figure;
hold on
ft_plot_mesh(sourcemodel_warp, 'vertexcolor', 'red')
ft_plot_headmodel(headmodel)
ft_plot_sens(grad);

%% COMMON FILTER ANALYSIS

load(fullfile(data_path, 'headmodel.mat'))
load(fullfile(data_path, 'warped_leadfield.mat'))
load(fullfile(data_path, 'grad.mat'))

% pow

cfg = [];
cfg.channel = 'MEGGRAD';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'powandcsd';
cfg.keeptrials = 'no';
cfg.foi = frequency;
cfg.tapsmofrq = 4;
cfg.pad = 'nextpow2';

powcsd_all = ft_freqanalysis(cfg, data_t_win);

% common filter

cfg = [];
cfg.method = 'dics';
cfg.sourcemodel = warped_leadfield;
cfg.headmodel = headmodel;
cfg.channel = 'MEGGRAD';
cfg.frequency = frequency;
cfg.grad = grad;
cfg.dics.keepfilter = 'yes';

source_all = ft_sourceanalysis(cfg, powcsd_all);

%% SOURCE ANALYSIS (common filter)         

load(['/home/lau/matlab/fieldtrip/template/' ...
    'sourcemodel/standard_sourcemodel3d10mm']); % loads "sourcemodel"
load(fullfile(data_path, 'grad.mat'))

events = {18 92 102 256};
n_events = length(events);

sources_common = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'dics';
    cfg.sourcemodel = warped_leadfield;
    cfg.headmodel = headmodel;
    cfg.channel = 'MEGGRAD';
    cfg.frequency = frequency;
    cfg.grad = grad;
    cfg.sourcemodel.filter = source_all.avg.filter;

    sources_common{event_index} = ft_sourceanalysis(cfg, ...
                                                    powcsds{event_index});
    sources_common{event_index}.pos = sourcemodel.pos; % template
    
end

%% SOURCE CONTRASTS

contrasts = {[1 4] [2 4] [3 4] [1 2] [1 3] [2 3]};
n_contrasts = 3; %length(contrasts);

source_contrasts = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    contrast = contrasts{contrast_index};
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = '(x1-x2) / (x1+x2)';
    
    source_contrasts{contrast_index} = ft_math(cfg, ...
                                               sources_common{contrast(1)}, ...
                                               sources_common{contrast(2)});
    source_contrasts{contrast_index}.coordsys = 'mni';
                                 
end

savename = fullfile(data_path, ['source_contrasts_' band_string '_' ...
                                time_string]);

save(savename, 'source_contrasts')

%% INTERPOLATE CONTRASTS TO TEMPLATE BRAIN

load('standard_mri.mat') %% loads "mri"

n_contrasts = length(source_contrasts);

source_contrast_ints = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts

    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;

    source_contrast_ints{contrast_index} = ft_sourceinterpolate(cfg, ...
                                           source_contrasts{contrast_index}, ...
                                           mri);    
end

savename = fullfile(data_path, ['source_contrast_ints_' band_string '_' ...
                                time_string]);

save(savename, 'source_contrast_ints')

%% LOAD IF NECESSARY

savename = fullfile(data_path, ['source_contrast_ints_' band_string '_' ...
                                time_string]);
load(savename)

%% AND GET NEG IF NECESSARY

n_contrasts = 3;
source_contrast_ints_neg = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    cfg = [];
    cfg.operation = 'x1 * -1';
    cfg.parameter = 'pow';
    
    source_contrast_ints_neg{contrast_index} = ft_math(cfg, ...
                                        source_contrast_ints{contrast_index});
    source_contrast_ints_neg{contrast_index}.anatomy = ...
                   source_contrast_ints{contrast_index}.anatomy;
    source_contrast_ints_neg{contrast_index}.coordsys = ...
                   source_contrast_ints{contrast_index}.coordsys;
    source_contrast_ints_neg{contrast_index}.dim = ...
                   source_contrast_ints{contrast_index}.dim;
    source_contrast_ints_neg{contrast_index}.transform = ...
                   source_contrast_ints{contrast_index}.transform;
    source_contrast_ints_neg{contrast_index}.unit = ...
                   source_contrast_ints{contrast_index}.unit;
    source_contrast_ints_neg{contrast_index}.inside = ...
                   source_contrast_ints{contrast_index}.inside;
                                            
end

%% PLOT INTERPOLATED CONTRASTS ON TEMPLATE BRAIN

close all
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));
                             
names = {'0%' '5%' '15%'};

for contrast_index = 1:1%3
    name = names{contrast_index};

    cfg = [];
    cfg.funparameter = 'pow';
    cfg.location = [-34 -65 -36]; 
%     cfg.location = [-46 -59 -30]; 
%     cfg.location = [8 -69 -44];
%     cfg.location = [44 -39 -34];
    cfg.funcolorlim = [0.06 0.14];
    cfg.opacitylime = cfg.funcolorlim;
    cfg.atlas = atlas;
%     cfg.roi = {'Cerebellum_10_R'}
    cfg.crosshair = 'no';

    ft_sourceplot(cfg, source_contrast_ints{contrast_index});
%     ft_sourceplot(cfg, source_contrast_ints_neg{contrast_index});

%     title(['Jitter: ' name])
%     print(fullfile(figures_path, ['sourceplot_' name '_' ...
%                                    band_string '_' time_string '.jpg']), ...
%                                     '-djpeg', '-r300');
                                 
end