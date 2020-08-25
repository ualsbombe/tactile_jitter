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

subject_index = 5;

subjects = {'0001' '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' '20190930_000000' '20190930_000000' ...
         '20191008_000000' '20191008_000000'};
subject = subjects{subject_index};
date = dates{subject_index};

project_name = 'MINDLAB2019_MEG-CerebellarClock';
raw_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', project_name, ...
                    'raw', ...
                    subject, date, 'MEG', ...
                    '001.tactile_jitter_raw', 'files');
data_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', project_name, ...
                      'scratch', 'tactile_jitter', 'MEG', subject, date);
figures_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                      project_name, 'scratch', 'tactile_jitter', 'figures', ...
                      subject, date);
fieldtrip_path = '/home/lau/matlab/fieldtrip';     

addpath(fieldtrip_path)
ft_defaults

%% DEFINE TRIALS

filenames = {'tactile_jitter_raw.fif'
             'tactile_jitter_raw-1.fif'
             'tactile_jitter_raw-2.fif'
             'tactile_jitter_raw-3.fif'
            };
        
n_files = length(filenames);
split_files = cell(1, n_files);

for file_index = 1:n_files

    filename = filenames{file_index};
    
    cfg = [];
    cfg.dataset = fullfile(raw_path, filename);
    cfg.trialdef.prestim = 1.500;
    cfg.trialdef.poststim = 1.500;
%     cfg.trialdef.eventvalue = {18 92 102 2001:2005}; % 0002 0003
    cfg.trialdef.eventvalue = {18 92 102 112:2:140}; % 0001 0004 0005
    cfg.trialdef.eventtype = 'STI101';
%     cfg.trialfun = 'ft_trialfun_ns_fix'; % 0002 0003

    cfg = ft_definetrial(cfg);  
    
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf Inf];

    data = ft_preprocessing(cfg);

    cfg = [];
    cfg.resamplefs = 250; % Hz

    data = ft_resampledata(cfg, data);
    
    split_files{file_index} = data;
    
end

cfg = [];

preprocessed_data_tfr = ft_appenddata(cfg, split_files{:});

% non_stim_triggers = 2001:2005; % 0002 0003
non_stim_triggers = 112:2:140; % 0001 0004 0005

for trigger = non_stim_triggers
    
    for trial_index = 1:length(preprocessed_data_tfr.trialinfo)
        
        trial = preprocessed_data_tfr.trialinfo(trial_index);
        if trigger == trial
            preprocessed_data_tfr.trialinfo(trial_index) = 256;
        end
    end
end
    
save(fullfile(data_path, 'preprocessed_data_tfr'), 'preprocessed_data_tfr')

%% SELECT DATA %% bug 1152

cfg = [];
cfg.channel = 'MEGGRAD';

data_select = ft_selectdata(cfg, preprocessed_data_tfr);

%% CLEAN DATA 

cfg = [];
cfg.method = 'summary';
% cfg.channel = 'MEGGRAD';

cleaned_data = ft_rejectvisual(cfg, data_select);

save(fullfile(data_path, 'cleaned_data'), 'cleaned_data');

%% TFR

load(fullfile(data_path, 'cleaned_data.mat'))

events = {18 92 102 256};
n_events = length(events);

tfrs = cell(1, n_events);

for event_index = 1:n_events

    event = events{event_index};
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
%     cfg.method = 'mtmconvol';
    cfg.method = 'wavelet';
    cfg.width = 7;
%     cfg.taper = 'hanning';
    cfg.toi = -1.500:0.020:1.500;
%     cfg.foi = 2:2:40;
    cfg.foilim = [1 100];
%     cfg.t_ftimwin = ones(size(cfg.foi)) * 0.500;
    cfg.trials = find(cleaned_data.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    tfrs{event_index} = ft_freqanalysis(cfg, cleaned_data);
    
end

save(fullfile(data_path, 'tfrs'), 'tfrs');

%% COMBINE TFRS

load(fullfile(data_path, 'tfrs'))

n_events = length(tfrs);
tfr_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tfr_cmbs{event_index} = ft_combineplanar(cfg, tfrs{event_index});
    
end

%% TFR CONTRASTS

contrasts = {[1 4] [2 4] [3 4]};
n_contrasts = length(contrasts);

tfr_contrasts = cell(1, n_contrasts);
tfr_cmb_contrasts = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    contrast = contrasts{contrast_index};
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2) / (x1+x2)';
    
    tfr_contrasts{contrast_index} = ft_math(cfg, ...
                                               tfrs{contrast(1)}, ...
                                               tfrs{contrast(2)});
                                           
    cfg = [];
    
    tfr_cmb_contrasts{contrast_index} = ft_combineplanar(cfg, ...
                                            tfr_contrasts{contrast_index});
                                 
end

save(fullfile(data_path, 'tfr_cmb_contrasts'), 'tfr_cmb_contrasts');

%% LOAD IF NECESSARY
load(fullfile(data_path, 'tfr_cmb_contrasts.mat'))

%% PLOT TFR

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.xlim = [-1.200 1.200];
% cfg.zlim = [-0.250 0.250];
cfg.zlim = [-0.250 0.250];
% cfg.ylim = [2 40];
cfg.ylim = [40 100];
% cfg.channel = {'all' '-MEG1012+1013'};

n_contrasts = length(tfr_cmb_contrasts);

for index = 1:n_contrasts
    figure
    ft_multiplotTFR(cfg, tfr_cmb_contrasts{index});
end

%% DEFAULTS

set(0, 'defaultaxesfontweight', 'bold', 'defaultaxesfontsize', 14)

%% PLOT SINGLE PLOTS

close all

load(fullfile(data_path, 'tfr_cmb_contrasts.mat'))
% channels = {'MEG2142+2143' 'MEG2132+2133'};
channels = {'MEG0412+0413' 'MEG1322+1323'};
n_channels = length(channels);
n_contrasts = length(tfr_cmb_contrasts);
names = {'0%' '5%' '15%'};

for contrast_index = 1:n_contrasts
    
    name = names{contrast_index};
    h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);

    for channel_index = 1:n_channels

        subplot(1, n_channels, channel_index);

        cfg = [];
        cfg.channel = channels{channel_index};
%         cfg.zlim = [-0.300 0.200];
        cfg.ylim = [40 100];

        ft_singleplotTFR(cfg, tfr_cmb_contrasts{contrast_index});

    end
        print(fullfile(figures_path, ['tfr_plot_' name '_' [channels{:}] ...
                                    '.jpg']), ...
                                    '-djpeg', '-r300');
   
end  

%% TOPO PLOT 3D

load(fullfile(data_path, 'headmodel.mat'))
load(fullfile(data_path, 'tfr_cmb_contrasts.mat'))

%%
close all

cfg = [];
cfg.frequency = 15;
cfg.latency = [0.600 0.800];

temp = ft_selectdata(cfg, tfr_cmb_contrasts{1});
temp.grad = ft_convert_units(temp.grad, 'm');

temp.powspctrm = squeeze(mean(temp.powspctrm, 3));

figure; hold on
view(0, 0)
ft_plot_topo3d(temp.grad.chanposold(1:3:306, :), temp.powspctrm(1:102), ...
                'facealpha', 0.9);
ft_plot_headmodel(headmodel)
c = colorbar;

%% PLOT TOPO PLOTS

close all

load(fullfile(data_path, 'tfr_cmb_contrasts.mat'))
n_contrasts = length(tfr_cmb_contrasts);
names = {'0%' '5%' '15%'};

for contrast_index = 1:n_contrasts
    
    name = names{contrast_index};
    figure
    cfg = [];
    cfg.ylim = [80 80];
    cfg.xlim = [0.600 0.900];
    cfg.zlim = [-0.150 0.150];
    cfg.colorbar = 'yes';

    ft_topoplotTFR(cfg, tfr_cmb_contrasts{contrast_index});

%     print(fullfile(figures_path, ['topo_plot_' name ...
%                                 '.jpg']), ...
%                                 '-djpeg', '-r300');
   
end  

%% TFR NON CONTRASTED

close all

n_events = length(tfr_cmbs);

for event_index = 1:n_events
    
    cfg = [];
    cfg.baseline = 'yes';
    cfg.baselinetype = 'relative';
    cfg.layout = 'neuromag306cmb.lay';
%     cfg.zlim = [0.600 1.400];
    cfg.ylim = [2 40];
%     cfg.ylim = [40 100];
    
    figure
    ft_multiplotTFR(cfg, tfr_cmbs{event_index});
    
end