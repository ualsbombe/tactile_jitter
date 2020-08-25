%%  CLEAR AND ADDPATHS (KI)

clear variables
restoredefaultpath

% sd = 'NatMEG_0245/190809/';
sd = 'NatMEG_0177/190619/';
% sd = 'NatMEG_0521/190808/';

mount_path = '/home/lau/mounts/compute';
project_path = '/home/lau/analyses/tactile_prediction_jitter';

raw_path = fullfile(mount_path, 'archive', '20068_tactile_prediction', ...
                    'MEG', sd);
data_path = fullfile(mount_path, project_path, 'data', sd, 'global_ISI');
       
fieldtrip_path = '/home/lau/matlab/fieldtrip';     
figures_path = fullfile(mount_path, project_path,  'figures', ...
                sd(1:11), 'global_ISI');
            
mr_path = fullfile(data_path, '..');
            
addpath(fieldtrip_path);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)        

%% DEFINE TRIALS

% filenames = {'tactile_jitter_raw_tsss_mc.fif'
%              'tactile_jitter_raw_tsss_mc-1.fif'
%              'tactile_jitter_raw_tsss_mc-2.fif'
%              'tactile_jitter_raw_tsss_mc-3.fif'
%             };
        
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
    cfg.trialdef.prestim = 1.159;
    cfg.trialdef.poststim = 1.241;
%     cfg.trialdef.eventvalue = {295 551 1063, [328:332 584:588 1096:1100]};
    cfg.trialdef.eventvalue = {18 112:2:140};
    cfg.trialdef.eventtype = 'STI101';

    cfg = ft_definetrial(cfg);
    
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf Inf];
%     cfg.dftfilter = 'yes';
%     cfg.dftfreq = 50:50:300;
%     cfg.lpfilter = 'yes';
%     cfg.lpfreq = 70;
    
    data = ft_preprocessing(cfg);
    
%     cfg = [];
%     cfg.offset = -41; % samples
%    
%     data = ft_redefinetrial(cfg, data);

    cfg = [];
    cfg.resamplefs = 250; % Hz

    data = ft_resampledata(cfg, data);
    
    split_files{file_index} = data;
    
end

cfg = [];

preprocessed_data_tfr = ft_appenddata(cfg, split_files{:});

% non_stim_triggers = [328:332 584:588 1096:1100];
non_stim_triggers = 112:2:140;

for trigger = non_stim_triggers
    
    for trial_index = 1:length(preprocessed_data_tfr.trialinfo)
        
        trial = preprocessed_data_tfr.trialinfo(trial_index);
        if trigger == trial
            preprocessed_data_tfr.trialinfo(trial_index) = 2048;
        end
    end
end

save(fullfile(data_path, 'preprocessed_data_tfr'), 'preprocessed_data_tfr');


%% CLEANED DATA

load(fullfile(data_path, 'preprocessed_data_tfr.mat'))

cfg = [];
cfg.method = 'summary';
cfg.channel = 'MEGGRAD';

cleaned_data = ft_rejectvisual(cfg, preprocessed_data_tfr);

save(fullfile(data_path, 'cleaned_data'), 'cleaned_data');

%% TIMELOCKED

load(fullfile(data_path, 'cleaned_data.mat'))

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;

filtered_data = ft_preprocessing(cfg, cleaned_data);
% events = {295 551 1063};
events = {18 2048};
n_events = length(events);

timelockeds = cell(1, n_events);

for event_index = 1:n_events
    
    event = events{event_index};
    
    cfg = [];
    cfg.trials = find(filtered_data.trialinfo == event);
    cfg.preproc.demean = 'yes';
    cfg.preproc.baselinewindow = [-0.200 0];
    cfg.latency = [-0.200 0.750];
    timelockeds{event_index} = ft_timelockanalysis(cfg, filtered_data);
    
end

%% COMBINE PLANARS

n_events = 2;
tl_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tl_cmbs{event_index} = ft_combineplanar(cfg, timelockeds{event_index});
    
end

%% PLOT

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.xlim = [-0.200 0.750];

% ft_multiplotER(cfg, timelockeds{:});
ft_multiplotER(cfg, tl_cmbs{:});

%% POWER ANALYSIS

load(fullfile(data_path, 'cleaned_data.mat'))

events = {295 551 1063};
n_events = length(events);
psds = cell(1, n_events);

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;

preprocessed_data = ft_preprocessing(cfg, cleaned_data);

for event_index = 1:n_events

    event = events{event_index};

    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foi = 3:100;
    cfg.trials = find(preprocessed_data.trialinfo == event);

    psds{event_index} = ft_freqanalysis(cfg, preprocessed_data);
    
end

%% PLOT CHANNEL

figure
hold on

for i = 1:3

    channel = 'MEG0422';

    cfg = [];
    cfg.channel = channel;

    this_data = ft_selectdata(cfg, psds{i});


plot(this_data.freq, this_data.powspctrm)

end

%% TFR

% load(fullfile(data_path, 'preprocessed_data_tfr.mat'))
load(fullfile(data_path, 'cleaned_data.mat'))

% events = {295 551 1063 2048};
events = {18 2048};
n_events = length(events);

tfrs = cell(1, n_events);

for event_index = 1:n_events

    event = events{event_index};
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.toi = -1.200:0.020:1.200;
    cfg.foi = 2:2:100;
    cfg.t_ftimwin = ones(size(cfg.foi)) * 0.500;
    cfg.trials = find(cleaned_data.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    tfrs{event_index} = ft_freqanalysis(cfg, cleaned_data);
    
end

save(fullfile(data_path, 'tfrs'), 'tfrs');

%% COMBINE TFRS

load(fullfile(data_path, 'tfrs'))

n_events = 4;
n_events = 2;
tfr_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tfr_cmbs{event_index} = ft_combineplanar(cfg, tfrs{event_index});
    
end

%% TFR CONTRASTS

% contrasts = {[1 4] [2 4] [3 4]};
contrasts = {[1 2]};
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


%% PLOT

% load(fullfile(data_path, 'tfr_cmb_contrasts.mat'))

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.layout = 'neuromag306mag.lay';
% cfg.baseline = [-0.400 Inf];
% cfg.baselinetype = 'relative';
cfg.xlim = [-1.000 1.000];
cfg.ylim = [2 40];
cfg.zlim = [-0.1500 0.1500];
% % ft_multiplotTFR(cfg, tfr_cmbs{1});
% ft_multiplotTFR(cfg, tfrs{3});

figure
ft_multiplotTFR(cfg, tfr_cmb_contrasts{1});
% figure
% ft_multiplotTFR(cfg, tfr_cmb_contrasts{2});
% figure
% ft_multiplotTFR(cfg, tfr_cmb_contrasts{3});

%% LOAD HEADMODEL

load(fullfile(mr_path, 'headmodel.mat'))

%% TOPO PLOT 3D

cfg = [];
cfg.frequency = 15;
cfg.latency = [-0.100 0.100];

temp = ft_selectdata(cfg, tfr_cmb_contrasts{1});
temp.grad = ft_convert_units(temp.grad, 'm');

temp.powspctrm = squeeze(mean(temp.powspctrm, 3));

figure; hold on
ft_plot_topo3d(temp.grad.chanposold(1:3:306, :), temp.powspctrm(1:102), ...
                'facealpha', 0.9);
ft_plot_headmodel(headmodel)

%% SINGLE PLOT

names = {'0%' '5%' '15%'};

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.channel = 'MEG0422+0423';
cfg.zlim = [-0.750 0.750];
cfg.colorbar = 'yes';
cfg.colorbartext = 'Power relative to non-stimulation';
cfg.title = ' ';

for contrast_index = 1:3

    name = names{contrast_index};
    figure;
    ft_singleplotTFR(cfg, tfr_cmb_contrasts{contrast_index});
    xlabel('Time (s)')
    ylabel('Frequency (Hz')
    title(['Global ISI; Jitter: ' name])
    print(fullfile(figures_path, ['tfr_' cfg.channel '_' name '.jpg']), ...
            '-djpeg', '-r300')
end