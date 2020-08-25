%%  CLEAR AND ADDPATHS

clear variables
restoredefaultpath

sd = 'NatMEG_0245/190809/';
% sd = 'NatMEG_0177/190619/';

raw_path = ['/archive/20068_tactile_prediction/MEG/' sd];
data_path = ['/home/lau/analyses/tactile_prediction_jitter/data/' ...
            sd 'global_ISI'];
fieldtrip_path = '/home/lau/matlab/fieldtrip';     
figures_path = ['/home/lau/analyses/tactile_prediction_jitter/figures/' ...
                sd(1:11) '/global_ISI'];


addpath(fieldtrip_path);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)        

%% DEFINE TRIALS
        
filenames = {'tactile_jitter_raw_tsss_mc.fif'
             'tactile_jitter_raw_tsss_mc-1.fif'
             'tactile_jitter_raw_tsss_mc-2.fif'
             'tactile_jitter_raw_tsss_mc-3.fif'
            };        
        
n_files = length(filenames);
split_files = cell(1, n_files);

for file_index = 1:n_files

    filename = filenames{file_index};
    
    cfg = [];
    cfg.dataset = fullfile(raw_path, filename);
    cfg.trialdef.prestim = 0.159;
    cfg.trialdef.poststim = 0.641;
    cfg.trialdef.eventvalue = {19 295 551 1063, [328:332 584:588 1096:1100]};
    cfg.trialdef.eventtype = 'STI101';

    cfg = ft_definetrial(cfg);
    
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf 0];
    
    data = ft_preprocessing(cfg);
    
    cfg = [];
    cfg.offset = -41; % samples
   
    data = ft_redefinetrial(cfg, data);

    cfg = [];
    cfg.resamplefs = 250; % Hz

    data = ft_resampledata(cfg, data);
    
    split_files{file_index} = data;
    
end

cfg = [];

preprocessed_data_tl = ft_appenddata(cfg, split_files{:});

non_stim_triggers = [328:332 584:588 1096:1100];

for trigger = non_stim_triggers
    
    for trial_index = 1:length(preprocessed_data_tl.trialinfo)
        
        trial = preprocessed_data_tl.trialinfo(trial_index);
        if trigger == trial
            preprocessed_data_tl.trialinfo(trial_index) = 2048;
        end
    end
end

save(fullfile(data_path, 'preprocessed_data_tl'), 'preprocessed_data_tl');

%% filter

load(fullfile(data_path, 'preprocessed_data_tl.mat'))

cfg = [];
cfg.lpfreq = 30;
cfg.lpfilter = 'yes';

filtered_data = ft_preprocessing(cfg, preprocessed_data_tl);

%% timelock

events = [19 295 551 1063 2048];
n_events = length(events);

tl = cell(1, n_events);
tl_filtered = cell(1, n_events);

for event_index = 1:n_events
    
    event = events(event_index);
    
    cfg = [];
    cfg.trials = find(preprocessed_data_tl.trialinfo == event);
    
    tl{event_index} = ft_timelockanalysis(cfg, preprocessed_data_tl);
    tl_filtered{event_index} = ft_timelockanalysis(cfg, filtered_data);
   
    cfg = [];
    
    tl{event_index} = ft_combineplanar(cfg, tl{event_index});
    tl_filtered{event_index} = ft_combineplanar(cfg, tl_filtered{event_index});
    
end

%% plot

cfg = [];
cfg.layout = 'neuromag306mag.lay';

figure
ft_multiplotER(cfg, tl{:});
figure
ft_multiplotER(cfg, tl_filtered{:});

%% POWER ANALYSIS

events = [19 295 551 1063 2048];
n_events = length(events);

psd = cell(1, n_events);
psd_filtered = cell(1, n_events);

for event_index = 1:n_events
    
    event = events(event_index);
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.foi = 3:100;
    
    psd{event_index} = ft_freqanalysis(cfg, tl{event_index});
    psd_filtered{event_index} = ft_freqanalysis(cfg, tl_filtered{event_index});
    
end

%% PLOT POWER

figure

for event_index = 2:4%n_events
    
    cfg = [];
    cfg.channel = 'MEG0422+0423';
    
    this_data = ft_selectdata(cfg, psd{event_index});
    this_data_filtered = ft_selectdata(cfg, psd_filtered{event_index});
    
    subplot(1, 2, 1); hold on
    plot(this_data.freq, this_data.powspctrm);
    
    subplot(1, 2, 2); hold on
    plot(this_data_filtered.freq, this_data_filtered.powspctrm);
    
end