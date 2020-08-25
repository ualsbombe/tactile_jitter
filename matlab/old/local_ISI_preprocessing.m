%%  CLEAR AND ADDPATHS

clear variables
restoredefaultpath

% sd = 'NatMEG_0245/190809/';
% sd = 'NatMEG_0177/190619/';
sd = 'NatMEG_0521/190808/';

raw_path = ['/archive/20068_tactile_prediction/MEG/' sd];
data_path = ['/home/lau/analyses/tactile_prediction_jitter/data/' ...
            sd 'local_ISI'];
fieldtrip_path = '/home/lau/matlab/fieldtrip';     
figures_path = ['/home/lau/analyses/tactile_prediction_jitter/figures/' ...
                sd(1:11) '/local_ISI'];


addpath(fieldtrip_path);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)

%% DEFINE TRIALS AND PREPROCESS

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
    cfg.trialfun = 'local_ISI_trialfun';
    cfg.trialdef.first_value  = 21;
    cfg.trialdef.second_value = 22;
    cfg.trialdef.new_values = [100 200 300];
    cfg.trialdef.non_stim_values = [328:332 584:588 1096:1100];
    cfg.trialdef.non_stim_type   = 'STI101';
    cfg.trialdef.prestim = 1.159;
    cfg.trialdef.poststim = 1.241;
    
    cfg = ft_definetrial(cfg);
    
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf Inf];    
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

preprocessed_data_tfr = ft_appenddata(cfg, split_files{:});

non_stim_triggers = [328:332 584:588 1096:1100];

for trigger = non_stim_triggers
    
    for trial_index = 1:length(preprocessed_data_tfr.trialinfo)
        
        trial = preprocessed_data_tfr.trialinfo(trial_index);
        if trigger == trial
            preprocessed_data_tfr.trialinfo(trial_index) = 2048;
        end
    end
end

save(fullfile(data_path, 'preprocessed_data_tfr'), 'preprocessed_data_tfr');

%% CLEAN DATA

load(fullfile(data_path, 'preprocessed_data_tfr.mat'))

cfg = [];
cfg.method = 'summary';
cfg.channel = 'MEGGRAD';

cleaned_data = ft_rejectvisual(cfg, preprocessed_data_tfr);

save(fullfile(data_path, 'cleaned_data'), 'cleaned_data');

%%  TFR

% load(fullfile(data_path, 'preprocessed_data_tfr.mat'))
load(fullfile(data_path, 'cleaned_data.mat'))

events = {100 200 300 2048};
n_events = length(events);

tfrs = cell(1, n_events);

for event_index = 1:n_events

    event = events{event_index};
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEGGRAD';
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.toi = -1.200:0.020:1.200;
    cfg.foi = 2:2:40;
    cfg.t_ftimwin = ones(size(cfg.foi)) * 0.500;
    cfg.trials = find(cleaned_data.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    tfrs{event_index} = ft_freqanalysis(cfg, cleaned_data);
    
end

save(fullfile(data_path, 'tfrs'), 'tfrs');

%% COMBINE TFRS

load(fullfile(data_path, 'tfrs.mat'))

n_events = 4;
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

save(fullfile(data_path, 'tfr_cmb_contrasts.mat'), 'tfr_cmb_contrasts')

%% TFR CONTRAST LAYER TWO

contrasts = {[1 2] [1 3] [2 3]};
n_contrasts = length(contrasts);

tfr_contrasts_l2 = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    contrast = contrasts{contrast_index};
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1 - x2';
    
    tfr_contrasts_l2{contrast_index} = ft_math(cfg, ...
                                           tfr_cmb_contrasts{contrast(1)}, ...
                                           tfr_cmb_contrasts{contrast(2)});
                                          
                                 
end


%% PLOT

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.layout = 'neuromag306mag.lay';
% cfg.baseline = [-0.400 Inf];
% cfg.baselinetype = 'relative';
% cfg.xlim = [-0.400 Inf];
cfg.zlim = [-0.250 0.250];

% ft_multiplotTFR(cfg, tfr_cmbs{1});
% ft_multiplotTFR(cfg, tfrs{3});
figure
ft_multiplotTFR(cfg, tfr_cmb_contrasts{1});
figure
ft_multiplotTFR(cfg, tfr_cmb_contrasts{2});
figure
ft_multiplotTFR(cfg, tfr_cmb_contrasts{3});

% figure
% ft_multiplotTFR(cfg, tfr_contrasts_l2{1});
% figure
% ft_multiplotTFR(cfg, tfr_contrasts_l2{2});
% figure
% ft_multiplotTFR(cfg, tfr_contrasts_l2{3});

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
    title(['Local ISI; Jitter: ' name])
    print(fullfile(figures_path, ['tfr_' cfg.channel '_' name '.jpg']), ...
            '-djpeg', '-r300')
end

%% TEMPORAL SPECTRAL EVOLUTIONN

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [8 12];

alpha_data = ft_preprocessing(cfg, preprocessed_data_tfr);

cfg.bpfreq = [15 30];

beta_data = ft_preprocessing(cfg, preprocessed_data_tfr);

cfg = [];
cfg.parameter = 'trial';
cfg.operation = 'abs';

alpha_data = ft_math(cfg, alpha_data);
beta_data  = ft_math(cfg, beta_data);

cfg = [];

tse_alpha = ft_timelockanalysis(cfg, alpha_data);
tse_beta  = ft_timelockanalysis(cfg, beta_data);

cfg = [];
cfg.layout = 'neuromag306mag.lay';
cfg.xlim = [-0.500 0.750];

figure
ft_multiplotER(cfg, tse_alpha);
figure
ft_multiplotER(cfg, tse_beta);
