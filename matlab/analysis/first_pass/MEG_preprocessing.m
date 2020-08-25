%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath


subject_index = 2;

subjects = {'0007' '0008'};
dates = {'20200120_000000' '20200120_000000'};
subject = subjects{subject_index};
date = dates{subject_index};
home_dir = '/home/lau';

project_name = 'cerebellar_clock';

raw_path = fullfile(home_dir, 'projects', project_name, 'raw', ...
                    subject, date, 'MEG', ...
                    '001.tactile_jitter_raw', 'files');
data_path = fullfile(home_dir, 'projects', project_name, ...
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
%     cfg.trialdef.eventvalue = {18 92 102 112:2:140}; % 0001 0004 0005
    cfg.trialdef.eventvalue = {18};
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
% 
% % non_stim_triggers = 2001:2005; % 0002 0003
% non_stim_triggers = 112:2:140; % 0001 0004 0005
% 
% for trigger = non_stim_triggers
%     
%     for trial_index = 1:length(preprocessed_data_tfr.trialinfo)
%         
%         trial = preprocessed_data_tfr.trialinfo(trial_index);
%         if trigger == trial
%             preprocessed_data_tfr.trialinfo(trial_index) = 256;
%         end
%     end
% end
    
% save(fullfile(data_path, 'preprocessed_data_tfr'), 'preprocessed_data_tfr')

%% TFR

% load(fullfile(data_path, 'cleaned_data.mat'))

events = {18};
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
    cfg.foilim = [1 40];
%     cfg.t_ftimwin = ones(size(cfg.foi)) * 0.500;
    cfg.trials = find(preprocessed_data_tfr.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    tfrs{event_index} = ft_freqanalysis(cfg, preprocessed_data_tfr);
    
end

% save(fullfile(data_path, 'tfrs'), 'tfrs');

%% COMBINE TFRS

% load(fullfile(data_path, 'tfrs'))

n_events = length(tfrs);
tfr_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tfr_cmbs{event_index} = ft_combineplanar(cfg, tfrs{event_index});
    
end

%% PLOT TFR

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.xlim = [-1.200 1.200];
% cfg.zlim = [-0.250 0.250];
cfg.zlim = [0.700 1.300];
% cfg.ylim = [2 40];
% cfg.ylim = [40 100];
% cfg.channel = {'all' '-MEG1012+1013'};
cfg.baselinetype = 'relative';
cfg.baseline = 'yes';

n_contrasts = length(tfr_cmbs);

for index = 1:n_contrasts
    figure
    ft_multiplotTFR(cfg, tfr_cmbs{index});
end