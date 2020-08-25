function [dummy] = beamformer(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
fieldtrip_path = input.fieldtrip_path;
subject = input.subject;
meg_date = input.meg_date;
mr_date = input.mr_date;
input_files = input.input_files;
template = input.template;
clean_based_on = input.clean_based_on;
frequency = input.frequency;
toilim = input.toilim;
smoothing = input.smoothing;
output_file = input.output_file;
overwrite = input.overwrite;

%% LOAD FILES
full_meg_date = [meg_date '_000000'];

mr_path = fullfile(load_path, subject, mr_date);
meg_path = fullfile(load_path, subject, full_meg_date);

meg_data = load(fullfile(meg_path, input_files{1}));
headmodel = load(fullfile(mr_path, input_files{2}));
warped_leadfield = load(fullfile(meg_path, input_files{3}));
% load sourcemodel
load(fullfile(fieldtrip_path, 'template', 'sourcemodel', template));

%% CLEAN

cleaned_data = clean_based_on_indices(meg_path, clean_based_on, meg_data);
grad = cleaned_data.grad;
grad = ft_convert_units(grad, 'm');

%% CUT IN TIME

cfg = [];
cfg.toilim = toilim;

data_t_win = ft_redefinetrial(cfg, cleaned_data);
clear cleaned_data meg_data

%% ESTIMATE POW AND CSD

events = {1 3}; %% change this
n_events = length(events);

powcsds = cell(1, n_events);

for event_index = 1:n_events
    
    event = events{event_index};

    cfg = [];
    cfg.trials = find(data_t_win.trialinfo == event);
    cfg.channel = 'MEGGRAD'; %% change this
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'powandcsd';
    cfg.keeptrials = 'no';
    cfg.foi = frequency; 
    cfg.tapsmofrq = smoothing; %
    cfg.pad = 'nextpow2';

    powcsds{event_index} = ft_freqanalysis(cfg, data_t_win);
    
end

%% SOURCE ANALYSIS

events = {1 3}; %% change this
n_events = length(events);

beamformers = struct();

for event_index = 1:n_events
    
    event = events{event_index};
    event_indices = [];
 
    for trigger = event
        event_indices = [event_indices 
                    find(data_t_win.trialinfo == trigger)];  %#ok<AGROW>
    end
 
    field_name = ['event_' num2str(event(1))];

    cfg = [];
    cfg.method = 'dics';
    cfg.sourcemodel = warped_leadfield;
    cfg.headmodel = headmodel;
    cfg.channel = 'MEGGRAD';
    cfg.frequency = frequency;
    cfg.grad = grad;

    beamformers.(field_name) = ft_sourceanalysis(cfg, ...
                                                    powcsds{event_index});
    beamformers.(field_name).pos = sourcemodel.pos; % template
    
end

%% SAVE

fullpath = fullfile(meg_path, output_file);

save_struct(beamformers, fullpath, overwrite);