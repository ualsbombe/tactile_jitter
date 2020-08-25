function [dummy] = tfr(input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
subject = input.subject;
date = input.date;
input_file = input.input_file;
channel_sets = input.channel_sets;
clean_based_on = input.clean_based_on;
events = input.events;
overwrite = input.overwrite;

%% LOAD DATA

full_date = [date '_000000'];
subject_path = fullfile(load_path, subject, full_date);
disp(['Loading: ' input_file ' for subject: ' subject]);
data = load(fullfile(subject_path, input_file));

%% CLEAN BASED ON INDICES (make this a function in itself?)

cleaned_data = clean_based_on_indices(subject_path, clean_based_on, data);

%% DO TFR
%#ok<*AGROW>
n_events = length(events);
tfrs = struct();
 
% find indices
for event_index = 1:n_events
    event = events{event_index};
    event_indices = [];
 
    for trigger = event
        event_indices = [event_indices 
                         find(cleaned_data.trialinfo == trigger)]; 
    end
 
    field_name = ['event_' num2str(event(1))];
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = channel_sets;
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = 1:1:40;
    cfg.t_ftimwin = 5 ./ cfg.foi;
    cfg.toi = -1.500:0.020:1.500;
    cfg.trials = event_indices;
    cfg.pad = 'nextpow2';

    tfrs.(field_name) = ft_freqanalysis(cfg, cleaned_data);
     
end
%% SAVE
 
filename = 'tfrs.mat';
fullpath = fullfile(subject_path, filename);
 
save_struct(tfrs, fullpath, overwrite);