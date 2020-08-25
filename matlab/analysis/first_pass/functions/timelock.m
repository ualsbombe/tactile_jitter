function [dummy] = timelock(input)
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
lowpass_filter = input.lowpass_filter;
highpass_filter = input.highpass_filter;
output_file = input.output_file;
overwrite = input.overwrite;

%% LOAD DATA

full_date = [date '_000000'];
subject_path = fullfile(load_path, subject, full_date);
disp(['Loading: ' input_file ' for subject: ' subject]);
data = load(fullfile(subject_path, input_file));

%% CLEAN BASED ON INDICES (make this a function in itself?)

cleaned_data = clean_based_on_indices(subject_path, clean_based_on, data);

%% DO TIMELOCKED
%#ok<*AGROW>

n_events = length(events);
timelockeds = struct();
 
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
    cfg.trials = event_indices;
    cfg.channel = channel_sets;
    if ~isempty(lowpass_filter)
        cfg.preproc.lpfreq = lowpass_filter;
        cfg.preproc.lpfilter = 'yes';
    end
    if ~isempty(highpass_filter)
        cfg.preproc.hpfreq = highpass_filter;
        cfg.preproc.hpfilter = 'yes';
    end

    timelockeds.(field_name) = ft_timelockanalysis(cfg, cleaned_data);
     
end

 %% SAVE
 if isempty(lowpass_filter)
    filename = output_file;
 else
     filename = ['lowpass_' num2str(lowpass_filter) '_Hz_' output_file ];
 end
 
 if isempty(highpass_filter)
    filename = filename; %#ok<ASGSL>
 else
     filename = ['highpass_' num2str(highpass_filter) '_Hz_' ...
         filename];
 end
 fullpath = fullfile(subject_path, filename);
 
 save_struct(timelockeds, fullpath, overwrite);