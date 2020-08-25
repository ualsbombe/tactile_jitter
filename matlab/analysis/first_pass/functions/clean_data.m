function [dummy] = clean_data(input)
%CLEAN DATA Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';
load_path = input.load_path;
subject = input.subject;
date = input.date;
input_file = input.input_file;
channel_sets = input.channel_sets;
rejected_trials_file = input.rejected_trials_file;
output_type = input.output_type;
overwrite = input.overwrite;


%% LOAD DATA

full_date = [date '_000000'];
subject_path = fullfile(load_path, subject, full_date);
disp(['Loading: ' input_file ' for subject: ' subject]);
data = load(fullfile(subject_path, input_file));

%% LOOP THROUGH CHANNELS

n_channel_sets = length(channel_sets);

for channel_set_index = 1:n_channel_sets
    
    channel_set = channel_sets{channel_set_index};
    
    cfg = [];
    cfg.channel = channel_set;
    
    temp = ft_selectdata(cfg, data);
    
    cfg = [];
    cfg.method = 'summary';
    cfg.keeptrial = 'nan';

    temp = ft_rejectvisual(cfg, temp);
    n_trials = length(temp.trial);
    
    % write removed trial indices file
    removed_trial_indices = [];
    for trial_index = 1:n_trials
        if isnan(temp.trial{trial_index}(1))
            removed_trial_indices = [removed_trial_indices ...
                                                trial_index]; %#ok<AGROW>
        end
    end
    % write a file for each channel type
    filename = [rejected_trials_file '_' channel_set output_type];
    fullpath = fullfile(subject_path, filename);
    
    
    if ~exist(fullpath, 'file') || overwrite
        dlmwrite(fullpath, removed_trial_indices, 'delimiter', ',');
        disp(['Wrote: ' fullpath]);
    else
        disp([fullpath ' already exists, not overwriting'])
    end
    
end

