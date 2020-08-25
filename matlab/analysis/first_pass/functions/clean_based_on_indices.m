function [cleaned_data] = clean_based_on_indices(subject_path, ...
                                                 channel_sets, data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

trials_to_remove = []; %#ok<*AGROW>
n_clean_sets = length(channel_sets);

for clean_index = 1:n_clean_sets
    
    cleaning_set =  channel_sets{clean_index};
    fullpath = fullfile(subject_path, ['rejected_trials_' cleaning_set ...
                                        '.csv']);
    try
        indices = dlmread(fullpath, ',');
    catch ME
        if strcmp(ME.identifier, 'MATLAB:textscan:EmptyFormatString')
            indices = [];
        else
            error(['File with indices to be cleaned (probably)' ...
                    'does not exist']);
        end
    end
    trials_to_remove = [trials_to_remove indices];

end

% find unique trials
trials_to_remove = unique(trials_to_remove);

% clean
n_trials = length(data.trial);
trials = 1:n_trials;
trials(trials_to_remove) = []; % remove trials

cfg = [];
cfg.trials = trials;

cleaned_data = ft_selectdata(cfg, data);

end

