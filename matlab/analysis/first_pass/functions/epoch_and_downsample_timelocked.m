function [dummy] = epoch_and_downsample_timelocked(input)


dummy = 'just_a_dummy_variable';

raw_path = input.raw_path;
subject = input.subject;
date = input.date;
split_recording = input.split_recording;
input_file = input.input_file;
downsampling_Hz = input.downsampling_Hz;
output_type = input.output_type;
save_path = input.save_path;
overwrite = input.overwrite;

%% FIND ALL SPLITS (SEPARATE RECORDINGS)                              
                              
full_date = [date '_000000'];
subject_path = fullfile(raw_path, subject, full_date, 'MEG');
n_splits = split_recording + 1;
paths = cell(1, n_splits);
for split_index = 1:n_splits
    if n_splits > 1
        split_indicator = ['_' num2str(split_index)];
    else
        split_indicator = '';
    end
    file_path = [sprintf('%03d', split_index) '.tactile_jitter_raw' ...
                 split_indicator];
    paths{split_index} = fullfile(subject_path, file_path, 'files');
end

%% FIND ALL FILES: SPLITS DUE TO 2GB LIMIT

n_paths = length(paths);
full_paths = {};
for path_index = 1:n_paths
    path = paths{path_index};
    n_files = length(dir(path)) - 2; % subtracting . and ..
    filenames = cell(1, n_files);
    for file_index = 1:n_files
        
        if n_paths > 1
            split_indicator = ['_' num2str(path_index)];
        else
            split_indicator =  '';
        end
        
        if file_index == 1
            filename = [input_file split_indicator '.fif'];
        else
            filename = [input_file split_indicator ...
                        '-'  num2str(file_index - 1) '.fif'];
        end
        
        filenames{file_index} = filename;
        full_path = fullfile(path, filename);
        full_paths{length(full_paths) + 1} = full_path; %#ok<AGROW>
    end
end

%% DEFINE AND PREPROCESS

n_files = length(full_paths); % new n_files
split_files = cell(1, n_files);

for file_index = 1:n_files

    filename = full_paths{file_index};
    
    cfg = [];
    cfg.dataset = filename;
    cfg.trialdef.prestim = 0.200;
    cfg.trialdef.poststim = 0.750;
    cfg.trialdef.eventvalue = {1 3 5  ...
                               23 25 27 ...
                               71 73 75 ...
                               81 83 85 ...
                               87 89 91 ...
                               97 99 101 ...
                               18 28 38 ...
                               48:2:76};
    cfg.trialdef.eventtype = 'STI101';

    cfg = ft_definetrial(cfg);  
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf 0];

    data = ft_preprocessing(cfg);
    
    cfg = [];
    cfg.resamplefs = downsampling_Hz;

    data = ft_resampledata(cfg, data);
    
    split_files{file_index} = data;
    
    % get grad from first recording
    if file_index == 1
        grad = data.grad;
    end
    
end

cfg = [];

data = ft_appenddata(cfg, split_files{:});

non_stim_triggers = 48:2:76;


for trigger = non_stim_triggers
    
    for trial_index = 1:length(data.trialinfo)
        
        trial = data.trialinfo(trial_index);
        if trigger == trial
            data.trialinfo(trial_index) = 256;
        end
    end
end

data.grad = grad; %% might not be the best for split recordings

%% SAVE
warning('we got to here')
filename = ['tactile_jitter_demeaned_for_timelocked_downsampled_' ...
            num2str(downsampling_Hz) '_Hz' output_type '.mat'];
fullpath = fullfile(save_path, subject, full_date, filename);

save_struct(data, fullpath, overwrite);

end
