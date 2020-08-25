function [tfr_contrast] = quick_and_dirty(input)
                                      
raw_path = input.raw_path;
subject = input.subject;
date = input.date;
split_recording = input.split_recording;
input_file = input.input_file;
save_path = input.save_path;

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
    cfg.trialdef.prestim = 1.500;
    cfg.trialdef.poststim = 1.500;
    cfg.trialdef.eventvalue = {18 48:2:76};
    cfg.trialdef.eventtype = 'STI101';

    cfg = ft_definetrial(cfg);  
    cfg.demean = 'yes';
    cfg.baselinewindow = [-Inf Inf];

    data = ft_preprocessing(cfg);

    cfg = [];
    cfg.resamplefs = 250; % Hz

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


%% TFRs

events = {18 256};
n_events = length(events);

tfrs = cell(1, n_events);

for event_index = 1:n_events

    event = events{event_index};
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
    cfg.method = 'wavelet';
    cfg.width = 7;
    cfg.toi = -1.500:0.020:1.500;
    cfg.foilim = [1 40];
    cfg.trials = find(data.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    tfrs{event_index} = ft_freqanalysis(cfg, data);
    tfrs{event_index}.grad = grad;
    
end

%% COMBINE TFRS

tfr_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tfr_cmbs{event_index} = ft_combineplanar(cfg, tfrs{event_index});
    
end

%% OMISSION RATIO

cfg = [];
cfg.operation = 'x1 / x2';
cfg.parameter = 'powspctrm';

tfr_contrast = ft_math(cfg, tfr_cmbs{1}, tfr_cmbs{2});

%% SAVE PATH

full_save_path = fullfile(save_path, subject, full_date, ...
    'quick_and_dirty.mat');

save(full_save_path, 'tfr_cmbs', 'tfr_contrast');
