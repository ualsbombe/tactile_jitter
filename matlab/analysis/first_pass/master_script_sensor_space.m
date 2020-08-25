%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

% find user
[~, user] = system('uname -n');
user = strtrim(user); % remove trailing and leading whitespaces
if strcmp(user, 'lau') % local
    addpath /home/lau/matlab/fieldtrip/
    project_path = fullfile('/', 'home', 'lau', 'projects', ...
                            'cerebellar_clock');
elseif strcmp(user, 'hyades02') % server
    addpath /users/lau/matlab/fieldtrip/
    project_path = fullfile('/', 'projects', ...
                            'MINDLAB2019_MEG-CerebellarClock');
end
ft_defaults

raw_path = fullfile(project_path, 'raw');
save_path = fullfile(project_path, 'scratch', 'tactile_jitter', 'MEG');
figures_path = fullfile(project_path, 'scratch', 'tactile_jitter', ...
                        'figures');
script_path = fullfile(project_path, 'scripts', 'tactile_jitter', ...
                       'matlab', 'analysis');
addpath(fullfile(script_path, 'functions'))                   

%% SUBJECTS

subjects = {
    '0001' '0002' '0004' '0005' '0006' '0007' ...
    '0008' '0009' '0010' '0011' '0012' '0013' ...
    '0014' '0015' '0016' '0017' '0018' '0019' ...
    '0020' '0021' '0022' '0023' '0024' '0025' ...
    '0026' '0027' '0028' '0029' '0030' '0031'
    };
dates = {
    '20200124' '20200121' '20200121' '20200128' '20200120' '20200120' ...
    '20200120' '20200121' '20200122' '20200122' '20200122' '20200124' ...
    '20200127' '20200128' '20200128' '20200129' '20200129' '20200129' ...
    '20200131' '20200131' '20200131' '20200203' '20200203' '20200203' ...
    '20200204' '20200204' '20200204' '20200205' '20200205' '20200205'
         };
     
split_recordings = {
                    0 0 0 0 1 0 ...
                    0 0 1 0 0 0 ...
                    0 0 0 0 0 0 ...
                    0 0 0 0 0 0 ...
                    0 0 0 0 0 0
                    };
                
n_subjects = length(subjects);                

%% SELECT REPOSITORY

select_repository_version

%% CLUSTER (quick and dirty)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 9 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).split_recording = ...
                                    split_recordings{subject_index};
    cluster_input(cluster_index).input_file = 'tactile_jitter_raw';
    cluster_input(cluster_index).save_path = save_path;
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@quick_and_dirty, cluster_input);
result = jobresult(jobid);

%% CLUSTER (plot quick and dirty)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).input_file = 'quick_and_dirty';
    cluster_input(cluster_index).figures_path = figures_path;
end
n_clusters = length(cluster_input);   

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@plot_quick_and_dirty, cluster_input);
result = jobresult(jobid);

%% CLUSTER (epoch and downsample)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1%:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).split_recording = ...
                                    split_recordings{subject_index};
    cluster_input(cluster_index).input_file = 'tactile_jitter_raw';
    cluster_input(cluster_index).downsampling_Hz = 250;
    cluster_input(cluster_index).save_path = save_path;
    cluster_input(cluster_index).output_type = '_epo';
    cluster_input(cluster_index).overwrite = true;
end
n_clusters = length(cluster_input);  

% CALL cluster

clusterconfig('wait', true, 'slot', 4, 'queue', 'all.q')
jobid = job2cluster(@epoch_and_downsample, cluster_input);
result = jobresult(jobid);


%% CLUSTER (epoch and downsample for timelocked)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).split_recording = ...
                                    split_recordings{subject_index};
    cluster_input(cluster_index).input_file = 'tactile_jitter_raw';
    cluster_input(cluster_index).downsampling_Hz = 250;
    cluster_input(cluster_index).save_path = save_path;
    cluster_input(cluster_index).output_type = '_epo';
    cluster_input(cluster_index).overwrite = true;
end
n_clusters = length(cluster_input);  

% CALL cluster

clusterconfig('wait', true, 'slot', 4, 'queue', 'all.q')
jobid = job2cluster(@epoch_and_downsample_timelocked, cluster_input);
result = jobresult(jobid);

%% CLEAN DATA (cannot be done on cluster, because it requires interaction)

cluster_input = struct();
cluster_index = 0;
for subject_index = 1:n_subjects
    
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject = subjects{subject_index};
    cluster_input(cluster_index).date = dates{subject_index};
    cluster_input(cluster_index).input_file = ...
                    'tactile_jitter_demeaned_downsampled_250_Hz_epo.mat';
    cluster_input(cluster_index).channel_sets = {'MEGMAG' 'MEGGRAD' ...
                                                 'EOG' 'EMG'};
    cluster_input(cluster_index).rejected_trials_file = ...
                                                    'rejected_trials';
    cluster_input(cluster_index).output_type = '.csv';
    cluster_input(cluster_index).overwrite = false;
    
end
n_clusters = length(cluster_input);

%% TIME FREQUENCY REPRESENTATION

cluster_input = struct();
cluster_index = 0;
for subject_index = 1:n_subjects
    
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject = subjects{subject_index};
    cluster_input(cluster_index).date = dates{subject_index};
    cluster_input(cluster_index).input_file = ...
                    'tactile_jitter_demeaned_downsampled_250_Hz_epo.mat';
    cluster_input(cluster_index).channel_sets = {'MEGGRAD' 'MEGMAG'};
    cluster_input(cluster_index).clean_based_on = {'MEGGRAD' 'MEGMAG'};
    cluster_input(cluster_index).events = {1 3 5 ...
                                           23 25 27 ...
                                           [71 87] [73 89] [75 91] ...
                                           [81 97] [83 99] [85 101] ...
                                           18 28 38 ...
                                           256};
    cluster_input(cluster_index).overwrite = true;
    
end
n_clusters = length(cluster_input);

clusterconfig('wait', false, 'slot', 2, 'queue', 'all.q')
jobid = job2cluster(@tfr, cluster_input);
result = jobresult(jobid);

%% TIMELOCKED (based on tfr preprocessed, for "untimelocking")

cluster_input = struct();
cluster_index = 0;
for subject_index = 1:n_subjects
    
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject = subjects{subject_index};
    cluster_input(cluster_index).date = dates{subject_index};
    cluster_input(cluster_index).input_file = ...
                    'tactile_jitter_demeaned_downsampled_250_Hz_epo.mat';
    cluster_input(cluster_index).channel_sets = {'MEGGRAD'};
    cluster_input(cluster_index).clean_based_on = {'MEGGRAD'};
    cluster_input(cluster_index).events = {
                                           1 3 5 ...
                                           23 25 27 ...
                                           [71 87] [73 89] [75 91] ...
                                           [81 97] [83 99] [85 101] ...
                                           18 28 38 ...
                                           256
                                           };
    cluster_input(cluster_index).lowpass_filter = [];
    cluster_input(cluster_index).highpass_filter = [];                                   
    cluster_input(cluster_index).output_file = ...
                                        'timelocked_for_untimelocked.mat';
    cluster_input(cluster_index).overwrite = true;
    
end
n_clusters = length(cluster_input);

clusterconfig('wait', true, 'slot', 4, 'queue', 'all.q')
jobid = job2cluster(@timelock, cluster_input);
result = jobresult(jobid);


%% TIMELOCKED (for final analysis)

cluster_input = struct();
cluster_index = 0;
for subject_index = 1:n_subjects
    
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject = subjects{subject_index};
    cluster_input(cluster_index).date = dates{subject_index};
    cluster_input(cluster_index).input_file = ...
    'tactile_jitter_demeaned_for_timelocked_downsampled_250_Hz_epo.mat';
    cluster_input(cluster_index).channel_sets = {'MEGGRAD' 'MEGMAG'};
    cluster_input(cluster_index).clean_based_on = {'MEGGRAD' 'MEGMAG'};
    cluster_input(cluster_index).events = {
                                           1 3 5 ...
                                           23 25 27 ...
                                           [71 87] [73 89] [75 91] ...
                                           [81 97] [83 99] [85 101] ...
                                           18 28 38 ...
                                           256
                                           };
    cluster_input(cluster_index).lowpass_filter = 40; %Hz
    cluster_input(cluster_index).highpass_filter = 1; %Hz
    cluster_input(cluster_index).output_file = ...
                                        'timelocked.mat';
    cluster_input(cluster_index).overwrite = true;
    
end
n_clusters = length(cluster_input);

clusterconfig('wait', true, 'slot', 4, 'queue', 'all.q')
jobid = job2cluster(@timelock, cluster_input);
result = jobresult(jobid);

%% CALL FUNCTIONS (locally)

% quick_and_dirty: quick and dirty processing to see if there's any
    % evidence of a cerebellar activity for omissions

for cluster_index = 1:n_clusters
    
%     tfr_contrast = quick_and_dirty(cluster_input(cluster_index));    
%     epoch_and_downsample(cluster_input(cluster_index));
%     plot_quick_and_dirty(cluster_input(cluster_index));
%     clean_data(cluster_input(cluster_index));
    tfr(cluster_input(cluster_index));
%     timelock(cluster_input(cluster_index));
                           
end
                           
%% PLOT

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.xlim = [0.000 1.500];
cfg.zlim = [0.700 1.300];
% cfg.ylim = [9 40];
% cfg.channel = {'MEG' '-MEG2542+2543'};

ft_multiplotTFR(cfg, tfr_contrast);

%% QUICK AND DIRTY GRAND AVERAGE

ga_cell = {};

for subject_index = [1:8 10:(n_subjects-1)]
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    
    fullpath = fullfile(save_path, subject, full_date, ...
                        'quick_and_dirty.mat');
    
    try
        disp(['Loading subject: ' subject])
        load(fullpath) %% loads tfr_contrast
        ga_cell{length(ga_cell) + 1} = tfr_contrast; %#ok<*SAGROW>
    catch
        disp(['Subject ' subject ' could not be loaded'])
    end
    
end

cfg = [];

ga = ft_freqgrandaverage(cfg, ga_cell{:});

%% PLOT GA

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [0.8 1.2];
cfg.channel = {'MEG' '-MEG0522+0523' '-MEG0712+0713' '-MEG0922+0923' ...
               '-MEG1012+1013'};
cfg.showlabels = 'yes';           

ft_multiplotTFR(cfg, ga);


%% GA TIMELOCKEDS

ga_cell = {};

for subject_index = 1:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, ...
                        'timelocked.mat');
    disp(['Loading subject: ' subject])          
    tl = load(fullpath);
    ga_cell{length(ga_cell) + 1} = tl.event_18;
end

cfg = [];
ga = ft_timelockgrandaverage(cfg, ga_cell{:});