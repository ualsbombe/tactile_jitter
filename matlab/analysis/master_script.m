%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

% find user
[~, user] = system('uname -n');
user = strtrim(user); % remove trailing and leading whitespaces
if strcmp(user, 'lau') % local
    fieldtrip_path = '/home/lau/matlab/fieldtrip/';
    project_path = fullfile('/', 'home', 'lau', 'projects', ...
                            'cerebellar_clock');
elseif strcmp(user, 'hyades02') % server
    fieldtrip_path = '/users/lau/matlab/fieldtrip/';
    project_path = fullfile('/', 'projects', ...
                            'MINDLAB2019_MEG-CerebellarClock');
end
addpath(fieldtrip_path);
ft_defaults

raw_path = fullfile(project_path, 'raw');
scratch_path = fullfile(project_path, 'scratch', 'tactile_jitter', 'MEG');
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

%% CREATE INPUT STRUCTURE; PARCELLATE NIFTI

input = struct();
cluster_index = 0; % counter

for subject_index = 1:n_subjects
    
    cluster_index = cluster_index + 1;
    subject = subjects{subject_index};
    date = [dates{subject_index} '_000000'];
    
    input(cluster_index).input_path = fullfile(scratch_path, subject, ...
                                               date, 'hilbert');
                                           
    input(cluster_index).input_file = 'tactile_jitter_morph-vl.nii';
    input(cluster_index).h_freqs = {'14'};
    input(cluster_index).l_freqs = {'30'};
%     input(cluster_index).h_freqs = {'4' '8' '25' '60'};
%     input(cluster_index).l_freqs = {'7' '12' '35' '85'};
%     input(cluster_index).h_freqs = {'4' '8' '14' '25' '60'};
%     input(cluster_index).l_freqs = {'7' '12' '30' '35' '85'};
    input(cluster_index).spacing = '7.5';
    input(cluster_index).reg = '0.0';
    input(cluster_index).weight_norm = 'unit-noise-gain';
    input(cluster_index).channel_type = 'grad';
    input(cluster_index).contrasts = {{'o0' 'o15'} {'o0' 'o5'} ...
                                     {'o5' 'o15'}};
    input(cluster_index).tmin = -0.400;
    input(cluster_index).tmax = 1.100;
    input(cluster_index).ISI_type = 'normal';

    % input from fieldtrip templates
    input(cluster_index).mri = fullfile(fieldtrip_path, 'template', ...
                                        'anatomy', 'single_subj_T1.nii');
    input(cluster_index).atlas = fullfile(fieldtrip_path, 'template', ...
                                          'atlas', 'aal', ...
                                          'ROI_MNI_V4.nii');

    % input to ft_sourceparcellate                   
    input(cluster_index).parcel_method = 'mean';                   

    % output name
    input(cluster_index).overwrite = false;
    input(cluster_index).output_path = input(cluster_index).input_path;
    input(cluster_index).output_file = 'tactile_jitter_parcel.mat';
    
end
n_clusters = length(input);

clusterconfig('wait', false, 'slot', 2, 'queue', 'long.q')
jobid = job2cluster(@parcellate_nifti, input(1:end));
result = jobresult(jobid);

%% CREATE INPUT STRUCTURE; PARCELLATE NIFTI STAT

input = struct();

% input from MNE-Python
input.type = 'stcs';
input.input_path = fullfile(scratch_path, 'grand_averages', ...
                            'statistics', 'hilbert', input.type, ...
                            'common_filter', 'contrasts');
input.input_file = 'tactile_jitter_morph-vl.nii';
input.h_freq = 14;
input.l_freq = 30;
input.spacing = '7.5';
input.reg = '0.0';
input.weight_norm = 'unit-noise-gain';
input.channel_type = 'mag_and_grad';
input.contrasts = {'o0' 'o15'};
input.p_threshold = 0.05;
input.tmin = -0.400;
input.tmax = 1.100;
input.ISI_type = 'normal';
input.stat_function = 't-test';

% input from fieldtrip templates
input.mri = fullfile(fieldtrip_path, 'template', 'anatomy', ...
                    'single_subj_T1.nii');
input.atlas = fullfile(fieldtrip_path, 'template', 'atlas', 'aal', ...
                       'ROI_MNI_V4.nii');

% input to ft_sourceparcellate                   
input.parcel_method = 'max';                   

% output name
input.overwrite = false;
input.output_path = input.input_path;
input.output_file = 'tactile_jitter_parcel.mat';

%% RUN LOCALLY

parcellate_nifti(input(1));
parcellate_nifti_stat(input);

%% CLUSTER CALL
clusterconfig('wait', false, 'slot', 2, 'queue', 'all.q')
jobid = job2cluster(@parcellate_nifti_stat, input);
result = jobresult(jobid);