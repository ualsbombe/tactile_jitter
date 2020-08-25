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
save_path = fullfile(project_path, 'scratch', 'tactile_jitter', 'MEG');
figures_path = fullfile(project_path, 'scratch', 'tactile_jitter', ...
                        'figures');
script_path = fullfile(project_path, 'scripts', 'tactile_jitter', ...
                       'matlab', 'analysis');
addpath(fullfile(script_path, 'functions'))                   

%% SUBJECTS

subjects = {
    '0001' %'0002' '0004' ...
    %'0005' %'0006' '0007' ...
%     '0008' '0009' '0010' ...
%     '0011' '0012' '0013' ...
%     '0014' '0015' '0016' ...
%     '0017' '0018' '0019' ...
%     '0020' '0021' '0022' ...
%     '0023' '0024' '0025' ...
%     '0026' '0027' '0028' ...
%     '0029' '0030' '0031'
    };
mr_dates = {
    '20180515_132117' '20191015_104445' '20191015_112257' ...
    '20191015_121553' 
        };
    
meg_dates = {
    '20200124' '20200121' '20200121' '20200128' %'20200120' '20200120' ...
%     '20200120' '20200121' '20200122' '20200122' '20200122' '20200124' ...
%     '20200127' '20200128' '20200128' '20200129' '20200129' '20200129' ...
%     '20200131' '20200131' '20200131' '20200203' '20200203' '20200203' ...
%     '20200204' '20200204' '20200204' '20200205' '20200205' '20200205'
         };    
     
                
n_subjects = length(subjects);

%% SELECT REPOSITORY

select_repository_version

%% CLUSTER (source model)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).fieldtrip_path = fieldtrip_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).mr_date     = mr_dates{subject_index};
    cluster_input(cluster_index).meg_date     = meg_dates{subject_index};
    cluster_input(cluster_index).input_file = 'mri_resliced';
    cluster_input(cluster_index).template = 'standard_sourcemodel3d10mm';
    cluster_input(cluster_index).output_file = 'warped_sourcemodel.mat';
    cluster_input(cluster_index).save_path = save_path;
    cluster_input(cluster_index).overwrite = false;
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@create_sourcemodel, cluster_input);
result = jobresult(jobid);

%% CHECK SOURCE MODEL

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;    
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).mr_date     = mr_dates{subject_index};
    cluster_input(cluster_index).meg_date     = meg_dates{subject_index};
    cluster_input(cluster_index).input_files = { ...
                                        'mesh_scalp' 'warped_sourcemodel'};
    check_sourcemodel(cluster_input(cluster_index));
end
n_clusters = length(cluster_input);

%% CLUSTER (lead field)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;    
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).mr_date     = mr_dates{subject_index};
    cluster_input(cluster_index).meg_date     = meg_dates{subject_index};
    cluster_input(cluster_index).input_files = { ...
                            'single_shell_headmodel' 'warped_sourcemodel'};
    cluster_input(cluster_index).normalize = 'yes';
    cluster_input(cluster_index).output_file = ...
        'warped_leadfield_normalized.mat';
    cluster_input(cluster_index).overwrite = false;
end
n_clusters = length(cluster_input);          %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', true, 'slot', 1, 'queue', 'short.q')
jobid = job2cluster(@create_leadfield, cluster_input);
result = jobresult(jobid);

%% CLUSTER (beamformer contrast)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).fieldtrip_path = fieldtrip_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).mr_date     = mr_dates{subject_index};
    cluster_input(cluster_index).meg_date     = meg_dates{subject_index};
    cluster_input(cluster_index).input_files = {
        'tactile_jitter_demeaned_downsampled_250_Hz_epo.mat' ...
        'single_shell_headmodel.mat' ...
        'warped_leadfield.mat'};
    cluster_input(cluster_index).template = 'standard_sourcemodel3d10mm';
    cluster_input(cluster_index).clean_based_on = {'MEGGRAD'};
    cluster_input(cluster_index).frequency = 15;
    cluster_input(cluster_index).smoothing = 7;
    cluster_input(cluster_index).toilim = [-0.300 -0.150];
    cluster_input(cluster_index).output_file = 'beamformers';
    cluster_input(cluster_index).overwrite = false;
end
n_clusters = length(cluster_input);          %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@common_filter_beamformer, cluster_input);
result = jobresult(jobid);

%% CLUSTER (beamformer)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1 %19:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).fieldtrip_path = fieldtrip_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).mr_date     = mr_dates{subject_index};
    cluster_input(cluster_index).meg_date     = meg_dates{subject_index};
    cluster_input(cluster_index).input_files = {
        'tactile_jitter_demeaned_downsampled_250_Hz_epo.mat' ...
        'single_shell_headmodel.mat' ...
        'warped_leadfield_normalized.mat'};
    cluster_input(cluster_index).template = 'standard_sourcemodel3d10mm';
    cluster_input(cluster_index).clean_based_on = {'MEGGRAD'};
    cluster_input(cluster_index).frequency = 16;
    cluster_input(cluster_index).smoothing = 7;
    cluster_input(cluster_index).toilim = [0.500 1.000];
    cluster_input(cluster_index).output_file = 'beamformers_normalized';
    cluster_input(cluster_index).overwrite = false;
end
n_clusters = length(cluster_input);          %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', false, 'slot', 1, 'queue', 'short.q')
jobid = job2cluster(@beamformer, cluster_input);
result = jobresult(jobid);

%% CALL FUNCTIONS (locally)

% quick_and_dirty: quick and dirty processing to see if there's any
    % evidence of a cerebellar activity for omissions

for cluster_index = 1%:n_clusters
    
%     create_sourcemodel(cluster_input(cluster_index));
    create_leadfield(cluster_input(cluster_index));
    beamformer(cluster_input(cluster_index));
%     common_filter_beamformer(cluster_input(cluster_index));
    
end