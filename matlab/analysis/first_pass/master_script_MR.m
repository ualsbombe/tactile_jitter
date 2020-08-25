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
    '0001' '0002' '0004' ...
    '0005' %'0006' '0007' ...
%     '0008' '0009' '0010' ...
%     '0011' '0012' '0013' ...
%     '0014' '0015' '0016' ...
%     '0017' '0018' '0019' ...
%     '0020' '0021' '0022' ...
%     '0023' '0024' '0025' ...
%     '0026' '0027' '0028' ...
%     '0029' '0030' '0031'
    };
dates = {
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

%% CLUSTER (read mri and convert to mat)

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 2:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).recording = 't1_mprage_3D_sag_fatsat';
    cluster_input(cluster_index).save_path = save_path;
    cluster_input(cluster_index).output_file = 'mri.mat';
    cluster_input(cluster_index).overwrite = false;
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>

% CALL cluster

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@convert_mri_to_mat, cluster_input);
result = jobresult(jobid);

%% ALIGN WITH FIDUCIALS

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 4:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).input_file = 'mri.mat';
    cluster_input(cluster_index).output_file = ...
                                'mri_realigned_fiducial.mat';
    cluster_input(cluster_index).overwrite = true;
    realign_mri_to_fiducials(cluster_input(cluster_index));
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>

%% ALIGN WITH HEAD POINTS

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 2:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).meg_date = meg_dates{subject_index};
    cluster_input(cluster_index).input_file = 'mri_realigned_fiducial.mat';
    cluster_input(cluster_index).output_file = ...
                                           'mri_realigned_head_points.mat';
    cluster_input(cluster_index).overwrite = true;
    realign_mri_to_head_points(cluster_input(cluster_index));
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>

%% CHECK FINAL ALIGNMENT

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 2:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).raw_path = raw_path;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).meg_date = meg_dates{subject_index};
    cluster_input(cluster_index).input_file = ...
                                           'mri_realigned_head_points.mat';
    check_realignment(cluster_input(cluster_index));
end
n_clusters = length(cluster_input);                 %#ok<*NASGU>


%% RESLICE, SEGMENT, CREATE MESHES AND HEADMODEL

cluster_input = struct();
cluster_index = 0; %% just a counter
for subject_index = 1%:n_subjects
    cluster_index = cluster_index + 1;
    cluster_input(cluster_index).load_path = save_path;
    cluster_input(cluster_index).subject  = subjects{subject_index};
    cluster_input(cluster_index).date     = dates{subject_index};
    cluster_input(cluster_index).input_file = ...
                                           'mri_realigned_head_points.mat';
    cluster_input(cluster_index).output_files = ...
                                           {'mri_resliced.mat' ...
                                            'mri_segmented.mat' ...
                                            'mesh_brain.mat' ...
                                            'mesh_scalp.mat' ...
                                            'single_shell_headmodel.mat'
                                            };
    cluster_input(cluster_index).overwrite = true;
    
end
n_clusters = length(cluster_input);    

clusterconfig('wait', true, 'slot', 1)
jobid = job2cluster(@reslice_segment_and_create_headmodel, cluster_input);
result = jobresult(jobid);


%% CALL FUNCTIONS (locally)


for cluster_index = 1:n_clusters
    
%     convert_mri_to_mat(cluster_input(cluster_index));
%     realign_mri_to_fiducials(cluster_input(cluster_index));
%     realign_mri_to_head_points(cluster_input(cluster_index));
%     check_realignment(cluster_input(cluster_index));
    reslice_segment_and_create_headmodel(cluster_input(cluster_index));
end