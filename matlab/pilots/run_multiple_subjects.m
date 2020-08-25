%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

[~, user] = system('uname -n');
if strcmp(user(1:(end-1)), 'lau-LIFEBOOK-E744')
    home_dir = '/home/lau/';
elseif strcmp(user, 'hyades02')
    home_dir = '/users/lau/';
else
    disp('"User" has not been specified yet')
end

project_name = 'MINDLAB2019_MEG-CerebellarClock';
fieldtrip_path = '/home/lau/matlab/fieldtrip';     

addpath(fieldtrip_path)
ft_defaults


subjects = {'0001' '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' '20190930_000000' '20190930_000000' ...
         '20191008_000000' '20191008_000000'};

n_subjects = length(subjects);

%% RUN SUBJECT LOOP

% toilim = [0.600 0.900];
% frequency = 80;
% tapsmofrq = 20;

% cerebellum and insula preparation?
% toilim = [-0.700 -0.300];
% frequency = 10;
% tapsmofrq = 4;

% cerebellum error?
% toilim = [0.600 0.800]; 
% frequency = 15;
% tapsmofrq = 4;

for subject_index = 1:n_subjects
     
    subject = subjects{subject_index};
    date = dates{subject_index};

    data_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                          project_name, 'scratch', 'tactile_jitter', ...
                          'MEG', subject, date);
                      
%     source_analysis_function(data_path, toilim, frequency, tapsmofrq);
    gamma_oscillations_function(data_path);
    
end
