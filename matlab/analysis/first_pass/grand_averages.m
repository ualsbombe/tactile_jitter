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
ga_path = fullfile(save_path, 'grand_averages');
figures_path = fullfile(project_path, 'scratch', 'tactile_jitter', ...
                        'figures');
script_path = fullfile(project_path, 'scripts', 'tactile_jitter', ...
                       'matlab', 'analysis');
addpath(fullfile(script_path, 'functions'))                   

%% SUBJECTS

subjects = {
    '0001' '0002' '0004' '0005' '0006' '0007' ...
    '0008' '0009' '0010' '0011' '0012' '0013' ...
    ...'0014' ...
    '0015' '0016' '0017' '0018' '0019' ...
    '0020' '0021' '0022' '0023' '0024' '0025' ...
    '0026' '0027' '0028' '0029' '0030' '0031'
    };
dates = {
    '20200124' '20200121' '20200121' '20200128' '20200120' '20200120' ...
    '20200120' '20200121' '20200122' '20200122' '20200122' '20200124' ...
    ...'20200127' ...
    '20200128' '20200128' '20200129' '20200129' '20200129' ...
    '20200131' '20200131' '20200131' '20200203' '20200203' '20200203' ...
    '20200204' '20200204' '20200204' '20200205' '20200205' '20200205'
         };
     
split_recordings = {
                    0 0 0 0 1 0 ...
                    0 0 1 0 0 0 ...
                    ...0 ...
                    0 0 0 0 0 ...
                    0 0 0 0 0 0 ...
                    0 0 0 0 0 0
                    };
                
n_subjects = length(subjects);

%% GA TIMELOCKEDS

n_events = 16;
ga_cell = cell(n_events, n_subjects);
filename = 'highpass_1_Hz_lowpass_40_Hz_timelocked.mat';

for subject_index = 1:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, filename);
    disp(['Loading subject: ' subject])          
    tfrs = load(fullpath);
    names = fieldnames(tfrs);
    n_names = length(names);
    for name_index = 1:n_names
        name = names{name_index};
        ga_cell{name_index, subject_index} = tfrs.(name);
    end
    
end

gas = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    this_ga = ft_timelockgrandaverage(cfg, ...
                                                  ga_cell{event_index, :});
    this_ga.cfg.previous = []; %% drastically reduces size
    gas{event_index} = this_ga;
end

fullpath = fullfile(ga_path, filename);
save(fullpath, 'gas')

%% GA CONTRASTS

ga_cell = cell(1, n_subjects);
counter = 1;

for subject_index = 1:n_subjects %[1:11 13:24 26:n_subjects]%    26:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, 'tfrs');
    disp(['Loading subject: ' subject])          
    tfrs = load(fullpath);
    
    cfg = [];
    cfg.operation = '(x1 - x2) ./ (x1 + x2)';
    cfg.parameter = 'powspctrm';
    
    cmb_1 = ft_combineplanar([], tfrs.event_18);
    cmb_2 = ft_combineplanar([], tfrs.event_256);
    
    ga_cell{counter} = ft_math(cfg, cmb_1, cmb_2);
    counter = counter + 1;
    
end

cfg = [];

ga = ft_freqgrandaverage(cfg, ga_cell{:});
ga.cfg.previous = [];

fullpath = fullfile(ga_path, 'omission_vs_non-stimulation');
save(fullpath, 'ga')


%% LOAD SUBJECTS

ga_cell = cell(2, n_subjects);
counter = 1;

for subject_index = 1:n_subjects %[1:11 13:24 26:n_subjects]%    26:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, 'tfrs');
    disp(['Loading subject: ' subject])          
    tfrs = load(fullpath);
       
    cmb_1 = ft_combineplanar([], tfrs.event_18);
    cmb_2 = ft_combineplanar([], tfrs.event_256);
    
    ga_cell{1, counter} = cmb_1;
    ga_cell{2, counter} = cmb_2;
    counter = counter + 1;
    
end

%% LOAD SUBJECTS 2

ga_cell_2 = cell(2, n_subjects);
counter = 1;

for subject_index = 1:n_subjects %[1:11 13:24 26:n_subjects]%    26:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, 'tfrs');
    disp(['Loading subject: ' subject])          
    tfrs = load(fullpath);
       
    cmb_1 = ft_combineplanar([], tfrs.event_28);
    cmb_2 = ft_combineplanar([], tfrs.event_256);
    
    ga_cell_2{1, counter} = cmb_1;
    ga_cell_2{2, counter} = cmb_2;
    counter = counter + 1;
    
end


%% LOAD SUBJECTS 3

ga_cell_3 = cell(2, n_subjects);
counter = 1;

for subject_index = 1:n_subjects %[1:11 13:24 26:n_subjects]%    26:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    fullpath = fullfile(save_path, subject, full_date, 'tfrs');
    disp(['Loading subject: ' subject])          
    tfrs = load(fullpath);
       
    cmb_1 = ft_combineplanar([], tfrs.event_38);
    cmb_2 = ft_combineplanar([], tfrs.event_256);
    
    ga_cell_3{1, counter} = cmb_1;
    ga_cell_3{2, counter} = cmb_2;
    counter = counter + 1;
    
end
%% temporary crao

ga_cell = ga_cell_2;


%% STATS

% one_cell = cell(1, n_subjects);
% 
% for subject_index = 1:n_subjects
%     
%     temp = ga_cell{subject_index};
%     temp.powspctrm = ones(size(temp.powspctrm));
%     one_cell{subject_index} = temp;
%     
% end

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_depsamplesT';
cfr.correctm = 'fdr';
cfg.alpha = 0.05;
cfg.design(1, :) = [ones(1, n_subjects) 2*ones(1, n_subjects)];
cfg.design(2, :) = [1:n_subjects 1:n_subjects];
cfg.ivar = 1;
cfg.uvar = 2;

stat = ft_freqstatistics(cfg, ga_cell{1, :}, ga_cell{2, :});

%% plot

figure

cfg = [];
% cfg.parameter = 'stat';
% cfg.maskparameter = 'mask';
% cfg.maskstyle = 'outline';
cfg.layout = 'neuromag306cmb.lay';
% cfg.baseline = 'yes';
% cfg.baselinetype = 'relative';
cfg.zlim = [-0.03 0.03];
% cfg.channel = {'MEG' '-MEG1012+1013'};

% ft_multiplotTFR(cfg, stat);
ft_multiplotTFR(cfg, ga);


% ft_multiplotTFR(cfg, ft_freqgrandaverage([], ga_cell{1, :}));
% figure
% ft_multiplotTFR(cfg, ft_freqgrandaverage([], ga_cell{2, :}));

%% cluster

cfg = [];
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
% cfg.channel =  {'MEG1832+1833 MEG2012+2013 MEG2022+2023 MEG2242+MEG2243'};

neighbours = ft_prepare_neighbours(cfg);

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.alpha = 0.05;
cfg.design(1, :) = [ones(1, n_subjects) 2*ones(1, n_subjects)];
cfg.design(2, :) = [1:n_subjects 1:n_subjects];
cfg.ivar = 1;
cfg.uvar = 2;
cfg.neighbours = neighbours;
% cfg.frequency = [6 15];
% cfg.latency = [0.000 1.000];
% cfg.channel = {'MEG1832+1833 MEG2012+2013 MEG2022+2023 MEG2242+MEG2243'};

cfg.alpha = 0.25;
cfg.clusteralpha = 0.05;
cfg.numrandomization = 50;

stat = ft_freqstatistics(cfg, ga_cell{1, :}, ga_cell{2, :});