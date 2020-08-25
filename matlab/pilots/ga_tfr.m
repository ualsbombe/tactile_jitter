%% CLEAR AND ADD PATHS

clear variables
restoredefaultpath

subjects = {'0001' ...
            '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' ...
    '20190930_000000' '20190930_000000' ...
         '20191008_000000' '20191008_000000'};
n_subjects = length(subjects);
n_contrasts = 3;
n_conditions = 4;
project_name = 'MINDLAB2019_MEG-CerebellarClock';
fieldtrip_path = '/home/lau/matlab/fieldtrip';     

addpath(fieldtrip_path)
ft_defaults

%% COLLECT DATA

ga_cell = cell(n_subjects, n_contrasts);

for subject_index = 1:n_subjects

    subject = subjects{subject_index};
    date = dates{subject_index};

    data_path = fullfile('/', 'home', 'lau', 'mounts', 'hyades', 'projects', ...
                          project_name, 'scratch', ...
                          'tactile_jitter', 'MEG', subject, date);
                                          
    load(fullfile(data_path, 'tfr_cmb_contrasts'))
    for contrast_index = 1:n_contrasts
        temp = tfr_cmb_contrasts{contrast_index};
        ga_cell{subject_index, contrast_index} = temp;
    end
end

% stat_cell = cell(n_subjects, n_conditions);
% 
% for subject_index = 1:n_subjects
% 
%     subject = subjects{subject_index};
%     date = dates{subject_index};
% 
%     data_path = fullfile('/', 'home', 'lau', 'mounts', 'hyades', 'projects', ...
%                           project_name, 'scratch', ...
%                           'tactile_jitter', 'MEG', subject, date);
%                                            
%     load(fullfile(data_path, 'tfrs'))
%     for condition_index = 1:n_conditions
%         temp = ft_combineplanar([], tfrs{condition_index});
%         stat_cell{subject_index, condition_index} = temp;
%     end
% end


%% GRAND AVERAGES

GAs = cell(1, n_contrasts);
cfg = [];

for contrast_index = 1:n_contrasts
    GAs{contrast_index} = ft_freqgrandaverage(cfg, ...
                                            ga_cell{:, contrast_index});
                                        
end

%% DEFAULTS

set(0, 'defaultaxesfontweight', 'bold', 'defaultaxesfontsize', 14)

%% PLOT GAs

% close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.ylim = [2 40];
cfg.ylim = [40 100];
cfg.zlim = [-0.150 0.150];
cfg.zlim = [-0.050 0.050];
cfg.colorbar = 'yes';
cfg.colorbartext = 'Relative power';
cfg.fontsize = 14;
cfg.fontweight = 'bold';

figure
ft_multiplotTFR(cfg, GAs{1});
figure
ft_multiplotTFR(cfg, GAs{2});
figure
ft_multiplotTFR(cfg, GAs{3});

%% TOPO PLOT 3D

load('standard_singleshell.mat') % loads vol
headmodel = ft_convert_units(vol, 'm');
headmodel.bnd.pos(:, 3) = headmodel.bnd.pos(:, 3) + 0.040;

data = GAs{1};

cfg = [];
cfg.frequency = [10 10]; %[15 15];
cfg.latency = [-0.700 -0.300]; %[0.600 0.800]; % cerebellar
cfg.frequency = [15 15];
cfg.latency = [0.600 0.800];

data = ft_selectdata(cfg, data);

data.powspctrm = squeeze(mean(data.powspctrm, 3));
data.powspctrm = squeeze(mean(data.powspctrm, 2));

temp.grad = ft_convert_units(temp.grad, 'm'); % last subject's grad

figure; hold on
view(0, 0)
ft_plot_topo3d(temp.grad.chanposold(1:3:306, :), data.powspctrm(1:102), ...
                'facealpha', 0.9);
% ft_plot_sens(temp.grad);
ft_plot_headmodel(headmodel);
c = colorbar;
caxis([-0.10 0.10])


%% STATS MONTECARLO

cfg = [];
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';

neighbours = ft_prepare_neighbours(cfg);

events = [1 4];
temp = [stat_cell(:, 1); stat_cell(:, 4)];

cfg = [];
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.numrandomization = 100;
cfg.alpha = 0.05;
cfg.clusteralpha = 0.05;
cfg.statistic = 'depsamplesT';
cfg.design(1, :) = [1:n_subjects 1:n_subjects];
cfg.design(2, :) = [ones(1, n_subjects), 2*ones(1, n_subjects)];
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = neighbours;
cfg.frequency = [2 30];
cfg.latency = [0.400 1.000];

stat = ft_freqstatistics(cfg, temp{:});

%% STATS MONTECARLO GAMMA

cfg = [];
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';

neighbours = ft_prepare_neighbours(cfg);

events = [3 4];
temp = [stat_cell(:, 1); stat_cell(:, 4)];

cfg = [];
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.numrandomization = 100;
cfg.alpha = 0.05;
cfg.clusteralpha = 0.05;
cfg.statistic = 'depsamplesT';
cfg.design(1, :) = [1:n_subjects 1:n_subjects];
cfg.design(2, :) = [ones(1, n_subjects), 2*ones(1, n_subjects)];
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = neighbours;
cfg.frequency = [40 100];
cfg.latency = [0.600 0.900];

stat_gamma = ft_freqstatistics(cfg, temp{:});

%% STATS ANALYTIC

tests = {[1 4] [2 4] [3 4]};
n_tests = length(tests);
stats = cell(1, n_tests);

for test_index = 1:n_tests
    comparison = tests{test_index};
    temp = [stat_cell(:, comparison(1)), stat_cell(:, comparison(2))];

    cfg = [];
    cfg.method = 'analytic';
    cfg.correctm = 'no';
    cfg.alpha = 0.05;
    cfg.statistic = 'depsamplesT';
    cfg.design(1, :) = [1:n_subjects 1:n_subjects];
    cfg.design(2, :) = [ones(1, n_subjects), 2*ones(1, n_subjects)];
    cfg.uvar = 1;
    cfg.ivar = 2;

    stats{test_index} = ft_freqstatistics(cfg, temp{:});
    
end

%% STATS PLOT

% GAs{1}.mask = stats{1}.mask;
% GAs{1}.mask = stat.mask;

cfg = [];
cfg.frequency = [2 30];
cfg.latency = [0.400 1.000];

% cfg.frequency = [40 100];
% cfg.latency = [0.600 0.900];

temp = ft_selectdata(cfg, GAs{1});
temp.mask = stat.mask;
% temp.mask = stat_gamma.mask;


cfg = [];
% cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [2 30];
cfg.xlim = [0.400 1.000];
% cfg.ylim = [40 100];
% cfg.xlim = [0.600 0.900];
% cfg.zlim = [-0.100 0.100];

figure
ft_multiplotTFR(cfg, temp);
% ft_multiplotTFR(cfg, stats{1});
%% HACK OF STAT

temp = GAs{1};%stats{1};
temp.mask = stats{1}.stat > 2.04 | stats{1}.stat < -2.04;

cfg = [];
% cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
cfg.layout = 'neuromag306cmb.lay';
cfg.ylim = [2 40];
% cfg.zlim = [-3 3];

figure
ft_multiplotTFR(cfg, temp);