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

subjects = {'0001' '0002' '0003' '0004' '0005'};
dates = {'20190923_000000' '20190930_000000' '20190930_000000' ...
         '20191008_000000' '20191008_000000'};
n_subjects = length(subjects);
n_contrasts = 3;
project_name = 'MINDLAB2019_MEG-CerebellarClock';
fieldtrip_path = '/home/lau/matlab/fieldtrip';     

data_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                          project_name, 'scratch', ...
                          'tactile_jitter', 'MEG');
                      
figures_path = fullfile(home_dir, 'mounts', 'hyades', 'projects', ...
                          project_name, 'scratch', ...
                          'tactile_jitter', 'figures', 'grand_averages');
                    
addpath(fieldtrip_path)
ft_defaults

%% BAND AND TIME STRING

% cerebellum and insula preparation?
% toilim = [-0.700 -0.300];
% frequency = 10;

% cerebellum error?
toilim = [0.600 0.800]; 
frequency = 15;

% insula (L) and somatosensory cortices for regular only
% toilim = [-0.600 -0.400];
% frequency = 60;

% Insula (R) and cerebellum for irregular
% toilim = [0.600 0.900];
% frequency = 80;

band_string = [num2str(frequency) '_Hz'];
time_string = [num2str(toilim(1)*1e3) '_to_' num2str(toilim(2)*1e3) '_ms'];

filename = ['source_contrasts_' band_string '_' time_string '.mat'];

%% COLLECT DATA

load(fullfile(fieldtrip_path, 'template', 'sourcemodel', ...
              'standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel

ga_cell = cell(n_subjects, n_contrasts);

for subject_index = 1:n_subjects

    subject = subjects{subject_index};
    date = dates{subject_index};


    subject_path = fullfile(data_path, subject, date);
                      
                      
    load(fullfile(subject_path, filename))
    for contrast_index = 1:n_contrasts
        temp = source_contrasts{contrast_index};
        temp.pos = template_grid.pos;
        ga_cell{subject_index, contrast_index} = temp;
    end
end

%% GRAND AVERAGES

GAs = cell(1, n_contrasts);
cfg = [];


for contrast_index = 1:n_contrasts
    GAs{contrast_index} = ft_sourcegrandaverage(cfg, ...
                                            ga_cell{:, contrast_index});
                                        
end

%% INTERPOLATING

load('standard_mri.mat') %% loads "mri"

n_contrasts = 3;

GA_beam_ints = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts

    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;

    GA_beam_ints{contrast_index} = ft_sourceinterpolate(cfg, ...
                                           GAs{contrast_index}, ...
                                           mri);    
end

savename = ['GA_beam_ints_' band_string '_' time_string];
save(fullfile(data_path, 'grand_averages', savename), 'GA_beam_ints');

%% LOAD IF NECESSARY

savename = ['GA_beam_ints_' band_string '_' time_string];
load(fullfile(data_path, 'grand_averages', savename))

%% MINUS PLOTS

GA_beam_ints_neg = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    cfg = [];
    cfg.operation = 'x1 * -1';
    cfg.parameter = 'pow';
    
    GA_beam_ints_neg{contrast_index} = ft_math(cfg, ...
                                                GA_beam_ints{contrast_index});
    GA_beam_ints_neg{contrast_index}.anatomy = ...
                   GA_beam_ints{contrast_index}.anatomy;
    GA_beam_ints_neg{contrast_index}.coordsys = ...
                   GA_beam_ints{contrast_index}.coordsys;
    GA_beam_ints_neg{contrast_index}.dim = ...
                   GA_beam_ints{contrast_index}.dim;
    GA_beam_ints_neg{contrast_index}.transform = ...
                   GA_beam_ints{contrast_index}.transform;
    GA_beam_ints_neg{contrast_index}.unit = ...
                   GA_beam_ints{contrast_index}.unit;
    GA_beam_ints_neg{contrast_index}.inside = ...
                   GA_beam_ints{contrast_index}.inside;
                                            
end

%% DEFAULTS

set(0, 'defaultaxesfontweight', 'bold', 'defaultaxesfontsize', 14)

%% SOURCE PLOT

close all

atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));                            
cfg = [];
cfg.funparameter = 'pow';
cfg.atlas = atlas;
cfg.funcolorlim = [0.03 0.06];
cfg.location = [46 -49 -30];
% cfg.opacitylim = cfg.funcolorlim;
% cfg.roi = {'Cerebellum_6_R' 'Cerebellum_7_R' 'Insula_R';};
% cfg.roi = 'Insula_L';
names = {'0%' '5%' '15%'};

for contrast_index = 1:1%n_contrasts
    
    name = names{contrast_index};
    ft_sourceplot(cfg, GA_beam_ints{contrast_index});
%     ft_sourceplot(cfg, GA_beam_ints_neg{contrast_index});
%     print(fullfile(figures_path, ['sourceplot_' name '_' ...
%                                    band_string '_' time_string '.jpg']), ...
%                                     '-djpeg', '-r300');
                                
end                                
    
    
% ft_sourceplot(cfg, GA_beam_ints_neg{1});
% ft_sourceplot(cfg, GA_beam_ints_neg{2});
% ft_sourceplot(cfg, GA_beam_ints_neg{3});

%% PARCELLATE (LABEL GROUPS)

tissue_groups = [];
tissue_groups.somatomotor = {'Precentral_L' 'Postcentral_L' ...
                             'Precentral_R' 'Postcentral_R'
                             };

tissue_groups.cerebellum  = {'Cerebellum_Crus1_L' 'Cerebellum_Crus2_L' ...
                             'Cerebellum_6_L' ...
                             'Cerebellum_Crus1_R' 'Cerebellum_Crus2_R' ...
                             'Cerebellum_6_R'
                             };
                         
tissue_groups.hippocampus = {'Insula_L'... 
                             'Insula_R' };
                                                 
tissue_groups.parietal    = {'Parietal_Sup_L' 'Parietal_Inf_L' ...
                              'Precuneus_L' ...
                             'Parietal_Sup_R' 'Parietal_Inf_R' ...
                             'Precuneus_R' ...
                             };
                                              
tissue_groups.midbrain    = {'Thalamus_L' ...
                             'Thalamus_R'};


%% PARCELLATION

load(['/home/lau/matlab/fieldtrip/template/' ...
    'sourcemodel/standard_sourcemodel3d10mm']); % loads "sourcemodel"
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));
                             
group_names = fieldnames(tissue_groups);
n_group_names = length(group_names);
n_contrasts = 3;
bar_matrices = struct();
std_bar_matrices = struct();

for group_name_index = 1:n_group_names

    group_name = group_names{group_name_index};
    labels = tissue_groups.(group_name);
    n_labels = length(labels);
    bar_matrix = zeros(n_contrasts, n_labels);
    std_bar_matrix = zeros(n_contrasts, n_labels);

    for contrast_index = 1:n_contrasts

        cfg = [];
        cfg.interpmethod = 'nearest';
        cfg.parameter = 'tissue';

        int_atlas = ft_sourceinterpolate(cfg, atlas, ...
                                            GAs{contrast_index});
        int_atlas.pos = sourcemodel.pos;

        cfg = [];
        cfg.parameter = 'pow';
        cfg.parcellation = 'tissue';
        cfg.method = 'mean';

        parcellation = ft_sourceparcellate(cfg, ...
                                GAs{contrast_index}, int_atlas);
        
        cfg.method = 'std';
        
        std_parcellation = ft_sourceparcellate(cfg, ...
                                 GAs{contrast_index}, int_atlas);
                             
        for label_index = 1:n_labels
            label = labels{label_index};
            tissue_index = strcmp(atlas.tissuelabel, label);
            bar_matrix(contrast_index, label_index) = ...
                                                parcellation.pow(tissue_index);
            std_bar_matrix(contrast_index, label_index) = ...
                                             std_parcellation.pow(tissue_index);

        end

                             
            
    end
    
    bar_matrices.(group_name) = bar_matrix;
    std_bar_matrices.(group_name) = std_bar_matrix;
    
end

%% PARCELLATION BAR PLOTS

close all

group_names = fieldnames(tissue_groups);
n_group_names = length(group_names);

for group_name_index = 1:n_group_names
    
    group_name = group_names{group_name_index};
    bar_matrix = bar_matrices.(group_name)';
    std_matrix = std_bar_matrices.(group_name)';
    upper_matrix = bar_matrix + std_matrix;
    lower_matrix = bar_matrix - std_matrix;
    
    n_groups = size(bar_matrix, 1);
    n_bars = size(bar_matrix, 2);
    group_width = min(0.8, n_bars/(n_bars + 1.5));

    labels = tissue_groups.(group_name);
    n_labels = length(labels);

    mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on
    
    bar(1:n_labels, bar_matrix);
    for bar_index = 1:n_bars
        x = (1:n_groups) - group_width/2 + ...
                    (2*bar_index-1) * group_width / (2*n_bars);
        errorbar(x, bar_matrix(:, bar_index), std_matrix(:, bar_index), 'k.')
    end
    legend({'0%' '5%' '15%'})
    set(gca, 'xtick', 1:n_labels)
    set(gca, 'xticklabel', labels)
    set(gca, 'XTickLabelRotation', 90)
    set(gca, 'TickLabelInterpreter', 'none')
    title('Omission relative to non-stimulation')
    ylabel([upper(parcellation.cfg.method(1)) parcellation.cfg.method(2:end) ...
           ' Power increase'])
%     ylim([-0.05 0.05])
    xlabel('ROI')
    
%     print(fullfile(figures_path, band_string, [group_name '.jpg']), ...
%                             '-djpeg', '-r300');
end
