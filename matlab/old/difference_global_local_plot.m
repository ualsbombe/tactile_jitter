%%  CLEAR AND ADDPATHS

clear variables
restoredefaultpath

sd = 'NatMEG_0245/190809/';
% sd = 'NatMEG_0177/190619/';

raw_path = ['/archive/20068_tactile_prediction/MEG/' sd];
data_path = ['/home/lau/analyses/tactile_prediction_jitter/data/' ...
            sd];
fieldtrip_path = '/home/lau/matlab/fieldtrip';     
figures_path = ['/home/lau/analyses/tactile_prediction_jitter/figures/' ...
                sd(1:11) '/difference_global_local'];
            
mr_path = fullfile(data_path, '..');

addpath(fieldtrip_path);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)

%% SET TOILIM AND FREQUENCY

toilim = [0.000 0.500];
frequency = 10;

if frequency == 10
    band_string = 'alpha';
elseif frequency == 20
    band_string = 'beta';
else
    error('Frequency not defined')
end

time_string = [num2str(toilim(1)*1e3) '_to_' num2str(toilim(2)*1e3) '_ms'];
band_time_string = [band_string '_' time_string];

%% LOAD

load(fullfile(data_path, 'global_ISI', ['bar_matrices_' band_time_string]))
% load(fullfile(data_path, 'global_ISI', ['source_contrast_ints_' ...
%                                             band_time_string]))
load(fullfile(data_path, 'global_ISI', ['source_contrasts_' band_time_string]))
load(fullfile(data_path, 'global_ISI', 'tfr_cmb_contrasts.mat'))
global_bar_matrices = bar_matrices;
% global_source_contrast_ints = source_contrast_ints;
global_source_contrasts = source_contrasts;
global_tfrs = tfr_cmb_contrasts;

load(fullfile(data_path, 'local_ISI', ['bar_matrices_' band_time_string]))
% load(fullfile(data_path, 'local_ISI', ['source_contrast_ints_' ...
%                                                         band_time_string]))
load(fullfile(data_path, 'local_ISI', ['source_contrasts_' band_time_string]))
load(fullfile(data_path, 'local_ISI', 'tfr_cmb_contrasts.mat'))

local_bar_matrices = bar_matrices;
% local_source_contrast_ints = source_contrast_ints;
local_source_contrasts = source_contrasts;
local_tfrs = tfr_cmb_contrasts;

clear bar_matrices source_contrast_ints source_contrasts tfr_cmb_contrasts

%% DIFFERENCES MATRICES

diff_bar_matrices = struct();
names = fieldnames(global_bar_matrices);
n_names = length(names);

for name_index = 1:n_names
    
    name = names{name_index};
    diff_bar_matrices.(name) = global_bar_matrices.(name) - ...
                               local_bar_matrices.(name);
                           
end

%% DIFFERENCES SOURCES (not already interpolated)

load('standard_mri.mat') %% loads "mri"

n_events = length(global_source_contrasts);
diff_ints = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = 'x1 - x2';
    
    diff = ft_math(cfg, ...
                          global_source_contrasts{event_index}, ...
                          local_source_contrasts{event_index});

    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;  
    
    diff_ints{event_index} = ft_sourceinterpolate(cfg, diff, mri);
                                                  
end

%% TFR DIFFERENCES

n_events = length(global_tfrs);
diff_tfrs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1 - x2';
    
    diff_tfrs{event_index} = ft_math(cfg, global_tfrs{event_index}, ...
                                          local_tfrs{event_index});

end                                      
%% DIFFERENCES SOURCES (already interpolated)
% 
% n_events = length(global_source_contrast_ints);
% diff_ints = cell(1, n_events);
% diff_ints_reverse = cell(1, n_events); % workaround for missing negative ...
%                                        % plotting utility of ft_sourceplot
% 
% for event_index = 1:n_events
%     
%     cfg = [];
%     cfg.parameter = 'pow';
%     cfg.operation = 'x1 - x2';
%     
%     diff_ints{event_index} = ft_math(cfg, ...
%                                   global_source_contrast_ints{event_index}, ...
%                                   local_source_contrast_ints{event_index});
%     diff_ints{event_index}.anatomy = ...
%                 global_source_contrast_ints{event_index}.anatomy;
%     diff_ints{event_index}.coordsys = ...
%                 global_source_contrast_ints{event_index}.coordsys;
%             
%     diff_ints_reverse{event_index} = ft_math(cfg, ...
%                                   local_source_contrast_ints{event_index}, ...
%                                   global_source_contrast_ints{event_index});
%     diff_ints_reverse{event_index}.anatomy = ...
%                 global_source_contrast_ints{event_index}.anatomy;
%     diff_ints_reverse{event_index}.coordsys = ...
%                 global_source_contrast_ints{event_index}.coordsys;
%                                       
% end

%% PLOT DIFF TFRs

cfg = [];
cfg.layout = 'neuromag306cmb.lay';

figure
ft_multiplotTFR(cfg, diff_tfrs{3});

%% label grouping

tissue_groups = [];
tissue_groups.somatomotor = {'Precentral_L' 'Postcentral_L' ...
                             'Supp_Motor_Area_L' ...        
                             'Precentral_R' 'Postcentral_R' ...
                             'Supp_Motor_Area_R' ...
                             };

tissue_groups.cerebellum  = {'Cerebellum_Crus1_L' 'Cerebellum_Crus2_L' ...
                             'Cerebellum_3_L' 'Cerebellum_4_5_L' ...
                             'Cerebellum_6_L' 'Cerebellum_7b_L' ...
                             'Cerebellum_8_L' 'Cerebellum_9_L' ...
                             'Cerebellum_10_L' ...
                             'Cerebellum_Crus1_R' 'Cerebellum_Crus2_R' ...
                             'Cerebellum_3_R' 'Cerebellum_4_5_R' ...
                             'Cerebellum_6_R' 'Cerebellum_7b_R' ...
                             'Cerebellum_8_R' 'Cerebellum_9_R' ...
                             'Cerebellum_10_R'
                             };
                         
tissue_groups.vermis      = {'Vermis_1_2' 'Vermis_3' 'Vermis_4_5' 'Vermis_6' ...
                             'Vermis_7' 'Vermis_8' 'Vermis_9' 'Vermis_10'};

tissue_groups.frontal     = {'Frontal_Sup_L' 'Frontal_Sup_Orb_L' ...
                             'Frontal_Mid_L' 'Frontal_Mid_Orb_L' ...
                             'Frontal_Inf_Oper_L' 'Frontal_Inf_Tri_L' ...
                             'Frontal_Inf_Orb_L' 'Rolandic_Oper_L' ...
                             'Olfactory_L' 'Frontal_Sup_Medial_L' ...
                             'Frontal_Med_Orb_L' 'Rectus_L' ...
                             'Frontal_Sup_R' 'Frontal_Sup_Orb_R' ...
                             'Frontal_Mid_R' 'Frontal_Mid_Orb_R' ...
                             'Frontal_Inf_Oper_R' 'Frontal_Inf_Tri_R' ...
                             'Frontal_Inf_Orb_R' 'Rolandic_Oper_R' ...
                             'Olfactory_R' 'Frontal_Sup_Medial_R' ...
                             'Frontal_Med_Orb_R' 'Rectus_R'};

tissue_groups.hippocampus = {'Insula_L' 'Cingulum_Ant_L' 'Cingulum_Mid_L' ...
                             'Cingulum_Post_L' 'Hippocampus_L' ...
                             'ParaHippocampal_L' 'Amygdala_L' ...
                             'Insula_R' 'Cingulum_Ant_R' 'Cingulum_Mid_R' ...
                             'Cingulum_Post_R' 'Hippocampus_R' ...
                             'ParaHippocampal_R' 'Amygdala_R'};
                         
tissue_groups.occipital   = {'Calcarine_L' 'Cuneus_L' 'Lingual_L' ...
                             'Occipital_Sup_L' 'Occipital_Mid_L' ...
                             'Occipital_Inf_L' 'Fusiform_L' ...
                             'Calcarine_R' 'Cuneus_R' 'Lingual_R' ...
                             'Occipital_Sup_R' 'Occipital_Mid_R' ...
                             'Occipital_Inf_R' 'Fusiform_R'};
                         
tissue_groups.parietal    = {'Parietal_Sup_L' 'Parietal_Inf_L' ...
                             'SupraMarginal_L' 'Angular_L' 'Precuneus_L' ...
                             'Paracentral_Lobule_L' ...
                             'Parietal_Sup_R' 'Parietal_Inf_R' ...
                             'SupraMarginal_R' 'Angular_R' 'Precuneus_R' ...
                             'Paracentral_Lobule_R'};
                         
tissue_groups.temporal    = {'Heschl_L' 'Temporal_Sup_L' ...
                             'Temporal_Pole_Sup_L' 'Temporal_Mid_L' ...
                             'Temporal_Pole_Mid_L' 'Temporal_Inf_L' ...
                             'Heschl_R' 'Temporal_Sup_R' ...
                             'Temporal_Pole_Sup_R' 'Temporal_Mid_R' ...
                             'Temporal_Pole_Mid_R' 'Temporal_Inf_R'};
                         
tissue_groups.midbrain    = {'Caudate_L' 'Putamen_L' 'Pallidum_L' ...
                             'Thalamus_L' ...
                             'Caudate_R' 'Putamen_R' 'Pallidum_R' ...
                             'Thalamus_R'};


%% BARPLOT

close all

names = fieldnames(global_bar_matrices);
n_names = length(names);

for name_index = 1:n_names
    
    name = names{name_index};
    diff_matrix = diff_bar_matrices.(name)';
%     std_matrix = std_bar_matrices.(name)';
%     upper_matrix = bar_matrix + std_matrix;
%     lower_matrix = bar_matrix - std_matrix;
    
    n_groups = size(diff_matrix, 1);
    n_bars = size(diff_matrix, 2);
    group_width = min(0.8, n_bars/(n_bars + 1.5));

    labels = tissue_groups.(name);
    n_labels = length(labels);

    mplot = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    hold on
    
    bar(1:n_labels, diff_matrix);
%     for bar_index = 1:n_bars
%         x = (1:n_groups) - group_width/2 + ...
%                     (2*bar_index-1) * group_width / (2*n_bars);
%         errorbar(x, diff_matrix(:, bar_index), std_matrix(:, bar_index), 'k.')
%     end
    legend({'0%' '5%' '15%'})
    set(gca, 'xtick', 1:n_labels)
    set(gca, 'xticklabel', labels)
    set(gca, 'XTickLabelRotation', 90)
    set(gca, 'TickLabelInterpreter', 'none')
    title([band_string ' power: Difference global minus local ISI'])
%     ylabel([upper(parcellation.cfg.method(1)) parcellation.cfg.method(2:end) ...
%            ' Power increase'])
    ylabel('Power increase')
    ylim([-0.10 0.10])
    xlabel('ROI')
    
    print(fullfile(figures_path, band_string, ...
                    [name '.jpg']), '-djpeg', '-r300');
    
end

%% SOURCEPLOT DIFFERENCE (global minus local)

atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));

names = {'0%' '5%' '15%'};
                             
cfg = [];
cfg.funparameter = 'pow';
cfg.funcolorlim = [0.00 0.10];
cfg.opacitylim = cfg.funcolorlim;
cfg.atlas = atlas;
cfg.roi = tissue_groups.vermis;

for contrast_index = 1:3

    name = names{contrast_index};
    ft_sourceplot(cfg, diff_ints{contrast_index});
    title(['Jitter: ' name])
%     print(fullfile(figures_path, ['sourceplot_' name '.jpg']), ...
%                                     '-djpeg', '-r300');

end

% %% SOURCEPLOT DIFFERENCE (reverse WORKAROUND) (local minus global)
% 
% atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
%                                  'template/atlas/aal/ROI_MNI_V4.nii'));
% 
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.funcolorlim = [0.03 0.10];
% cfg.atlas = atlas;
% 
% ft_sourceplot(cfg, diff_ints_reverse{3});