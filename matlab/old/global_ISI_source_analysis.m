%%  CLEAR AND ADDPATHS

clear variables
restoredefaultpath

% sd = 'NatMEG_0245/190809/';
sd = 'NatMEG_0177/190619/';
% sd = 'NatMEG_0521/190808/';

mount_path = '/home/lau/mounts/compute';
project_path = '/home/lau/analyses/tactile_prediction_jitter';

raw_path = fullfile(mount_path, 'archive', '20068_tactile_prediction', ...
                    'MEG', sd);
data_path = fullfile(mount_path, project_path, 'data', sd, 'global_ISI');
       
fieldtrip_path = '/home/lau/matlab/fieldtrip';     
figures_path = fullfile(mount_path, project_path,  'figures', ...
                sd(1:11), 'global_ISI');
            
mr_path = fullfile(data_path, '..');

addpath(fieldtrip_path);
ft_defaults

%% SET PLOT DEFAULTS

set(0, 'defaultaxesfontsize', 14, 'defaultaxesfontweight', 'bold', ...
    'defaultlinelinewidth', 3)

%% SET TOILIM AND FREQUENCY

toilim = [-0.150 0.150];
frequency = 20;

if frequency == 10
    band_string = 'alpha';
elseif frequency == 20
    band_string = 'beta';
else
    error('Frequency not defined')
end

time_string = [num2str(toilim(1)*1e3) '_to_' num2str(toilim(2)*1e3) '_ms'];
    
%% CHOOSE DATA

% load(fullfile(data_path, 'preprocessed_data_tfr.mat'));
load(fullfile(data_path, 'cleaned_data.mat'));

cfg = [];
cfg.toilim = toilim;
% cfg.toilim = [0.000 0.500]; % gave rise to right cerebell
% cfg.toilim = [0.000 0.250];% gave rise to S1 and left! cerebell

% data_t_win = ft_redefinetrial(cfg, preprocessed_data_tfr);
data_t_win = ft_redefinetrial(cfg, cleaned_data);

%% POW AND CSD

events = {295 551 1063 2048};
n_events = length(events);

powcsds = cell(1, n_events);

for event_index = 1:n_events
    
    event = events{event_index};

    cfg = [];
    cfg.trials = find(data_t_win.trialinfo == event);
    cfg.channel = 'MEGGRAD';
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'powandcsd';
    cfg.keeptrials = 'no';
%     cfg.foi = 10;
%     cfg.foi = 20;
    cfg.foi = frequency;
    cfg.tapsmofrq = 4;
    cfg.pad = 'nextpow2';

    powcsds{event_index} = ft_freqanalysis(cfg, data_t_win);
    
end

%% warped leadfield

load(fullfile(data_path, 'cleaned_data.mat'))
% load(fullfile(data_path, 'preprocessed_data_tfr.mat'))
load(fullfile(mr_path, 'sourcemodel_warp.mat'))
load(fullfile(mr_path, 'headmodel.mat'))
grad = cleaned_data.grad;
% grad = preprocessed_data_tfr.grad;
grad = ft_convert_units(grad, 'm');

cfg = [];
cfg.sourcemodel = sourcemodel_warp;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
% cfg.grad        = grad; %% where are the sensors?
cfg.channel = 'MEGGRAD';
cfg.normalize = 'yes';

% how do sources and sensors connect?
warped_leadfield = ft_prepare_leadfield(cfg, cleaned_data);
% warped_leadfield = ft_prepare_leadfield(cfg, preprocessed_data_tfr);

save(fullfile(data_path, 'warped_leadfield'), 'warped_leadfield')
save(fullfile(data_path, 'grad'), 'grad');

%% quality plot

figure;
hold on
% ft_plot_mesh(sourcemodel)
ft_plot_mesh(sourcemodel_warp, 'vertexcolor', 'red')
ft_plot_headmodel(headmodel)
ft_plot_sens(grad);

%% run common filter analysis

load(fullfile(mr_path, 'sourcemodel_and_leadfield.mat'))
load(fullfile(mr_path, 'headmodel.mat'))
load(fullfile(data_path, 'warped_leadfield.mat'))
load(fullfile(data_path, 'grad.mat'))

% pow

cfg = [];
cfg.channel = 'MEGGRAD';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'powandcsd';
cfg.keeptrials = 'no';
% cfg.foi = 10;
% cfg.foi = 20;
cfg.foi = frequency;
cfg.tapsmofrq = 4;
cfg.pad = 'nextpow2';

powcsd_all = ft_freqanalysis(cfg, data_t_win);

% common filter

cfg = [];
cfg.method = 'dics';
% cfg.sourcemodel = sourcemodel_and_leadfield;
cfg.sourcemodel = warped_leadfield;
cfg.headmodel = headmodel;
cfg.channel = 'MEGGRAD';
% cfg.frequency = 10;
% cfg.frequency = 20;
cfg.frequency = frequency;
cfg.grad = grad;
cfg.dics.keepfilter = 'yes';

source_all = ft_sourceanalysis(cfg, powcsd_all);

%% SOURCE ANALYSIS (common filter)         

load(['/home/lau/matlab/fieldtrip/template/' ...
    'sourcemodel/standard_sourcemodel3d10mm']); % loads "sourcemodel"
load(fullfile(data_path, 'grad.mat'))

events = {295 551 1063 2048};
n_events = length(events);

sources_common = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'dics';
%     cfg.sourcemodel = sourcemodel_and_leadfield;
    cfg.sourcemodel = warped_leadfield;
    cfg.headmodel = headmodel;
    cfg.channel = 'MEGGRAD';
%     cfg.frequency = 10;
%     cfg.frequency = 20;
    cfg.frequency = frequency;
    cfg.grad = grad;
    cfg.sourcemodel.filter = source_all.avg.filter;

    sources_common{event_index} = ft_sourceanalysis(cfg, powcsds{event_index});
    sources_common{event_index}.pos = sourcemodel.pos; % template
    
end

%% contrasts

contrasts = {[1 4] [2 4] [3 4] [1 2] [1 3] [2 3]};
n_contrasts = length(contrasts);

source_contrasts = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    contrast = contrasts{contrast_index};
    
    cfg = [];
    cfg.parameter = 'pow';
    cfg.operation = '(x1-x2) / (x1+x2)';
    
    source_contrasts{contrast_index} = ft_math(cfg, ...
                                               sources_common{contrast(1)}, ...
                                               sources_common{contrast(2)});
    source_contrasts{contrast_index}.coordsys = 'mni';
                                 
end

savename = fullfile(data_path, ['source_contrasts_' band_string '_' ...
                                time_string]);

save(savename, 'source_contrasts')

%% interpolate contrasts (template brain)

load('standard_mri.mat') %% loads "mri"

n_contrasts = length(source_contrasts);

source_contrast_ints = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts

    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;

    source_contrast_ints{contrast_index} = ft_sourceinterpolate(cfg, ...
                                           source_contrasts{contrast_index}, ...
                                           mri);    
end

savename = fullfile(data_path, ['source_contrast_ints_' band_string '_' ...
                                time_string]);

save(savename, 'source_contrast_ints')

%% plot interpolated contrasts (template brain)

close all

atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));
                             
names = {'0%' '5%' '15%'};

for contrast_index = 1
    name = names{contrast_index};

    cfg = [];
    cfg.funparameter = 'pow';
%     cfg.location = [10 -63 -22];
    cfg.funcolorlim = [0.03 0.20];
    cfg.opacitylime = cfg.funcolorlim;
    cfg.atlas = atlas;
%     cfg.roi = {'Cerebellum_10_R'}

    ft_sourceplot(cfg, source_contrast_ints{contrast_index});
    title(['Jitter: ' name])
%     print(fullfile(figures_path, band_string, ['sourceplot_' name '_' ...
%                                                 time_string '.jpg']), ...
%                                     '-djpeg', '-r300');
                                
end

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

%% parcellate

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

        source_contrasts{contrast_index}.pos = sourcemodel.pos;

        cfg = [];
        cfg.interpmethod = 'nearest';
        cfg.parameter = 'tissue';

        int_atlas = ft_sourceinterpolate(cfg, atlas, ...
                                            source_contrasts{contrast_index});
        int_atlas.pos = sourcemodel.pos;

        cfg = [];
        cfg.parameter = 'pow';
        cfg.parcellation = 'tissue';
        cfg.method = 'mean';

        parcellation = ft_sourceparcellate(cfg, ...
                                source_contrasts{contrast_index}, int_atlas);
        
        cfg.method = 'std';
        
        std_parcellation = ft_sourceparcellate(cfg, ...
                                 source_contrasts{contrast_index}, int_atlas);
                             
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

band_time_string = ['_' band_string '_' time_string];

save(fullfile(data_path, ['bar_matrices' band_time_string]), 'bar_matrices')
save(fullfile(data_path, ['std_bar_matrices' band_time_string]), ...
     'std_bar_matrices')

%% plot

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
    title('Global ISI, Alpha power: Omission relative to non-stimulation')
    ylabel([upper(parcellation.cfg.method(1)) parcellation.cfg.method(2:end) ...
           ' Power increase'])
    ylim([-0.20 0.20])
    xlabel('ROI')
    
%     print(fullfile(figures_path, band_string, [group_name '.jpg']), ...
%                             '-djpeg', '-r300');
end

%% interpolate contrasts (own brain)

load(fullfile(data_path, 'mri_resliced.mat'))

n_contrasts = length(source_contrasts);

source_contrast_ints = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts

    cfg = [];
    cfg.parameter = 'pow';
    cfg.downsample = 2;

    source_contrast_ints{contrast_index} = ft_sourceinterpolate(cfg, ...
                                           source_contrasts{contrast_index}, ...
                                           mri_resliced);    
end
 

%% plot interpolated contrasts (own brain)


cfg = [];
cfg.funparameter = 'pow';
cfg.location = [20 -50 3];
cfg.funcolorlim = [0.05 0.10];
cfg.opacitylime = cfg.funcolorlim;
cfg.maskparameter = 

ft_sourceplot(cfg, source_contrast_ints{1});

%% interpolate contrasts (template brain)

load('standard_mri.mat')



%% plot brain

cfg = [];
cfg.location = [-35 -5 43];

ft_sourceplot(cfg, mri_resliced);

%% plot non-interpolated


cfg = [];
cfg.funparameter = 'pow';

ft_sourceplot(cfg, source_contrasts{1})


%%%%%%%%%%%%%%%%%%%%


%% find peak values

load('standard_mri.mat') %% loads "mri"
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));

labels = {'Cerebellum_6_R' 'Thalamus_R' 'Hippocampus_R' 'Postcentral_L' ...
          'Insula_R'};
n_labels = length(labels);

max_locations = zeros(n_labels, 3);
min_locations = zeros(n_labels, 3);

for label_index = 1:1%n_labels
    
    label = labels{label_index};
    this_contrast = source_contrast_ints{1}; %% based on 1 for now (what is more optimal?) 
       
%     tissue_index = find(strcmp(atlas.tissuelabel, label));
%     this_contrast.pow(atlas.tissue ~= tissue_index) = 0;
    [max_pow, max_index] = max(this_contrast.pow);
    [min_pow, min_index] = min(this_contrast.pow);
    max_locations(label_index, :) = this_contrast.pos(max_index, :);
    min_locations(label_index, :) = this_contrast.pos(min_index, :);
    
end

%% PLOT MAX AND MIN AND SAVE
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));
% for label_index = 1:1%n_labels
%     label = labels{label_index};
%     tissue_index = find(strcmp(atlas.tissuelabel, label));
    for contrast_index = 1:1%3
        this_source = source_contrast_ints{contrast_index};
    
        cfg = [];
        cfg.funparameter = 'pow';
%         cfg.maskparameter = 'mask';
%         cfg.location = min_locations(label_index, :);
        cfg.funcolorlim = [0.05 0.20];
        cfg.opacitylime = cfg.funcolorlim;
        cfg.atlas = atlas;
        cfg.roi = labels;

        ft_sourceplot(cfg, this_source);
    end
% end

%% POW AND CSD

events = {295 551 1063 2048};
n_events = length(events);

powcsds = cell(1, n_events);

for event_index = 1:n_events
    
    event = events{event_index};

    cfg = [];
    cfg.trials = find(data_t_win.trialinfo == event);
    cfg.channel = 'MEGGRAD';
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'powandcsd';
    cfg.keeptrials = 'no';
    cfg.foi = 10;
    cfg.tapsmofrq = 4;

    powcsds{event_index} = ft_freqanalysis(cfg, data_t_win);
    
end

%% SOURCE ANALYSIS

load(fullfile(data_path, 'sourcemodel_and_leadfield.mat'))
load(fullfile(data_path, 'headmodel.mat'))
load(fullfile(data_path, 'warped_leadfield.mat'))

grad = ft_read_sens(fullfile(data_path, 'tactile_jitter-ave.fif'), ...
                   'senstype', 'meg');
grad = ft_convert_units(grad, 'm');               

events = {295 551 1063 2048};
n_events = length(events);

sources = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'dics';
%     cfg.sourcemodel = sourcemodel_and_leadfield;
    cfg.sourcemodel = warped_leadfield;
    cfg.headmodel = headmodel;
    cfg.channel = 'MEGGRAD';
    cfg.frequency = 20;
    cfg.grad = grad;

    sources{event_index} = ft_sourceanalysis(cfg, powcsds{event_index});
%     sources{event_index}.coordsys = 'mni';
    
end

%% plot
% atlas = ft_read_atlas('~/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');

cfg = [];
cfg.funparameter = 'pow';
% cfg.atlas = atlas;

ft_sourceplot(cfg, sources{1});

%% interpolate

load(fullfile(data_path, 'mri_resliced.mat'))

events = {295 551 1063 2048};
n_events = length(events);

source_ints = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.parameter = 'pow';

    source_ints{event_index} = ft_sourceinterpolate(cfg, ...
                                                    sources{event_index}, ...
                                                    mri_resliced);
    
end

%% plot interpolated

cfg = [];
cfg.funparameter = 'pow';
% cfg.funcolorlim = [0.5e-23 2e-23];

ft_sourceplot(cfg, source_ints{4});


%% jan mathijs solution (https://mailman.science.ru.nl/pipermail/fieldtrip/2016-July/010674.html)

atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                                 'template/atlas/aal/ROI_MNI_V4.nii'));
load(fullfile(data_path, 'mri_resliced.mat'))

cfg = [];

normalised_mri = ft_volumenormalise(cfg, mri);

% 

cfg = [];
cfg.parameter = 'tissue';

mri_atlas = ft_sourceinterpolate(cfg, atlas, normalised_mri);

%% add masks

n_contrasts = length(source_contrasts);
labels = {'Cerebellum_6_R'};
n_labels = length(labels);

for contrast_index = 1:n_contrasts
    
    if isfield(source_contrasts{contrast_index}, 'mask')
        source_contrasts{contrast_index} = ...
            rmfield(source_contrasts{contrast_index}, 'mask');
    end
    
    source_contrasts{contrast_index}.mask = ...
        zeros(size(source_contrasts{contrast_index}.pow));
    
    for label_index = 1:n_labels
        
        label = labels{label_index};
        
        source_contrasts{contrast_index}.mask = ...
            source_contrasts{contrast_index}.mask + ...
            mri_atlas.tissue == strcmp(atlas.tissuelabel, label);
    end
end
