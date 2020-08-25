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

n_subjects = length(subjects); 

%% LOAD PATHS

spacing = '7.5';
contrast = {'o0' 'o15'};
h_freq = '14';
l_freq = '30';
input = 'tactile_jitter_parcel.mat';
output = 'tactile_jitter_parcel.png';
reg = '0.0';
weight_norm = 'unit-noise-gain';
channel_type = 'grad';
tmin = -0.400;
tmax = 1.100;
ISI_type = 'normal';

paths = get_beamformer_paths(spacing, contrast, h_freq, l_freq, input, ...
                             output, reg, weight_norm, channel_type, ...
                             tmin, tmax, ISI_type);
load_paths = paths(1, 1:2);

%% PLOT SUBJECTS

close all
n_times = 1501;
time = linspace(-0.400, 1.100, n_times);
% n_labels = length(parcel_time_series.label);
% ylims = [min(min(parcel_time_series.pow)), max(max(parcel_time_series.pow))];
xlims = [-0.400 1.100];
n_subjects = 30;
comparisons = zeros(n_subjects, n_labels, n_times);
condition_1 = comparisons;
condition_2 = comparisons;
    
for subject_index = 1:n_subjects
    figure_index = 0;
    subject = subjects{subject_index};
    date = [dates{subject_index} '_000000'];
    path_1 = fullfile(scratch_path, subject, date, 'hilbert', ...
                      load_paths{1});
    path_2 = fullfile(scratch_path, subject, date, 'hilbert', ...
                      load_paths{2});
    disp(['Loading data for subject: ' subject]);
    data_1 = load(path_1);
    data_2 = load(path_2);
    n_labels = length(data_1.label);
    if subject_index == 1
        comparisons = zeros(n_subjects, n_labels, n_times);
    end


%     comparison = log10(data_1.pow ./ data_2.pow);
    comparison = (data_1.pow - data_2.pow) ./ data_2.pow;
    comparisons(subject_index, :, :) = comparison;
    condition_1(subject_index, :, :) = data_1.pow;
    condition_2(subject_index, :, :) = data_2.pow;
    for label_index = 1:n_labels
        
        plot_index = mod(label_index, 20);
        if plot_index == 1
            figure_index = figure_index + 1;
            if subject_index == 1
                figure('units', 'normalized', ...
                       'outerposition', [0 0 1 1], 'visible', 'off');
            else
                figure(figure_index)
                if subject_index == n_subjects
                    set(gcf, 'visible', 'on')
                end
            end
        end
        if plot_index == 0
            plot_index = 20;
        end

        subplot(4, 5, plot_index)
        hold on
        plot(time, comparison(label_index, :));
        if subject_index == 1
            title(data_1.label{label_index}, 'interpreter', 'none')

        elseif subject_index == n_subjects
            average = squeeze(mean(comparisons, 1));
            plot(time, average(label_index, :), 'b', 'linewidth', 8)
            this_ylim = ylim;
            plot([0 0], [this_ylim(1) this_ylim(2)], 'k--', 'linewidth', 3);
            xlim(xlims)
            plot([xlims(1) xlims(2)], [0 0], 'k--', 'linewidth', 3);
            ylim(this_ylim)
            
        end
        
    end

end


%     ylim(ylims)

%% PLOT GRAND AVERAGE THRESHOLDED BY STATS

close all
parcel_method = 'max';

path = ['/home/lau/projects/cerebellar_clock/scratch/tactile_jitter' ...
        '/MEG/grand_averages/statistics/hilbert' ...
        '/stcs/common_filter/contrasts'];
filename = [parcel_method '_t-test_p_0.05_normal_ISI_mag_reg_0.0_' ...
            'weight_unit-noise-gain_o0_o15_contrast_7.5_mm_' ...
            'hp_14_Hz_lp_30_Hz_-400_ms_to_1100_ms_' ...
            'tactile_jitter_parcel.mat'];

data = load(fullfile(path, filename));
n_labels = length(data.label);
n_times = 1501;
time = linspace(-0.400, 1.100, n_times);
ylims = [0 max(max(data.pow))];


for label_index = 1:n_labels
        
        plot_index = mod(label_index, 20);
        if plot_index == 1
            figure('units', 'normalized', ...
                   'outerposition', [0 0 1 1], 'visible', 'on');
     
        end
        if plot_index == 0
            plot_index = 20;
        end

        subplot(4, 5, plot_index)
        hold on
        plot(time, data.pow(label_index, :));
        title([parcel_method ' ' data.label{label_index}], ...
            'interpreter', 'none')
        ylim(ylims)
        xlim([time(1) time(end)])
        
end

%% A PRIORI AREAS

parcels = {
            'Insula_L'
            'Insula_R'
            'Cingulum_Mid_L'
            'Cingulum_Mid_R'
            'Thalamus_L'
            'Thalamus_R'
            'Putamen_L'
            'Putamen_R'
            'Postcentral_L'
            'Postcentral_R'
            'Parietal_Inf_L'
            'Parietal_Inf_R'
            };
        
title_names = {
                'Left Insula'
                'Right Insula'
                'Left Middle Cingulate Cortex'
                'Right Middle Cingulate Cortex'
                'Left Thalamus'
                'Right Thalamus'
                'Left Putamen'
                'Right Putamen'
                'Left SI'
                'Right SI'
                'Left Inferior Parietal Cortex'
                'Right Inferior Parietal Cortex'
                };
            
cerebellar_parcels = {
                       'Cerebellum_Crus1_L'
                       'Cerebellum_Crus1_R'
                       'Vermis_1_2'
                       'Cerebellum_Crus2_L'
                       'Cerebellum_Crus2_R'
                       'Vermis_3'
                       'Cerebellum_3_L'
                       'Cerebellum_3_R'
                       'Vermis_4_5'
                       'Cerebellum_4_5_L'
                       'Cerebellum_4_5_R'
                       'Vermis_6'
                       'Cerebellum_6_L'
                       'Cerebellum_6_R'
                       'Vermis_7'
                       'Cerebellum_7b_L'
                       'Cerebellum_7b_R'
                       'Vermis_8'
                       'Cerebellum_8_L'
                       'Cerebellum_8_R'
                       'Vermis_9'
                       'Cerebellum_9_L'
                       'Cerebellum_9_R'
                       'Vermis_10'
                       'Cerebellum_10_L'
                       'Cerebellum_10_R'
                     };
                 
cerebellar_title_names = {
                            'Left Cerebellum Crus 1'
                            'Right Cerebellum Crus 1'
                            'Vermis 1 & 2'
                            'Left Cerebellum Crus 2'
                            'Right Cerebellum Crus 2'
                            'Vermis 3'
                            'Left Cerebellum 3'
                            'Right Cerebellum 3'
                            'Vermis 4 & 5'
                            'Left Cerebellum 4 & 5'
                            'Right Cerebellum 4 & 5'
                            'Vermis 6'
                            'Left Cerebellum 6'
                            'Right Cerebellum 6'
                            'Vermis 7'
                            'Left Cerebellum 7b'
                            'Right Cerebellum 7b'
                            'Vermis 8'
                            'Left Cerebellum 8'
                            'Right Cerebellum 8'
                            'Vermis 9'
                            'Left Cerebellum 9'
                            'Right Cerebellum 9'
                            'Vermis 10'
                            'Left Cerebellum 10'
                            'Right Cerebellum 10'
                            };
            
            
%% PLOT A PRIORI AREAS
            
close all
parcel_method = 'median';

path = ['/home/lau/projects/cerebellar_clock/scratch/tactile_jitter' ...
        '/MEG/grand_averages/statistics/hilbert' ...
        '/stcs/common_filter/contrasts'];
filename = [parcel_method '_t-test_p_0.05_normal_ISI_mag_reg_0.0_' ...
            'weight_unit-noise-gain_o0_o15_contrast_7.5_mm_' ...
            'hp_14_Hz_lp_30_Hz_-400_ms_to_1100_ms_' ...
            'tactile_jitter_parcel.mat'];
        
figure_path = ['/home/lau/projects/cerebellar_clock/scratch' ...
                '/tactile_jitter' ...
               '/figures/grand_averages/hilbert' ...
               '/stcs'];

data = load(fullfile(path, filename));
n_times = 1501;
time = linspace(-0.400, 1.100, n_times);
n_parcels = length(parcels);
n_cerebellar_parcels = length(cerebellar_parcels);
h = figure('units', 'normalized', 'outerposition', [0 0 0.5 1]);
ylims = [0 max(max(data.pow))];

set(0, 'defaultaxeslinewidth', 3, 'defaultlinelinewidth', 3, ...
       'defaultaxesfontsize', 12, 'defaultaxesfontweight', 'bold')

for parcel_index = 1:n_parcels
    
    parcel = parcels{parcel_index};
    data_index = strcmp(parcel, data.label);
    subplot(6, 2, parcel_index);
    hold on
    plot(time * 1e3, data.pow(data_index, :));
    plot([0 0], ylims, 'k--')
    title(title_names{parcel_index})
    ylim(ylims)
    xlim([-400 400])
    if parcel_index > 10
        xlabel('Time (ms)')
    else
        set(gca, 'xtick', [])
    end
    if parcel_index == 1
        ylabel('Rel. increase')
    end
    if mod(parcel_index, 2) ~= 1
        set(gca, 'ytick', [])
    end
    
end

text(-700, 0.155, 'A priori areas', 'clipping', 'off', ...
    'fontsize', 20, 'fontweight', 'bold');

filename = 'o0_o15_beta_band_a_priori_areas.png';
% print(fullfile(figure_path, filename), '-dpng');

h = figure('units', 'normalized', 'outerposition', [0.5 0 0.5 1]);

for cerebellar_parcel_index = 1:n_cerebellar_parcels
    
    cerebellar_parcel = cerebellar_parcels{cerebellar_parcel_index};
    data_index = strcmp(cerebellar_parcel, data.label);
    subplot(9, 3, cerebellar_parcel_index);
    pos = get(gca, 'Position');
    pos(3) = 0.2;
    set(gca, 'Position', pos)
    hold on
    plot(time * 1e3, data.pow(data_index, :));
    plot([0 0], ylims, 'k--')
    title(cerebellar_title_names{cerebellar_parcel_index})
    ylim(ylims)
    xlim([-400 400])
    if cerebellar_parcel_index > 23
        xlabel('Time (ms)')
    else
        set(gca, 'xtick', [])
    end
    if cerebellar_parcel_index == 1
        ylabel('Rel. increase')
    end
    if mod(cerebellar_parcel_index, 3) ~= 1
        set(gca, 'ytick', [])
    end
    
end

text(-400, 0.225, 'Cerebellar areas', 'clipping', 'off', ...
    'fontsize', 20, 'fontweight', 'bold');

filename = 'o0_o15_beta_band_cerebellar_areas.png';
% print(fullfile(figure_path, filename), '-dpng');