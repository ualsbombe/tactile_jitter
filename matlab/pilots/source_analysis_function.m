
function [] = source_analysis_function(data_path, toilim, frequency, ...
                                       tapsmofrq)

band_string = [num2str(frequency) '_Hz'];
time_string = [num2str(toilim(1)*1e3) '_to_' num2str(toilim(2)*1e3) '_ms'];

%% CHOOSE DATA

load(fullfile(data_path, 'cleaned_data.mat'));

cfg = [];
cfg.toilim = toilim;

data_t_win = ft_redefinetrial(cfg, cleaned_data);

%% POW AND CSD

events = {18 92 102 256};
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
    cfg.foi = frequency;
    cfg.tapsmofrq = tapsmofrq;
    cfg.pad = 'nextpow2';

    powcsds{event_index} = ft_freqanalysis(cfg, data_t_win);
    
end
%% WARPED LEADFIELD

% load(fullfile(data_path, 'sourcemodel_warp.mat'))
% load(fullfile(data_path, 'headmodel.mat'))
% grad = cleaned_data.grad;
% grad = ft_convert_units(grad, 'm');
% 
% cfg = [];
% cfg.sourcemodel = sourcemodel_warp;    %% where are the sources?
% cfg.headmodel   = headmodel;      %% how do currents spread?
% cfg.grad        = grad; %% where are the sensors?
% cfg.channel = 'MEGGRAD';
% cfg.normalize = 'yes';
% 
% % how do sources and sensors connect?
% warped_leadfield = ft_prepare_leadfield(cfg, cleaned_data);
% 
% save(fullfile(data_path, 'warped_leadfield'), 'warped_leadfield')
% save(fullfile(data_path, 'grad'), 'grad');

% %% QUALITY PLOT
% 
% figure;
% hold on
% ft_plot_mesh(sourcemodel_warp, 'vertexcolor', 'red')
% ft_plot_headmodel(headmodel)
% ft_plot_sens(grad);

%% COMMON FILTER ANALYSIS

load(fullfile(data_path, 'headmodel.mat'))
load(fullfile(data_path, 'warped_leadfield.mat'))
load(fullfile(data_path, 'grad.mat'))

% pow

cfg = [];
cfg.channel = 'MEGGRAD';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'powandcsd';
cfg.keeptrials = 'no';
cfg.foi = frequency;
cfg.tapsmofrq = 4;
cfg.pad = 'nextpow2';

powcsd_all = ft_freqanalysis(cfg, data_t_win);

% common filter

cfg = [];
cfg.method = 'dics';
cfg.sourcemodel = warped_leadfield;
cfg.headmodel = headmodel;
cfg.channel = 'MEGGRAD';
cfg.frequency = frequency;
cfg.grad = grad;
cfg.dics.keepfilter = 'yes';

source_all = ft_sourceanalysis(cfg, powcsd_all);

%% SOURCE ANALYSIS (common filter)         

load(['/home/lau/matlab/fieldtrip/template/' ...
    'sourcemodel/standard_sourcemodel3d10mm']); % loads "sourcemodel"
load(fullfile(data_path, 'grad.mat'))

events = {18 92 102 256};
n_events = length(events);

sources_common = cell(1, n_events);

for event_index = 1:n_events

    cfg = [];
    cfg.method = 'dics';
    cfg.sourcemodel = warped_leadfield;
    cfg.headmodel = headmodel;
    cfg.channel = 'MEGGRAD';
    cfg.frequency = frequency;
    cfg.grad = grad;
    cfg.sourcemodel.filter = source_all.avg.filter;

    sources_common{event_index} = ft_sourceanalysis(cfg, ...
                                                    powcsds{event_index});
    sources_common{event_index}.pos = sourcemodel.pos; % template
    
end

%% SOURCE CONTRASTS

contrasts = {[1 4] [2 4] [3 4] [1 2] [1 3] [2 3]};
n_contrasts = 3; %length(contrasts);

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

%% INTERPOLATE CONTRASTS TO TEMPLATE BRAIN

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

% %% PLOT INTERPOLATED CONTRASTS ON TEMPLATE BRAIN
% 
% % close all
% 
% atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
%                                  'template/atlas/aal/ROI_MNI_V4.nii'));
%                              
% names = {'0%' '5%' '15%'};
% 
% for contrast_index = 1:1
%     name = names{contrast_index};
% 
%     cfg = [];
%     cfg.funparameter = 'pow';
% %     cfg.location = [40 -71 -42];
%     cfg.funcolorlim = [0.035 0.15];
%     cfg.opacitylime = cfg.funcolorlim;
%     cfg.atlas = atlas;
% %     cfg.roi = {'Cerebellum_10_R'}
% 
%     ft_sourceplot(cfg, source_contrast_ints{contrast_index});
%     title(['Jitter: ' name])
% %     print(fullfile(figures_path, band_string, ['sourceplot_' name '_' ...
% %                                                 time_string '.jpg']), ...
% %                                     '-djpeg', '-r300');
%                                 
% end

end