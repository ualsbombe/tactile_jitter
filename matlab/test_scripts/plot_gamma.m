restoredefaultpath

addpath /home/lau/matlab/fieldtrip/
ft_defaults

%% COMBINE TFRS

load('gamma_tfrs')


n_events = length(gamma_tfrs);
tfr_cmbs = cell(1, n_events);

for event_index = 1:n_events
    
    cfg = [];
    
    tfr_cmbs{event_index} = ft_combineplanar(cfg, gamma_tfrs{event_index});
    
end

%% TFR CONTRASTS

contrasts = {[1 4] [2 4] [3 4]};
n_contrasts = length(contrasts);

tfr_contrasts = cell(1, n_contrasts);
tfr_cmb_contrasts = cell(1, n_contrasts);

for contrast_index = 1:n_contrasts
    
    contrast = contrasts{contrast_index};
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2) / (x1+x2)';
    
    tfr_contrasts{contrast_index} = ft_math(cfg, ...
                                               gamma_tfrs{contrast(1)}, ...
                                               gamma_tfrs{contrast(2)});
                                           
    cfg = [];
    
    tfr_cmb_contrasts{contrast_index} = ft_combineplanar(cfg, ...
                                            tfr_contrasts{contrast_index});
                                 
end

% %% PLOT TFR
% 
% close all
% 
% cfg = [];
% cfg.layout = 'neuromag306cmb.lay';
% % cfg.xlim = [-1.200 1.200];
% cfg.zlim = [-0.250 0.250];
% cfg.ylim = [2 40];
% % cfg.ylim = [40 100];
% 
% 
% n_contrasts = length(tfr_cmb_contrasts);
% 
% for index = 1:1%n_contrasts
%     figure
%     ft_multiplotTFR(cfg, tfr_cmb_contrasts{index});
% end

%% SINGLE PLOT

close all

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
% cfg.xlim = [-1.200 1.200];
% cfg.zlim = [-0.500 0.500];
cfg.ylim = [40 100];
% cfg.channel = 'MEG0412+0413';
cfg.channel = 'MEG1322+1323';


n_contrasts = length(tfr_cmb_contrasts);

for index = 1:3%n_contrasts
    figure
    ft_singleplotTFR(cfg, tfr_cmb_contrasts{index});
end
