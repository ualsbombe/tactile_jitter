function [] = gamma_oscillations_function(data_path)

%% CHOOSE DATA

load(fullfile(data_path, 'cleaned_data.mat'));


%% GAMMA OSCILLATIONS

events = {18 92 102 256};
n_events = length(events);

gamma_tfrs = cell(1, n_events);

for event_index = 1:n_events

    event = events{event_index};
    
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'MEG';
    cfg.method = 'mtmconvol';
    cfg.toi = -1.500:0.020:1.500;
    cfg.foi = 40:2:100;
    cfg.t_ftimwin = 7 ./ cfg.foi;
    cfg.tapsmofrq = 0.2 .* cfg.foi;
    cfg.trials = find(cleaned_data.trialinfo == event);
    cfg.pad = 'nextpow2';
   
    gamma_tfrs{event_index} = ft_freqanalysis(cfg, cleaned_data);
    
end

save(fullfile(data_path, 'gamma_tfrs'), 'gamma_tfrs');