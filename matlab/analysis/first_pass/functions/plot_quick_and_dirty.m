function [dummy] = plot_quick_and_dirty(input)

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
subject = input.subject;
date = input.date;
input_file = input.input_file;
figures_path = input.figures_path;

%% GET PATH

full_date = [date '_000000'];
subject_path = fullfile(load_path, subject, full_date);
load(fullfile(subject_path, input_file)); %% loads tfr contrast


%% PLOT

close all

h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(1, 2, 1)

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.zlim = [0.7 1.3];
cfg.channel = 'MEG2542+2543';

ft_singleplotTFR(cfg, tfr_contrast);

subplot(1, 2, 2)

cfg = [];
cfg.layout = 'neuromag306cmb.lay';
cfg.channel = 'MEG2542+2543';

ft_singleplotTFR(cfg, tfr_contrast);

%% OUTPUT

fullpath = fullfile(figures_path, subject, full_date, ...
                    'quick_and_dirty.png');

print(h, fullpath, '-dpng', '-r300');
        