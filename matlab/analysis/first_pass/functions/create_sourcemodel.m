function [dummy] = create_sourcemodel(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
fieldtrip_path = input.fieldtrip_path;
subject = input.subject;
meg_date = input.meg_date;
mr_date = input.mr_date;
input_file = input.input_file;
template = input.template;
output_file = input.output_file;
overwrite = input.overwrite;

%% LOAD FILES

mr_path = fullfile(load_path, subject, mr_date);

mri_resliced = load(fullfile(mr_path, input_file));
% load sourcemodel
load(fullfile(fieldtrip_path, 'template', 'sourcemodel', template));


%% WARPED SOURCE MODEL

cfg = [];
cfg.warpmni = 'yes';
cfg.template = sourcemodel;
cfg.nonlinear = 'yes';
cfg.mri = mri_resliced;
cfg.unit = 'm';

warped_sourcemodel = ft_prepare_sourcemodel(cfg);

%% SAVE

full_meg_date = [meg_date '_000000'];
fullpath = fullfile(load_path, subject, full_meg_date, output_file);

save_struct(warped_sourcemodel, fullpath, overwrite);