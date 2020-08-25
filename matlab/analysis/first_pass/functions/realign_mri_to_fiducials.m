function [dummy] = realign_mri_to_fiducials(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
subject = input.subject;
date = input.date;
input_file = input.input_file;
output_file = input.output_file;
overwrite = input.overwrite;

%% LOAD FILE

subject_path = fullfile(load_path, subject, date);
mri = load(fullfile(subject_path, input_file));

%% REALIGN

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'neuromag';

mri_realigned_fiducial = ft_volumerealign(cfg, mri);

%% SAVE
fullpath = fullfile(load_path, subject, date, output_file);

save_struct(mri_realigned_fiducial, fullpath, overwrite);