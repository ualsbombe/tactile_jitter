function [dummy] = check_realignment(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

raw_path = input.raw_path;
load_path = input.load_path;
subject = input.subject;
date = input.date;
meg_date = input.meg_date;
input_file = input.input_file;


%% LOAD FILE

subject_path = fullfile(load_path, subject, date);
mri_check = load(fullfile(subject_path, input_file));

%% GET HEAD POINTS

full_meg_date = [meg_date '_000000'];
head_points_path = fullfile(raw_path, subject, full_meg_date, 'MEG');

dir_struct = dir(head_points_path);
n_folders = length(dir_struct);

for folder_index = 1:n_folders
    folder_name = dir_struct(folder_index).name;
    if startsWith(folder_name, '001')
        relevant_folder = folder_name;
    end
end

fullpath = fullfile(head_points_path, relevant_folder, 'files');
dir_struct = dir(fullpath);
first_file_index = 3;
filename = dir_struct(first_file_index).name;

headshape = ft_read_headshape(fullfile(fullpath, filename));

%% REALIGN

cfg = [];
cfg.method = 'headshape';
cfg.headshape.headshape = headshape;
cfg.coordsys = 'neuromag';

ft_volumerealign(cfg, mri_check);
