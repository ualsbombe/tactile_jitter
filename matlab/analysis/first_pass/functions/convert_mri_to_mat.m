function [dummy] = convert_mri_to_mat(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

raw_path = input.raw_path;
subject = input.subject;
date = input.date;
recording = input.recording;
save_path = input.save_path;
output_file = input.output_file;
overwrite = input.overwrite;

%% SUBJECT PATH AND FIND RELEVANT RECORDING

subject_path = fullfile(raw_path, subject, date, 'MR');
dir_struct = dir(subject_path);
n_folders = length(dir_struct);

for folder_index = 1:n_folders
    folder_name = dir_struct(folder_index).name;
    if endsWith(folder_name, recording)
        relevant_folder = folder_name;
    end
end

% find first file
recording_path = fullfile(subject_path, relevant_folder, 'files');
dir_struct = dir(recording_path);
first_file_index = 3; %% 1 and 2 are '.' and '..' respectively

filename = dir_struct(first_file_index).name;
filepath = fullfile(recording_path, filename);

mri = ft_read_mri(filepath);


%% SAVE
fullpath = fullfile(save_path, subject, date, output_file);

save_struct(mri, fullpath, overwrite);