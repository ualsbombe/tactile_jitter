function [dummy] = check_sourcemodel(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

raw_path = input.raw_path;
load_path = input.load_path;
subject = input.subject;
meg_date = input.meg_date;
mr_date = input.mr_date;
input_files = input.input_files;


%% LOAD FILES
full_meg_date = [meg_date '_000000'];

mr_path = fullfile(load_path, subject, mr_date);
meg_path = fullfile(load_path, subject, full_meg_date);

mesh_scalp = load(fullfile(mr_path, input_files{1}));
warped_sourcemodel = load(fullfile(meg_path, input_files{2}));

%% GET SENSORS

sensor_path = fullfile(raw_path, subject, full_meg_date, 'MEG');

dir_struct = dir(sensor_path);
n_folders = length(dir_struct);

for folder_index = 1:n_folders
    folder_name = dir_struct(folder_index).name;
    if startsWith(folder_name, '001')
        relevant_folder = folder_name;
    end
end

fullpath = fullfile(sensor_path, relevant_folder, 'files');
dir_struct = dir(fullpath);
first_file_index = 3;
filename = dir_struct(first_file_index).name;

sensors = ft_read_sens(fullfile(fullpath, filename));
sensors = ft_convert_units(sensors, 'm');

%% PLOT

figure
hold on
ft_plot_mesh(mesh_scalp)
ft_plot_sens(sensors)
ft_plot_mesh(warped_sourcemodel)