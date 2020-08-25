function [dummy] = create_leadfield(input)
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
normalize = input.normalize;
output_file = input.output_file;
overwrite = input.overwrite;


%% LOAD FILES
full_meg_date = [meg_date '_000000'];

mr_path = fullfile(load_path, subject, mr_date);
meg_path = fullfile(load_path, subject, full_meg_date);

headmodel = load(fullfile(mr_path, input_files{1}));
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
disp(fullfile(fullpath, filename));

sensors = ft_read_sens(fullfile(fullpath, filename));
sensors = ft_convert_units(sensors, 'm');

%% WARPED LEADFIELD

cfg = [];
cfg.sourcemodel = warped_sourcemodel;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.grad        = sensors; %% where are the sensors?
cfg.channel = 'MEGGRAD';
cfg.normalize = normalize;

% how do sources and sensors connect?
warped_leadfield = ft_prepare_leadfield(cfg);

%% SAVE

fullpath = fullfile(meg_path, output_file);

save_struct(warped_leadfield, fullpath, overwrite);
