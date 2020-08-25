function [dummy] = reslice_segment_and_create_headmodel(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% GET INPUT STRUCT

dummy = 'just_a_dummy_variable';

load_path = input.load_path;
subject = input.subject;
date = input.date;
input_file = input.input_file;
output_files = input.output_files;
overwrite = input.overwrite;

%% LOAD FILE

subject_path = fullfile(load_path, subject, date);
mri_realigned_head_points = load(fullfile(subject_path, input_file));

mri_realigned_head_points.coordsys = 'neuromag';

%% RESLICE

cfg = [];
cfg.resolution = 1;

mri_resliced = ft_volumereslice(cfg, mri_realigned_head_points);

%% SAVE RESLICED
fullpath = fullfile(load_path, subject, date, output_files{1});
save_struct(mri_resliced, fullpath, overwrite);

%% SEGMENT

cfg = [];
cfg.output = {'brain' 'skull' 'scalp'};

mri_segmented = ft_volumesegment(cfg, mri_resliced);

%% SAVE SEGMENTED
fullpath = fullfile(load_path, subject, date, output_files{2});
save_struct(mri_segmented, fullpath, overwrite);

%% CREATE MESHES

cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, mri_segmented);
mesh_brain = ft_convert_units(mesh_brain, 'm');

cfg.tissue = 'scalp';
cfg.numvertices = 1000;

mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
mesh_scalp = ft_convert_units(mesh_scalp, 'm');

%% SAVE MESHES

fullpath = fullfile(load_path, subject, date, output_files{3});
save_struct(mesh_brain, fullpath, overwrite);

fullpath = fullfile(load_path, subject, date, output_files{4});
save_struct(mesh_scalp, fullpath, overwrite);

%% HEAD MODEL (SINGLE SHELL)

cfg = [];
cfg.method = 'singleshell';

headmodel = ft_prepare_headmodel(cfg, mesh_brain);

%% SAVE RESLICED
fullpath = fullfile(load_path, subject, date, output_files{5});
save_struct(headmodel, fullpath, overwrite);