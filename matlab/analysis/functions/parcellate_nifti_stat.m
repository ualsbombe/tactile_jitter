function [dummy] = parcellate_nifti_stat(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dummy = 'just_a_dummy';

%% LOAD
filebegin = [input.stat_function '_p_' num2str(input.p_threshold) '_' ...
             input.ISI_type '_ISI_' ...
             input.channel_type '_reg_' input.reg '_weight_' ...
             input.weight_norm '_' input.contrasts{1} '_' ...
             input.contrasts{2} '_contrast_' input.spacing '_mm_hp_' ...
             num2str(input.h_freq) '_Hz_lp_' num2str(input.l_freq) ...
             '_Hz_' num2str(input.tmin * 1e3) ...
             '_ms_to_' num2str(input.tmax * 1e3) '_ms_'];
load_path = fullfile(input.input_path, [filebegin input.input_file]);
save_path = fullfile(input.output_path, [input.parcel_method '_' ...
                                         filebegin input.output_file]);

if ~exist(save_path, 'file') || input.overwrite

    beamformer = ft_read_mri(load_path);
    template_mri = ft_read_mri(input.mri);
    atlas = ft_read_atlas(input.atlas);

    beamformer.pow = beamformer.anatomy;
    template_mri.coordsys = 'mni';
    atlas.coordsys = 'mni';

    %% INTERPOLATE BEAMFORMER ONTO TEMPLATE MRI

    cfg = [];
    cfg.parameter = 'pow';

    interpolation = ft_sourceinterpolate(cfg, beamformer, template_mri);

    %% CREATE PARCEL TIME COURSES

    n_times = size(interpolation.pow, 4);
    pows = zeros(length(atlas.tissuelabel), n_times);

    for time_index = 1:n_times

        disp(['Running parcellation: ' num2str(time_index) ...
            ' out of: ' num2str(n_times)])
        
        temp = interpolation;
        temp.pow = temp.pow(:, :, :, time_index);
        temp.pos = temp.anatomy;
        atlas.pos = temp.pos; %% think about this

        cfg = [];
        cfg.method = input.parcel_method;
        cfg.parcellation = 'tissue';
        cfg.parameter = {'pow'};

        parcel = ft_sourceparcellate(cfg, temp, atlas);

        pows(:, time_index) = parcel.pow;

    end
    % use last copy of parcel to build full time series
    parcel_time_series = parcel;
    parcel_time_series.pow = pows;
    parcel_time_series.powdimord = 'chan_time';

    %% SAVE

    save(save_path, '-struct', 'parcel_time_series')
else
    disp(['Not overwriting: ' save_path])

end

end

