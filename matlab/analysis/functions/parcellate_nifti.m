function [dummy] = parcellate_nifti(input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dummy = 'just_a_dummy';

%% LOAD
n_contrasts = length(input.contrasts);
n_freqs = length(input.h_freqs);

template_mri = ft_read_mri(input.mri);
atlas = ft_read_atlas(input.atlas);
template_mri.coordsys = 'mni';
atlas.coordsys = 'mni';

for freq_index = 1:n_freqs
    h_freq = input.h_freqs{freq_index};
    l_freq = input.l_freqs{freq_index};
    disp(['Running frequency range: ' h_freq '-' l_freq ' Hz'])
    for contrast_index = 1:n_contrasts
        contrast = input.contrasts{contrast_index};
        disp(['Running contrast ' contrast{1} ' versus ' contrast{2}])
        paths = get_beamformer_paths(input.spacing, contrast, h_freq, ...
                                     l_freq, input.input_file, ...
                                     input.output_file, input.reg, ...
                                     input.weight_norm, ...
                                     input.channel_type, input.tmin, ...
                                     input.tmax, input.ISI_type);
        n_paths = length(paths);
        for path_index = [1:2 4:5]%n_paths
            load_path = fullfile(input.input_path, paths{1, path_index});
            save_path = fullfile(input.output_path, paths{2, path_index});
            if ~exist(save_path, 'file') || input.overwrite

                beamformer = ft_read_mri(load_path);
                beamformer.pow = beamformer.anatomy;

                %% INTERPOLATE BEAMFORMER ONTO TEMPLATE MRI

                cfg = [];
                cfg.parameter = 'pow';

                interpolation = ft_sourceinterpolate(cfg, beamformer, ...
                                                     template_mri);

                %% CREATE PARCEL TIME COURSES

                n_times = size(interpolation.pow, 4);
                pows = zeros(length(atlas.tissuelabel), n_times);

                for time_index = 1:n_times

                    disp(['Running frequency range: ' h_freq '-' l_freq ...
                          ' Hz'])
                    disp(['Running contrast ' contrast{1} ' versus ' ...
                          contrast{2}])
                    disp(['Running path ' num2str(path_index)])
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
                % preserve memory
                clear interpolation pows parcel_time_series temp ...
                      beamformer parcel
            else
                disp(['Not overwriting: ' save_path])
            end

        end
    end
end

end

