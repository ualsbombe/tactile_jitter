%% INITIALIZE

clear variables
restoredefaultpath

fieldtrip_path = '/home/lau/matlab/fieldtrip/';

addpath(fieldtrip_path)
ft_defaults

figure_path = fullfile('/home/lau/projects/cerebellar_clock/scratch/', ...
            'tactile_jitter', 'figures', 'grand_averages', ...
            'hilbert', 'stcs');


%% SET CHANNEL TYPE

filepath = ['/home/lau/projects/cerebellar_clock/scratch/' ...
            'tactile_jitter/MEG/grand_averages/statistics/hilbert/stcs' ...
            '/common_filter/contrasts'];
filename = ['t-test_p_0.05_normal_ISI_mag_reg_0.0_weight_unit-noise-gain' ...
            '_o0_o15_contrast_7.5_mm_hp_14_Hz_lp_30_Hz' ...
            '_-400_ms_to_400_ms_tactile_jitter_morph-vl.nii'];
        
% filename = ['t-test_p_0.05_normal_ISI_mag_reg_0.0_weight_unit-noise-gain' ...
%             '_s1_s2_contrast_7.5_mm_hp_4_Hz_lp_7_Hz' ...
%             '_-200_ms_to_600_ms_tactile_jitter_morph-vl.nii'];           

%% LOAD MRI AND ATLAS

mri = ft_read_mri(fullfile(fieldtrip_path, ...
                            'template', 'anatomy', 'single_subj_T1.nii'));
                        
mri.coordsys = 'mni';                        
atlas = ft_read_atlas(fullfile(fieldtrip_path, ...
                            'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));
% atlas = ft_read_atlas('/usr/local/fsl/data/atlases/Cerebellum/Cerebellum-MNIflirt-maxprob-thr0-1mm.nii.gz');
atlas.coordsys = 'mni';

beam = ft_read_mri(fullfile(filepath, filename));

%% INTERPOLATE FOR PLOT

temp = beam;
% time_index = 431;
time_index = 400;
temp.anatomy = temp.anatomy(:, :, :, time_index);
temp.dim = temp.dim(1:3);
temp.pow = temp.anatomy;

cfg = [];
cfg.parameter = 'pow';
cfg.downsample = 1;

plot_interp = ft_sourceinterpolate(cfg, temp, mri);

%% PLOT TIME POINT

cfg = [];
cfg.atlas = atlas;
cfg.funparameter = 'pow';
cfg.anaparameter = 'anatomy';
cfg.funcolorlim = [0.000 0.025];
% cfg.funcolorlim = [0.001 ...
%                    max(max(max(plot_interp.pow)))];
cfg.opacitylim = cfg.funcolorlim;
               
cfg.roi = {
%             'Rolandic_Oper_L' %% SII?
%             'Rolandic_Oper_R' %% SII?
%             'Hippocampus_L'
%             'Hippocampus_R'
            'Insula_L'
            'Insula_R'
            'Cingulum_Mid_L'
            'Cingulum_Mid_R'
            'Thalamus_L'
            'Thalamus_R'
            'Putamen_L'
            'Putamen_R'
            'Postcentral_L'
            'Postcentral_R'
            'Parietal_Inf_L'
            'Parietal_Inf_R'
%             'Cerebellum_Crus1_L'
%            'Cerebellum_Crus1_R'
%            'Vermis_1_2'
%            'Cerebellum_Crus2_L'
%            'Cerebellum_Crus2_R'
%            'Vermis_3'
%            'Cerebellum_3_L'
%            'Cerebellum_3_R'
%            'Vermis_4_5'
%            'Cerebellum_4_5_L'
%            'Cerebellum_4_5_R'
%            'Vermis_6'
%            'Cerebellum_6_L'
%            'Cerebellum_6_R'
%            'Vermis_7'
%            'Cerebellum_7b_L'
%            'Cerebellum_7b_R'
%            'Vermis_8'
%            'Cerebellum_8_L'
%            'Cerebellum_8_R'
%            'Vermis_9'
%            'Cerebellum_9_L'
%            'Cerebellum_9_R'
%            'Vermis_10'
%            'Cerebellum_10_L'
%            'Cerebellum_10_R'
            };
cfg.location = [-8 -10 8]; %% sub-thalamic view
% cfg.location = [0 -34 38]; % cortical view
cfg.crosshair = 'no';

ft_sourceplot(cfg, plot_interp);

fig = gcf;
set(fig, 'units', 'normalized', 'outerposition', [0 0 0.5 1]);
filename = 'theta_stimulation_sub_cortical_look.png';
% print(fullfile(figure_path, filename), '-dpng')