user=$(uname -n)

if [ $user == 'hyades02' ]
then
    filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/all_meg_paths.txt" `)
elif [ $user == 'lau' ]
then
    filecontent=( `cat "/home/lau/projects/cerebellar_clock/scripts/tactile_jitter/bash/all_meg_paths.txt" `)
fi

fig_string=''
MEG_path=/projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/MEG/


for path in "${filecontent[@]}"
do
    subject=${path:0:4}
    if [ $user == 'hyades02' ]
    then
        cd ${MEG_path}${path}/hilbert/connectivity/stcs
        # cd /projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/freesurfer/$subject/bem/
    elif [ $user == 'lau' ]
    then
        echo lala
        # cd /home/lau/projects/cerebellar_clock/scratch/tactile_jitter/figures/$path/hilbert/stcs/common_filter
    fi
    pwd
    ls
    # fig_string="$fig_string `pwd`/normal_ISI_full_baseline_mag_reg_0.0_weight_unit-noise-gain_o15_7.5_mm_common_filter_o0_o5_o15_hp_14_Hz_lp_30_Hz_-750_ms_to_750_ms_Right_Cerebellum_6_beamfomer_time_courses.png"

    # string=normal
    # for filename in *
    # do
    #
    #         if [ ! -d $filename ]  && [[ $filename != *"baseline"* ]] && [[ $filename == *"-ave.fif"* ]]
    #         then
    #             # echo $filename
    #             rm -v $filename
    #             # new_name=${filename/no_baselinemag/no_baseline_mag}
    #             # # echo $new_name
    #             # echo mv -v $filename $new_name
    #
    #             # new_name=${filename/ms_tactile/-400_ms_to_1100_ms_tactile}
    #             # echo mv -v $filename $new_name
    #         fi
    # done


    # cd /projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/freesurfer/${subject}/bem
    # pwd
    # ln -sfv ./watershed/${subject}_brain_surface brain_surface.surf
    # ln -sfv ./watershed/${subject}_inner_skull_surface inner_skull.surf
    # ln -sfv ./watershed/${subject}_outer_skin_surface outer_skin.surf
    # ln -sfv ./watershed/${subject}_outer_skull_surface outer_skull.surf
done
# xviewer $fig_string &
#
#
# cd /projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/MEG/grand_averages/hilbert
# pwd
#
# for filename in *
# do
#
#         if [ ! -d $filename ] && [[ $filename == *"z_transform_contrast-ave"* ]]
#             then
#              # echo $filename
#              new_name=${filename/Hz_proj/Hz_-400_ms_to_1100_ms_proj}
#              # echo $new_name
#              mv -v $filename $new_name
#
#              # new_name=${filename/ms_proj_tactile/range_-400_ms_to_1100_ms_proj_tactile}
#              # mv -v $filename $new_name
#         fi
# done
