export SUBJECTS_DIR=/scratch5/MINDLAB2019_MEG-CerebellarClock/tactile_jitter/freesurfer

filecontent=( `cat 'paths.txt' `)

for path in "${filecontent[@]}"

do

    mne_make_scalp_surfaces --subject ${path:0:4}

done
