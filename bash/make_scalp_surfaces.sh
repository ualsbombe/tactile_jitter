
make_scalp_surfaces () {
    export SUBJECTS_DIR=/projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/freesurfer
    echo $1
    subject=$1
    mne_make_scalp_surfaces --subject $subject
}

make_scalp_surfaces $1
