#!/bin/bash
#$ -S /bin/bash
export TERM vt100
# Make sure MNI and CFIN are in the path
PATH=/users/lau/MNE-2.7.0-3106-Linux-x86_64/bin:/usr/local/common/GridEngine/bin/lx-amd64:/usr/local/mni/bin:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/fsl/bin:/usr/local/freesurfer/mni/bin:/users/lau/miniconda3/bin:/users/lau/miniconda3/condabin:/users/lau/.local/bin:/usr/local/common/anaconda/bin:/usr/local/common/GridEngine/bin/lx-amd64:/usr/local/mni/bin:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/fsl/bin:/usr/local/freesurfer/mni/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/cfin/bin:/usr/local/fsl/bin:/usr/local/bin:/usr/local/mrtrix3/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/cfin/bin:/usr/local/fsl/bin:/usr/local/bin:/usr/local/mrtrix3/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/mni/bin:/usr/local/cfin/bin:.
# Change to current directory.
cd /projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash
# Run the programme with the parameters:
./make_scalp_surfaces.sh 0005 > subject_0005.txt
