#!/bin/bash

# exports
filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/mr_paths.txt" `)

for b in "${filecontent[@]}"
do
   subject=${b:0:4}
   SCRIPT=subject_$subject.sh
   OUTPUT=subject_$subject.txt

## Generate a single script
cat <<EOF > ${SCRIPT}
#!/bin/bash
#$ -S /bin/bash
export TERM vt100
# Make sure MNI and CFIN are in the path
PATH=$PATH:/usr/local/mni/bin:/usr/local/cfin/bin:.
# Change to current directory.
cd `pwd`
# Run the programme with the parameters:
./make_scalp_surfaces.sh ${subject} > $OUTPUT
EOF

   # Make the new script exectutable
   chmod u+x ${SCRIPT}

   # Finally, submit it to the cluster
   submit_to_cluster -q all.q -p MINDLAB2019_MEG-CerebellarClock ./${SCRIPT}

done


# filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/mr_paths.txt" `)
#
# for b in "${filecontent[@]}"
# do
#    subject=${b:0:4}
#    SCRIPT=subject_$subject.sh
#    OUTPUT=subject_$subject.txt
#
# ## Generate a single script
# cat <<EOF > ${SCRIPT}
# !/bin/bash
# $ -S /bin/bash
# export TERM vt100
# Make sure MNI and CFIN are in the path
# PATH=$PATH:/usr/local/mni/bin:/usr/local/cfin/bin:.
# Change to current directory.
# cd `pwd`
# Run the programme with the parameters:
# ./import_MRI.sh ${b} > $OUTPUT
# EOF
#
#    # Make the new script exectutable
#    chmod u+x ${SCRIPT}
#
#    # Finally, submit it to the cluster
#    submit_to_cluster -q short.q -p MINDLAB2019_MEG-CerebellarClock ./${SCRIPT}
#
# done
#
# filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/mr_paths.txt" `)
#
# for b in "${filecontent[@]}"
# do
#     subject=${b:0:4}
#     SCRIPT=subject_$subject.sh
#     OUTPUT=subject_$subject.txt
#
# # Generate a single script
# cat <<EOF > ${SCRIPT}
# #!/bin/bash
# #$ -S /bin/bash
# export TERM vt100
# # Make sure MNI and CFIN are in the path
# PATH=$PATH:/usr/local/mni/bin:/usr/local/cfin/bin:.
# # Change to current directory.
# cd `pwd`
# # Run the programme with the parameters:
# ./reconstruct.sh ${b} > $OUTPUT
# EOF
#
#     # Make the new script exectutable
#     chmod u+x ${SCRIPT}
#
#     # Finally, submit it to the cluster
#     submit_to_cluster -q long.q -p MINDLAB2019_MEG-CerebellarClock ./${SCRIPT}
#
# done

## MAXFILTERIJNG NON-SPLITS

filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/raw_meg_paths.txt" `)
raw_path="/projects/MINDLAB2019_MEG-CerebellarClock/raw/"
scratch_path="/projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/MEG/"
for subject_path in "${filecontent[@]}"
do
    cd ${raw_path}${subject_path}/MEG/001.tactile_jitter_raw/files
    pwd
    for filename in `ls *.fif* -p | grep -v / `
    do
        # echo $filename
        file_ending=${filename:18} ## 18
        output_name=${filename:0:18}_sss${file_ending}
        echo input: $filename
        echo output: $output_name
        full_output=${scratch_path}${subject_path}/${output_name}
#
    # Finally, submit it to the cluster
         submit_to_cluster -q maxfilter.q -p MINDLAB2019_MEG-CerebellarClock \
         -n 4 -p MINDLAB2019_MEG-CerebellarClock "/neuro/bin/util/maxfilter -f $filename -o $full_output"
     done
 done


## MAXFILTERING SPLITS

# filecontent=( `cat "/projects/MINDLAB2019_MEG-CerebellarClock/scripts/tactile_jitter/bash/meg_split_paths.txt" `)
# raw_path="/projects/MINDLAB2019_MEG-CerebellarClock/raw/"
# scratch_path="/projects/MINDLAB2019_MEG-CerebellarClock/scratch/tactile_jitter/MEG/"
# for subject_path in "${filecontent[@]}"
# do
#     cd ${raw_path}${subject_path}
#     pwd
#     for filename in `ls *.fif* -p | grep -v / `
#     do
#         # echo $filename
#         file_ending=${filename:20} ## 20
#         output_name=${filename:0:20}_sss${file_ending}
#         echo input: $filename
#         echo output: $output_name
#         full_output=${scratch_path}${subject_path:0:20}/${output_name}
#
#     # Finally, submit it to the cluster
#         submit_to_cluster -q maxfilter.q -p MINDLAB2019_MEG-CerebellarClock \
#         -n 4 -p MINDLAB2019_MEG-CerebellarClock "/neuro/bin/util/maxfilter -f $filename -o $full_output"
#     done
# done
