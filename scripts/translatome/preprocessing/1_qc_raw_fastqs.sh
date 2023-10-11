#!/bin/bash
#$ -N translatome_FastqQC
#$ -cwd
#$ -pe smp 6

module load fastqc

concatenated_raw_fastqs="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/concatenated_raw_fastqs/"
output_dir="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/fastqc"
adapter_sequences="/projectnb/mccall/sbandyadka/translatome_reanalysis/scripts/adapter_sequences.txt"

fastqfolder=$1
#echo $fastqfolder
outfiledir="$(echo $fastqfolder | cut -d'/' -f 7)"
echo $output_dir/$outfiledir
mkdir -p $output_dir/$outfiledir

for i in `find $fastqfolder -type f`
do 
    fastqc $i -o $output_dir/$outfiledir
done

