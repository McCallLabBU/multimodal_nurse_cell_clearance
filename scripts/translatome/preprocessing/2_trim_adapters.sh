#!/bin/bash
#$ -N translatome_FastqQC
#$ -cwd
#$ -pe smp 6


module load cutadapt

concatenated_raw_fastqs="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/concatenated_raw_fastqs/"
trimmed_fastqs_3keepSampleIndex="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/trimmed_fastqs_3keepSampleIndex"
trimmed_fastqs_3removeSampleIndex="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/trimmed_fastqs_3removeSampleIndex"
trimmed_fastqs_5keepSampleIndex="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/trimmed_fastqs_5keepSampleIndex"
trimmed_fastqs_5removeSampleIndex="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/trimmed_fastqs_5removeSampleIndex"
adapter_Seq="GATCGGAAGAGCACACGTCTGAACTCCAGTC"


mkdir -p $trimmed_fastqs_3keepSampleIndex
mkdir -p $trimmed_fastqs_3removeSampleIndex
mkdir -p $trimmed_fastqs_5keepSampleIndex
mkdir -p $trimmed_fastqs_5removeSampleIndex

#for i in `find $concatenated_raw_fastqs -type f`
#do 
#outfile="$(echo $i | cut -d'/' -f 8)"
#cutadapt -O 4 -a $adapter_Seq $i -o $trimmed_fastqs_3keepSampleIndex/$outfile ##treat as 3prime adapter
#cutadapt -O 4 -g $adapter_Seq $i -o $trimmed_fastqs_5keepSampleIndex/$outfile ##treat as 5prime adapter
#done

cutadapt -O 4 $concatenated_raw_fastqs/GC218.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGC -o $trimmed_fastqs_3removeSampleIndex/GC218.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GC220.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAA  -o $trimmed_fastqs_3removeSampleIndex/GC220.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GC226.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTG  -o $trimmed_fastqs_3removeSampleIndex/GC226.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS209.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAAT  -o $trimmed_fastqs_3removeSampleIndex/GS209.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS211.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATC  -o $trimmed_fastqs_3removeSampleIndex/GS211.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS217-225.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTA  -o $trimmed_fastqs_3removeSampleIndex/GS217-225.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P203.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAA  -o $trimmed_fastqs_3removeSampleIndex/P203.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P205.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCC  -o $trimmed_fastqs_3removeSampleIndex/P205.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P207.fastq -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCA  -o $trimmed_fastqs_3removeSampleIndex/P207.fastq


cutadapt -O 4 $concatenated_raw_fastqs/GC218.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGC -o $trimmed_fastqs_5removeSampleIndex/GC218.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GC220.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAA  -o $trimmed_fastqs_5removeSampleIndex/GC220.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GC226.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTG  -o $trimmed_fastqs_5removeSampleIndex/GC226.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS209.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAAT  -o $trimmed_fastqs_5removeSampleIndex/GS209.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS211.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATC  -o $trimmed_fastqs_5removeSampleIndex/GS211.fastq
cutadapt -O 4 $concatenated_raw_fastqs/GS217-225.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTA  -o $trimmed_fastqs_5removeSampleIndex/GS217-225.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P203.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAA  -o $trimmed_fastqs_5removeSampleIndex/P203.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P205.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCC  -o $trimmed_fastqs_5removeSampleIndex/P205.fastq
cutadapt -O 4 $concatenated_raw_fastqs/P207.fastq -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCA  -o $trimmed_fastqs_5removeSampleIndex/P207.fastq