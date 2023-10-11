#!/bin/bash
#$ -N salmonQuant
#$ -cwd
#$ -pe smp 15

salmon_path="/projectnb/mccall/sbandyadka/translatome_reanalysis/software/salmon-latest_linux_x86_64/bin"
#reference_fasta="/projectnb/mccall/sbandyadka/translatome_reanalysis/reference/BDGP6.32.104/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz"
reference_fasta="/projectnb/mccall/sbandyadka/translatome_reanalysis/reference/BDGP6.32.104/Drosophila_melanogaster.BDGP6.32.cdna.all.fa.gz"
#salmon_index="/projectnb/mccall/sbandyadka/translatome_reanalysis/reference/BDGP6.32.104/BDGP6.32.104.toplevel.index"
salmon_index="/projectnb/mccall/sbandyadka/translatome_reanalysis/reference/BDGP6.32.104/BDGP6.32.104.cdna.index"

## create salmon index

#$salmon_path/salmon index -t $reference_fasta -i $salmon_index

## quantify reads

quant_outputfiles="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/salmon_BDGP6_32_104"
samples="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/concatenated_raw_fastqs"

mkdir -p $quant_outputfiles

for i in `ls $samples/*.fastq`; do 
filename=$(basename -- "$i")
echo $quant_outputfiles/${filename}_quant
qsub -P mccall <<<  "$salmon_path/salmon quant -i $salmon_index -l A --gcBias -r $i --validateMappings -o $quant_outputfiles/${filename}_quant"
done