
input_folder="/projectnb/mccall/sbandyadka/translatome_reanalysis/data"
output_folder="/projectnb/mccall/sbandyadka/translatome_reanalysis/analysis2/concatenated_raw_fastqs"

mkdir -p $output_folder
for i in `ls $input_folder`;
do 
	cd $input_folder/$i;
	file="$(echo $input_folder/$i/*L001* | cut -d'_' -f 3 | cut -d'/' -f 2)";
	echo $file;
	cat ${input_folder}/${i}/*.fastq > ${output_folder}/${file}.fastq;
	cd $input_folder;

done



