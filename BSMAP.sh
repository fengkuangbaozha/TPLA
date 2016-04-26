#$ -cwd

dir=$1
file=$2
genome=$3

mkdir $dir/${file}${SGE_TASK_ID}
fastq-dump --split-3 $dir/${file}${SGE_TASK_ID}.sra -O $dir/${file}${SGE_TASK_ID}
DynamicTrim $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq -h 17 -d $dir/${file}${SGE_TASK_ID} 
#LengthSort $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq.trimmed -l 25 -d $dir/${file}${SGE_TASK_ID}
cutadapt -a AGATCGGAAGAG -f fastq $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq.trimmed -o $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq.trimmed.cut
LengthSort $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq.trimmed.cut -l 25 -d $dir/${file}${SGE_TASK_ID}
bsmap -a $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.fastq.trimmed.cut.single -d $genome -o $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.bsmap.bsp -v 2 -p 30 -S 1 2> $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.mapping.sta
methratio.py -r -d $genome -z $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.bsmap.bsp -o $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}_bsmap_report.txt 1> $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.count.sta
sed '1d' $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}root_bsmap_report.txt |awk -F"\t" '{OFS="\t"; print $1,"BSMAP","exon",$2,$2,".",$3,".",$4"; "$7"; "$8}' > $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.gtf 
#Rscript /psc/bioinformatics/sunyd/identify/script/methylation_gtf.R $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}_bsmap_report.txt $dir/${file}${SGE_TASK_ID}/${file}${SGE_TASK_ID}.gtf


# end of job script


