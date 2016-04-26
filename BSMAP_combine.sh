#$cwd

dir=$1
file=$2
genome=$3

#mkdir $dir/$file
##fastq-dump --split-3 $dir/${file}.sra -O $dir/${file}
#DynamicTrim $dir/${file}/${file}.fastq -h 17 -d $dir/${file} 
##LengthSort $dir/$file/${file}.fastq.trimmed -l 25 -d $dir/$file
#cutadapt -a AGATCGGAAGAG -f fastq $dir/$file/${file}.fastq.trimmed -o $dir/$file/${file}.fastq.trimmed.cut
#LengthSort $dir/$file/${file}.fastq.trimmed.cut -l 25 -d $dir/$file
bsmap -a $dir/$file/${file}.fastq.trimmed.cut.single -d $genome -o $dir/$file/${file}.bsmap.bsp -v 2 -p 30 -S 1 2> $dir/$file/${file}.mapping.sta
methratio.py -r -d $genome -z $dir/$file/${file}.bsmap.bsp -o $dir/$file/${file}_bsmap_report.txt > $dir/$file/${file}.count.sta
Rscript /psc/bioinformatics/sunyd/identify/script/methylation_gtf.R $dir/$file/${file}_bsmap_report.txt $dir/$file/${file}.gtf


# end of job script


