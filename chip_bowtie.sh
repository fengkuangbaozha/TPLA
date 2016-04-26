#$ -cwd

genome=$1
dir=$2
prefix=$3
file=${prefix}${SGE_TASK_ID}
mkdir $dir/$file
fastq-dump --split-3 $dir/${file}.sra -O $dir/${file}
DynamicTrim $dir/${file}/${file}.fastq -h 17 -d $dir/${file} 
LengthSort $dir/$file/${file}.fastq.trimmed -l 25 -d $dir/$file
bowtie -p 30 $genome $dir/$file/${file}.fastq.trimmed.single -S $dir/$file/${file}.sam
samtools view -bS $dir/$file/${file}.sam > $dir/$file/${file}.bam
rm $dir/${file}.fastq
# end of job script


