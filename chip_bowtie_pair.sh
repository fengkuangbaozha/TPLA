#$ -cwd

genome=$1
dir=$2
prefix=$3
file=${prefix}${SGE_TASK_ID}
mkdir $dir/$file
fastq-dump --split-3 $dir/${file}.sra -O $dir/${file}
DynamicTrim $dir/${file}/${file}_1.fastq -h 17 -d $dir/${file} 
DynamicTrim $dir/${file}/${file}_2.fastq -h 17 -d $dir/${file} 
LengthSort $dir/$file/${file}_1.fastq.trimmed $dir/$file/${file}_2.fastq.trimmed -l 25 -d $dir/$file
bowtie -p 30 --chunkmbs 200 $genome -1 $dir/$file/${file}_1.fastq.trimmed.paired1 -2 $dir/$file/${file}_1.fastq.trimmed.paired2 -S $dir/$file/${file}.sam
samtools view -bS $dir/$file/${file}.sam > $dir/$file/${file}.bam
rm $dir/${file}/${file}_*.fastq
# end of job script

