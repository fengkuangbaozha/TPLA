genome=$1
sample=$2
location=$3
name=$4
species=$5
adapter=$6
mkdir ${location}/${name}_bsmap
/psc/program/bin/bsmap -a ${sample}1_001.fastq.gz -b ${sample}2_001.fastq.gz -d ${genome}.fa -o ${location}/${name}_bsmap/accepted_hits.bam -v 4 -p 30 -S 1 2> ${location}/${name}_bsmap/align.summary
bam2fq.py -i ${location}/${name}_bsmap/accepted_hits.bam -o ${location}/${name}_bsmap/${name}_${species}
sed -n '/@/p' ${location}/${name}_bsmap/${name}_${species}.R1.fastq > ${location}/${name}_bsmap/r1.name
sed -n '/@/p' ${location}/${name}_bsmap/${name}_${species}.R2.fastq > ${location}/${name}_bsmap/r2.name
cat ${location}/${name}_bsmap/r1.name ${location}/${name}_bsmap/r2.name |sort -u |cut -d/ -f1 |sed 's/@//g' > ${location}/${name}_bsmap/${name}_${species}.name
awk -v n=$adapter '{OFS=""; print $1," 1:N:0:",n}' ${location}/${name}_bsmap/${name}_${species}.name > ${location}/${name}_bsmap/${name}_${species}.name.r1
awk -v n=$adapter '{OFS=""; print $1," 2:N:0:",n}' ${location}/${name}_bsmap/${name}_${species}.name > ${location}/${name}_bsmap/${name}_${species}.name.r2
gunzip ${sample}1_001.fastq.gz
gunzip ${sample}2_001.fastq.gz
perl ~/scripts/fastq-filter_extract_reads.pl -r ${location}/${name}_bsmap/${name}_${species}.name.r1 -f ${sample}1_001.fastq 1>${location}/${name}_bsmap/${name}_${species}.R1.fastq.more 2>${location}/${name}_bsmap/${name}_${species}.R1.fastq.more.not
perl ~/scripts/fastq-filter_extract_reads.pl -r ${location}/${name}_bsmap/${name}_${species}.name.r2 -f ${sample}2_001.fastq 1>${location}/${name}_bsmap/${name}_${species}.R2.fastq.more 2>${location}/${name}_bsmap/${name}_${species}.R2.fastq.more.not
mv ${location}/${name}_bsmap/${name}_${species}.R1.fastq.more ${location}/${name}_bsmap/${name}_${species}.R1.fastq
mv ${location}/${name}_bsmap/${name}_${species}.R2.fastq.more ${location}/${name}_bsmap/${name}_${species}.R2.fastq
gzip ${location}/${name}_bsmap/${name}_${species}.R1.fastq
gzip ${location}/${name}_bsmap/${name}_${species}.R2.fastq
