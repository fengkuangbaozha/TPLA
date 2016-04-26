location=$1
name=$2
species=$3
adapter=$4
sample=$5
sed -n '/@HWI-/p' ${location}/${name}_${species}.R1.fastq > ${location}/r1.name
sed -n '/@HWI-/p' ${location}/${name}_${species}.R2.fastq > ${location}/r2.name
cat ${location}/r1.name ${location}/r2.name |cut -d " " -f1 |cut -d/ -f1 |sed 's/@//g' |sort -u > ${location}/${name}_${species}.name
awk -v n=$adapter '{OFS=""; print $1," 1:N:0:",n}' ${location}/${name}_${species}.name > ${location}/${name}_${species}.name.r1
awk -v n=$adapter '{OFS=""; print $1," 4:N:0:",n}' ${location}/${name}_${species}.name > ${location}/${name}_${species}.name.r2
#gunzip ${sample}1_001.fastq.gz
#gunzip ${sample}4_001.fastq.gz
perl ~/scripts/fastq-filter_extract_reads.pl -r ${location}/${name}_${species}.name.r1 -f ${sample}1_001.fastq 1>${location}/${name}_${species}.R1.fastq.more 2>${location}/${name}_${species}.R1.fastq.more.not
perl ~/scripts/fastq-filter_extract_reads.pl -r ${location}/${name}_${species}.name.r2 -f ${sample}4_001.fastq 1>${location}/${name}_${species}.R4.fastq.more 2>${location}/${name}_${species}.R4.fastq.more.not
#mv ${location}/${name}_${species}.R1.fastq.more ${location}/${name}_${species}.R1.fastq
#mv ${location}/${name}_${species}.R4.fastq.more ${location}/${name}_${species}.R4.fastq
#gzip ${location}/${name}_${species}.R1.fastq
#gzip ${location}/${name}_${species}.R4.fastq
