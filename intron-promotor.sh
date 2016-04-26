file=$1
genome=$2
span=$3

echo "Usage:sh ~/sunyd/identify/script/intron-promotor.sh cuff228-CPC_left ~/sunyd/genome/tomato/ITAG2.4_genomic.fa"
#########make the intron, promotor region for lncRNAs gff file###########
subtractBed -a ${file}.gff3 -b ${file}.gtf |sed 's/transcript/intron/g' > ${file}.intron.gff3 
perl /psc/bioinformatics/sunyd/identify/script/promotor.pl ${file}.gff3 ${genome}-length.txt $span > ${file}.promotor.gff3
cat ${file}.gff3 ${file}.intron.gff3 ${file}.promotor.gff3 |sortBed |sed '/\ttranscript\t/d'> ${file}.chip.gff3

