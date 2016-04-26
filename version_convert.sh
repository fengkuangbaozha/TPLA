file=$1
genome1=$2
genome2=$3
script=/psc/home/sunyidan/scripts

gtf_to_fasta $file $genome1 ${file}.fa
blat $genome2 ${file}.fa ${file}.psl  ####
head -n 5 ${file}.psl > tmp1
sed '1,5d' ${file}.psl |sort -k2 -g > tmp2   ########only keep 0 mismatch and change the header
cat tmp1 tmp2 > ${file}.psl0

rm tmp*


