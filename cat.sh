sample=$1   ####cat file dir
catin=$2    ####the initial name of the cat files
sym=$3      ###whether cuff or ${sym} file
#catout=$2
#number2=`echo "$catin+11" |bc`
#number3=`echo "$catin+22" |bc`
#number4=`echo "$catin+33" |bc`

#cat SRR1138$catin/SRR1138${catin}_1.fastq sra/SRR1138${number2}_1.fastq sra/SRR1138${number3}_1.fastq sra/SRR1138${number4}_1.fastq > CSRR1138${catout}.fastq

for i in $sample/${catin}*;
do
	j=`echo $sample |rev |cut -d/ -f 1 |rev`;
	z=`echo $i |rev |cut -d/ -f 1 |rev`;
	cat $i/${z}_${sym}.fa >> $sample/${j}_${sym}.fa
done;
