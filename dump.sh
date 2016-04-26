dir=$1   ####dump file dir
sdir=$2    ####the initial name of the dump files dir
sample=$3
for i in $dir/${sdir}*/${sample}*.sra;
do
	j=`echo $i |rev |cut -d/ -f 2 |rev`;
#	z=`echo $i |rev |cut -d/ -f 3 |rev |cut -d "." -f 1`;
	cd $dir/$j;
	fastq-dump --split-files $i
done;
