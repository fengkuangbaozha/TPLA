gtflist=$1
gff=$2
genome=$3
dir=$4
name=$5
num=`echo "$name - 1" |bc`
###/psc/bioinformatics/sunyd/genome/Cucumber/Cucumber.gff3
/psc/home/sunyidan/tool/cuffcompare -r $gff -C -o $dir/cuff$name -i $gtflist  ##run cuffcompare to combine assembly transcripts
gtf_to_fasta $dir/cuff${name}.combined.gtf $genome $dir/cuff${name}.combined.fa
Rscript ~/sunyd/identify/script/lengthdistribution.R $dir/cuff${name}.combined.fa cuff${name}.combined
sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $dir/cuff${name}.combined.fa |cut -d " " -f 1,2 |sed -r 's/>[0-9]+ />/g' > $dir/cuff${name}.combinedonelinename.fa  ##change the fasta name and make it into one line
awk '{
  for (i = 1; i <= NF; i++) {
  if ($i ~ /gene_id|transcript_id|class_code|oId/){
  printf "%s %s", $i, $(i+1)
  }
  }
  print ""
  }' $dir/cuff${name}.combined.gtf |sed 's/;/\t/g' > $dir/cuff${name}.combined.gtfname    #########only extract gene name and class code, oId
awk '{
  for (i = 1; i <= NF; i++) {
  if ($i ~ /gene_id|transcript_id|class_code/){
  printf "%s %s", $i, $(i+1)
  }
  }
  print ""
  }' $dir/cuff${name}.combined.gtf |sort -u |cut -d " " -f 4 |sort |uniq -c > $dir/cuff${name}.combined.gtfclasscodeuniqcount      ##only extract class code and gene name
Rscript /psc/home/sunyidan/sunyd/identify/script/high_assemble.R $dir/cuff${name}.loci $dir/cuff${name}.tracking $dir/cuff${name}.combined.gtf $dir/cuff${name}.combinedhq.gtf $num  ####delete the low quality transcripts
awk '{
  for (i = 1; i <= NF; i++) {
  if ($i ~ /gene_id|transcript_id|class_code/){
  printf "%s %s", $i, $(i+1)
  }
  }
  print ""
  }' $dir/cuff${name}.combinedhq.gtf |sort -u |cut -d " " -f 4 |sort |uniq -c > $dir/cuff${name}.combinedhq.gtfclasscodeuniqcount  ##only extract class code and gene name
gtf_to_fasta $dir/cuff${name}.combinedhq.gtf $genome $dir/cuff${name}.combinedhq.fa   ###convert gtf into fasta file
sed -n '/>/p' $dir/cuff${name}.combinedhq.fa |wc
awk '{
 for (i = 1; i <= NF; i++) {
 if ($i ~ /gene_id|transcript_id|class_code/){
 printf "%s %s", $i, $(i+1)
 }
 }
 print ""
 }' $dir/cuff${name}.combinedhq.gtf |sort -u |wc   #######count transcript num of gtf file and check if it's the same as fasta file
sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $dir/cuff${name}.combinedhq.fa |cut -d " " -f 1,2 |sed -r 's/>[0-9]+ />/g' > $dir/cuff${name}.combinedhqonelinename.fa  ##change the fasta name and make it into one line
mkdir $dir/split_cuff
split -d -a 2 -l 20000 $dir/cuff${name}.combinedhqonelinename.fa split_cuff/cuff  ###split the file into small piecies
cd $dir/split_cuff
rename '0' '' cuff0*

