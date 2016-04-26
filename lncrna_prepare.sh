gtflist=$1
gff=$2
genome=$3
dir=$4
name=$5
#num=`echo "$name - 1" |bc`
#/psc/home/sunyidan/tool/cuffcompare -r $gff -o $dir/cuff${name} -i $gtflist  ##run cuffcompare to combine assembly transcripts
awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|nearest_ref|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}.combined.gtf |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/cuff${name}.combined.gtf-correspond        #######make the corresponding gene id and name and nearest gene for all del gtf
awk '{  for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}.combined.gtf |sort -u |cut -d " " -f 4 |sort |uniq -c > $dir/cuff${name}.combined.gtf-classcodecount      ##count the number of each classcode
awk -F" " '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$12,$11,$10,$13,$14,$17,$18,$19,$20,$21,$22,$23,$24}' $dir/cuff${name}.combined.gtf |sed 's/\t"/ "/g' |sed 's/;\t/; /g' > $dir/cuff${name}.combined.gtf.change
gtf2bed --do-not-sort < $dir/cuff${name}.combined.gtf.change > $dir/cuff${name}.combined.bed
fastaFromBed -tab -s -name -fi $genome -bed $dir/cuff${name}.combined.bed -fo $dir/cuff${name}.combined.fa
sort $dir/cuff${name}.combined.fa |awk -F'\t' 'NF>1{a[$1] = a[$1]""$2}END{for(i in a){print i"\t"a[i]}}' |sed 's/\t,/\t/g' |sed 's/^/>/g' |sed 's/\t/\n/g'> $dir/cuff${name}.combined.name.fa
#gtf_to_fasta $dir/cuff${name}.combined.gtf $genome $dir/cuff${name}.combined.fa
#sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $dir/cuff${name}.combined.fa |cut -d " " -f 1,2 |sed -r 's/>[0-9]+ />/g' > $dir/cuff${name}.combined.name.fa  ##change the fasta name and make it into one line

#####make the gtf file for cuffnorm,delete those ambigious assemblies#######
sed -n '/class_code \"\=\"/p' $dir/cuff${name}.combined.gtf > $dir/tmp1
sed -n '/class_code \"\.\"/p' $dir/cuff${name}.combined.gtf > $dir/tmp2.tmp
bedtools intersect -a $dir/tmp2.tmp -b $gff -wa |cut -d " " -f4 |sed 's/"//g' |sed 's/;//g' |sort -u > $dir/cuff${name}.dot.del
grep -Fvf $dir/cuff${name}.dot.del $dir/tmp2.tmp |sed 's/class_code \"\.\"/class_code \"u\"/g'> $dir/tmp2
wc $dir/cuff${name}.dot.del
rm $dir/tmp2.tmp*
sed -n '/class_code \"i\"/p' $dir/cuff${name}.combined.gtf > $dir/tmp3
sed 's/\"e\"/\"o\"/g' $dir/cuff${name}.combined.gtf |sed -n '/class_code \"o\"/p' > $dir/tmp4
sed -n '/class_code \"u\"/p' $dir/cuff${name}.combined.gtf > $dir/tmp5
sed -n '/class_code \"x\"/p' $dir/cuff${name}.combined.gtf > $dir/tmp6
cat $dir/tmp* > $dir/cuff${name}.combined.cuffnorm.gtf       ####delete those ambiguous classifications, get the gtf file for cuffnorm
rm $dir/tmp*
awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|nearest_ref|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}.combined.cuffnorm.gtf |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/cuff${name}.combined.cuffnorm.gtf-correspond       #######make the corresponding gene id and name and nearest gene for all del gtf
cut -f2 $dir/cuff${name}.combined.cuffnorm.gtf-correspond |sort -u > $dir/cuff${name}.combined.cuffnorm.gtf-correspond-name
perl /psc/home/sunyidan/scripts/extractFromFasta.pl $dir/cuff${name}.combined.name.fa list $dir/cuff${name}.combined.cuffnorm.gtf-correspond-name > $dir/cuff${name}.combined.cuffnorm.name.fa
#gtf_to_fasta $dir/cuff${name}.combined.cuffnorm.gtf $genome $dir/cuff${name}.combined.cuffnorm.fa
#sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $dir/cuff${name}.combined.cuffnorm.fa |cut -d " " -f 1,2 |sed -r 's/>[0-9]+ />/g' > $dir/cuff${name}.combined.cuffnorm.name.fa  ##change the fasta name and make it into one line

########prepare the high quality RNA gtf for identification##########
/psc/program/install/R-3.0.2/bin/Rscript /psc/home/sunyidan/sunyd/identify/script/fpkm_filter.R $dir/cuff${name}.tracking $dir/cuff${name}.combined.fpkm.name  ####delete the low quality transcripts
grep -Fwf $dir/cuff${name}.combined.fpkm.name $dir/cuff${name}.combined.cuffnorm.gtf |sort -u > $dir/cuff${name}.combinedhq.gtf
awk '{  for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}.combinedhq.gtf |sort -u |cut -d " " -f 4 |sort |uniq -c > $dir/cuff${name}.combinedhq.gtf-classcodecount      ##count the number of each classcode
cut -d " " -f4 $dir/cuff${name}.combinedhq.gtf |sed 's/"//g' |sed 's/;//g' > $dir/cuff${name}.combinedhq.gtf-coresspond-name
perl /psc/home/sunyidan/scripts/extractFromFasta.pl $dir/cuff${name}.combined.name.fa list $dir/cuff${name}.combinedhq.gtf-coresspond-name > $dir/cuff${name}.combinedhq.name.fa
sed -n '/>/p' $dir/cuff${name}.combinedhq.name.fa |wc
awk '{  for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}'  $dir/cuff${name}.combinedhq.gtf |sort -u |wc   #######count transcript num of gtf file and check if it's the same as fasta file
mkdir $dir/split_cuff
#split -d -a 2 -l 20000 $dir/cuff${name}.combinedhq.name.fa $dir/split_cuff/cuff  ###split the file into small piecies
#cd $dir/split_cuff
#rename '0' '' cuff0*

########make the cufflinks name correspond to genename and transname############
sed '/oId "CUFF/d' $dir/cuff${name}.combined.gtf |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|oId|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' > $dir/cuff${name}.combined.gtf-correspond-mrnaname     #########only extract oId= Csa and make corresponding geneid, genename,trans id,mrna name classcode in cuffcompare combined gtf, lack 3176 mrna compared to cucumber.gff3 file
sed -n '/=/p' $dir/cuff${name}.combined.gtf |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|oId|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' > $dir/cuff${name}.combined.gtf-correspond-genename   ##########only extract classcode= and make the corresponding genename and id, lack 16 genes
cut -f2 $dir/cuff${name}.combined.gtf-correspond-genename |sort > $dir/cuff${name}.combined.gtf-correspond-genename.tmp
grep -Fwf $dir/cuff${name}.combined.gtf-correspond-genename.tmp $dir/cuff${name}.combined.cuffnorm.gtf > $dir/cuff${name}.combined.gtf-gene
cut -f 2 $dir/cuff${name}.combined.gtf-correspond-genename |sort -u |wc           ###count the gene num
cut -f 2 $dir/cuff${name}.combined.gtf-correspond-mrnaname |sort -u |wc           ###count the mrna num
rm $dir/cuff${name}.combined.gtf-correspond-genename.tmp

##############extract all transcripts exon num#######
awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}.combined.cuffnorm.gtf |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' |sed 's/\t\t/\t/g' > $dir/cuff${name}.combined.gtfexontmp1
awk -F"\t" '{OFS="\t"; print $5-$4+1,$4,$5,$7,$1}' $dir/cuff${name}.combined.cuffnorm.gtf > $dir/cuff${name}.combined.gtfexontmp2
paste $dir/cuff${name}.combined.gtfexontmp1 $dir/cuff${name}.combined.gtfexontmp2 |sed 's/\t\t/\t/g' > $dir/cuff${name}.combined.gtf-exon
rm $dir/cuff${name}.combined.gtfexontmp*
cut -f 2,8 $dir/cuff${name}.combined.gtf-exon |sort -u > $dir/cuff${name}.combined.gtf-exon-chr
##############calculate transcript length info###############
/psc/program/install/R-3.0.2/bin/Rscript /psc/bioinformatics/sunyd/identify/script/seqcomp.R $dir/cuff${name}.combined.cuffnorm.name.fa $dir/cuff${name}.sequence_info.txt  ###do not forget to change the last line to the first,calculate the seq length and GC, AU info
/psc/program/install/R-3.0.2/bin/Rscript /psc/bioinformatics/sunyd/identify/script/lengthdistribution.R $dir/cuff${name}.combined.fa $dir/cuff${name}.combined   ###calculate length

##############extract the = kind of transcripts exon num#######
#sed -n '/oId "Csa/p' $dir/cuff${name}.combined.gtf |awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code|oId|gene_name/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/cuff${name}.combined.gtfexontmp1
#sed -n '/oId "Csa/p' $dir/cuff${name}.combined.gtf |awk -F"\t" '{OFS="\t"; print $5-$4,$4,$5,$7}' > $dir/cuff${name}.combined.gtfexontmp2
#paste $dir/cuff${name}.combined.gtfexontmp1 $dir/cuff${name}.combined.gtfexontmp2 > $dir/cuff${name}.combined.gtf-gene-exon
#rm $dir/cuff${name}.combined.gtfexontmp*

##############extract the lncRNA exon num#######
#awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/cuff${name}-CPC_left.gtf |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$6}' |sed 's/\t\t/\t/g' > $dir/cuff${name}-CPC_left.gtf-exontmp1
#awk -F"\t" '{OFS="\t"; print $5-$4,$4,$5,$7}' $dir/cuff${name}-CPC_left.gtf > $dir/cuff${name}-CPC_left.gtf-exontmp2
#paste $dir/cuff${name}-CPC_left.gtf-exontmp1 $dir/cuff${name}-CPC_left.gtf-exontmp2 > $dir/cuff${name}-CPC_left.gtf-exon
#rm $dir/cuff${name}-CPC_left.gtf-exontmp*


