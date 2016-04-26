dir=$1
name=$2    ###the cuffcompare merged combined gtf file

#####make the gtf file for cuffnorm,delete those ambigious assemblies#######
sed '/class_code \"\.\"/d' $dir/${name}.combined.gtf |sed '/class_code \"e\"/d' |sed '/class_code \"p\"/d' |sed '/class_code \"r\"/d' |sed '/class_code \"s\"/d' > ${name}.combined.gtf-cuffnorm       ####delete those ambiguous classifications, get the gtf file for cuffnorm
awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|nearest_ref|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' ${name}.combined.gtf-cuffnorm |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/${name}.combined.gtf-correspond-cuffnorm        #######make the corresponding gene id and name and nearest gene for all del gtf

########make the cufflinks name correspond to genename and transname############
sed -n '/oId "LOC/p' $dir/${name}.combined.gtf |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|oId|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' > $dir/${name}.combined.gtf-correspond-mrnaname     #########only extract oId= Csa and make corresponding geneid, genename,trans id,mrna name classcode in cuffcompare combined gtf, lack 3176 mrna compared to cucumber.gff3 file
sed -n '/=/p' $dir/${name}.combined.gtf |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|oId|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' > $dir/${name}.combined.gtf-correspond-genename   ##########only extract classcode= and make the corresponding genename and id, lack 16 genes
cut -f 2 $dir/${name}.combined.gtf-correspond-genename |sort -u |wc           ###count the gene num
cut -f 2 $dir/${name}.combined.gtf-correspond-mrnaname |sort -u |wc           ###count the mrna num
grep -Fwf $dir/split_cuff/cuff-CPC_left.name $dir/${name}.combined.gtf-correspond-cuffnorm > $dir/${name}-CPC_left.lncname  ###get the lncRNA name corresponding

##############extract all transcripts exon num#######
awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/${name}.combined.gtf-cuffnorm |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' |sed 's/\t\t/\t/g' > $dir/${name}.combined.gtfexontmp1
awk -F"\t" '{OFS="\t"; print $5-$4,$4,$5,$7,$1}' $dir/${name}.combined.gtf-cuffnorm > $dir/${name}.combined.gtfexontmp2
paste $dir/${name}.combined.gtfexontmp1 $dir/${name}.combined.gtfexontmp2 > $dir/${name}.combined.gtf-exon
rm $dir/${name}.combined.gtfexontmp*
cut -f 2,8 $dir/${name}.combined.gtf-exon |sort -u > $dir/${name}.combined.gtf-exon-chr
##############calculate transcript length info###############
Rscript ~/sunyd/identify/script/seqcomp.R $dir/${name}.combinedonelinename.fa $dir/lnclength.txt  ###do not forget to change the last line to the first,calculate the seq length and GC, AU info

##############extract the = kind of transcripts exon num#######
#sed -n '/oId "Csa/p' $dir/${name}.combined.gtf |awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code|oId|gene_name/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/${name}.combined.gtfexontmp1
#sed -n '/oId "Csa/p' $dir/${name}.combined.gtf |awk -F"\t" '{OFS="\t"; print $5-$4,$4,$5,$7}' > $dir/${name}.combined.gtfexontmp2
#paste $dir/${name}.combined.gtfexontmp1 $dir/${name}.combined.gtfexontmp2 > $dir/${name}.combined.gtf-gene-exon
#rm $dir/${name}.combined.gtfexontmp*

##############extract the lncRNA exon num#######
#awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/${name}-CPC_left.gtf |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$6}' |sed 's/\t\t/\t/g' > $dir/${name}-CPC_left.gtf-exontmp1
#awk -F"\t" '{OFS="\t"; print $5-$4,$4,$5,$7}' $dir/${name}-CPC_left.gtf > $dir/${name}-CPC_left.gtf-exontmp2
#paste $dir/${name}-CPC_left.gtf-exontmp1 $dir/${name}-CPC_left.gtf-exontmp2 > $dir/${name}-CPC_left.gtf-exon
#rm $dir/${name}-CPC_left.gtf-exontmp*


