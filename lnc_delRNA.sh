file=$1
sample=$2

bedtools intersect -a $file -b /psc/bioinformatics/sunyd/genome/Cucumber/Cucumber_v2i.hkRNA.gff3 |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$6,"HousekeepRNA"}' > ${sample}_hkRNA.name        ############identify structure RNAs and small RNA precursors
bedtools intersect -a $file -b /psc/bioinformatics/sunyd/genome/Cucumber/Cucumber_v2i.miRNA.gff3 |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$6,"SmallRNA"}' > ${sample}_miRNA.name
bedtools intersect -a $file -b /psc/bioinformatics/sunyd/genome/Cucumber/Cucumber_v2i.TE.gff3 |awk '{	for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$6,"TE"}' > ${sample}_TE.name
cat ${sample}_hkRNA.name ${sample}_miRNA.name ${sample}_TE.name > ${sample}_delRNA.name

