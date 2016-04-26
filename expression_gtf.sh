sample=$1
dir=$2
gtf=$3


cat <(grep -Fvf <(grep -Fwf $dir/${sample}-CPC_left.lncname.iuox.del.name $dir/${sample}-CPC_left.gtf |sed -n '/class_code \"=\"/p' |cut -d " " -f12 |sed 's/"//g' |sed 's/;//g' |sort -u |sed 's/^/transcript_id "/g') ${gtf}.gtf) <(grep -Fwf $dir/${sample}-CPC_left.lncname.iuox.del.name $dir/${sample}-CPC_left.gtf) |sed '/\tCDS\t/d' |sortBed > $dir/${sample}.expression.gene.gtf
awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/${sample}.expression.gene.gtf |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' |sed 's/\t\t/\t/g' > $dir/${sample}.expression.gene.gtf.tmp1
awk -F"\t" '{OFS="\t"; print $5-$4+1,$4,$5,$7,$1}' $dir/${sample}.expression.gene.gtf > $dir/${sample}.expression.gene.gtf.tmp2
paste $dir/${sample}.expression.gene.gtf.tmp1 $dir/${sample}.expression.gene.gtf.tmp2 |sed 's/\t\t/\t/g' |awk -F"\t" '{OFS="\t"; print$1,$2,"code",$3,$4,$5,$6,$7}' > $dir/${sample}.expression.gene.gtf-exon
rm $dir/${sample}.expression.gene.gtf.tmp*
perl ~/scripts/extractFromFasta.pl $dir/${sample}-CPC_left.fa list $dir/${sample}-CPC_left.lncname.iuox.del.name > $dir/${sample}-CPC_left.lncname.iuox.del.fa
cat $dir/${sample}.sequence_info.txt ${gtf}.cdna.sequence_info.txt > $dir/${sample}.sequence_info.all

