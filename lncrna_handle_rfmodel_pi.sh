dir=$1
name=$2
script=/psc/home/sunyidan/scripts

##########combine the lncrna_result from both identification methods #############
cat $dir/split_cuff/cuff_find*/*-lenfilterleft.fa > $dir/split_cuff/cuff-lenfilterleft.fa
cat $dir/split_cuff/cuff_find*/*-hmmer_left.fa > $dir/split_cuff/cuff-hmmer_left.fa
cat $dir/split_cuff/cuff_find*/*-protdb_blast_left.fa > $dir/split_cuff/cuff-protdb_blast_left.fa
cat $dir/split_cuff/cuff_find*/*-infer_left.fa > $dir/split_cuff/cuff-infer_left.fa
cat $dir/split_cuff/cuff_find*/*-CPC_left.fa > $dir/split_cuff/cuff-CPC_left.fa
sed -n '/>/p' $dir/split_cuff/cuff-CPC_left.fa |sed 's/>//g' > $dir/split_cuff/cuff-CPC_left.name1
cat $dir/split_cuff/cuff_find*/*-rfmodel_seqstruall_score.txt > $dir/split_cuff/cuff-rfmodel_seqstruall_score.txt
sed '/\t0.5\t/d' $dir/split_cuff/cuff-rfmodel_seqstruall_score.txt |cut -f 1 |sort -u > $dir/split_cuff/cuff-rfmodel_seqstruall.name
grep -Fvf $dir/${name}.dot.del $dir/split_cuff/cuff-CPC_left.name1 > $dir/split_cuff/cuff-CPC_left.name

########get the gtf file for lncRNA##########
grep -Fwf $dir/split_cuff/cuff-CPC_left.name $dir/${name}.combined.cuffnorm.gtf > $dir/${name}-CPC_left.gtf   ####generate lncRNA gtf file
awk '{  for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|gene_name|nearest_ref|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/${name}-CPC_left.gtf |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $2,$4,$10,$6,$8}' |sed 's/\t\t/\t/g' > $dir/${name}-CPC_left.lncname
#sed '/\t=\t/d' $dir/${name}-CPC_left.lncname > $dir/${name}-CPC_left.lncname.iuox
#sed '/\t=\t/d' $dir/${name}-CPC_left.lncname |sed '/\tj\t/d' > $dir/${name}-CPC_left.lncname.iuox
#sed '/\ti\t/d' $dir/${name}-CPC_left.lncname |sed '/\tu\t/d' |sed '/\to\t/d' |sed '/\tx\t/d'  > $dir/${name}-CPC_left.lncname.ej
gffread $dir/${name}-CPC_left.gtf -o $dir/${name}-CPC_left.gff3
sed 's/\ttranscript\t/\tgene\t/g' $dir/${name}-CPC_left.gff3 > $dir/${name}-CPC_leftnamechange.gff3
cut -d " " -f 2,4 $dir/${name}-CPC_left.gtf |sed 's/;//g' |sort -u |wc              ###check if the number is the same
wc $dir/split_cuff/cuff-CPC_left.name
cut -d " " -f 2,4 $dir/${name}-CPC_left.gtf |sed 's/;//g' |sort -u |cut -d " " -f 1 |uniq -c > cuff-CPC_left.gtfisocount   ##count the isoforms
awk '{  for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' $dir/${name}-CPC_left.gtf |sort -u |cut -d " " -f 4 |sort |uniq -c > $dir/${name}-CPC_left.gtf-classcodecount      ##count the number of each classcode
Rscript /psc/bioinformatics/sunyd/identify/script/lengthdistribution.R $dir/split_cuff/cuff-CPC_left.fa $dir/${name}-CPC_left   ###calculate length

###########extract the gene number for isoform lnRNAs############
#grep -Fwf $dir/split_cuff/cuff-rfmodel_seqstruall.name $dir/${name}.tracking |cut -f 2 |sort -u > $dir/${name}_gene.name     ###extract the gene number
#grep -Fwf $dir/${name}_gene.name $dir/${name}.loci |sed 's/\[/\t\.\t/g' |sed 's/\]/\t\.\t/g' |sed -e 's/\-[[:digit:]]/(&/g' |sed 's/(-/\t/g'> $dir/${name}_gene.loci  ###extract the gene info in the loci file
#awk -F"\t" '{OFS=""; print $2,"\t","Cufflinks","\t","exon","\t",$6,"\t",$7,"\t",$3,"\t",$4,"\t",$5,"\t","transcript_id"," ",$1 ";"," ","gene_name"," ",$8 ";"}' $dir/${name}_gene.loci |sed 's/ gene_name ;//g'> $dir/${name}_gene.gtf
#gffread $dir/${name}_gene.gtf -o $dir/${name}_gene.gff3
#sed 's/\ttranscript\t/\tgene\t/g' $dir/${name}_gene.gff3 > $dir/${name}_genenamechange.gff3
#sed -n '/\tgene\t/p' $dir/${name}_genenamechange.gff3 |wc
#sed '1,2d' $dir/${name}_genenamechange.gff3 > $dir/${name}_genenamechange.gff3tmp
#########use gffread to combine the isoforms#############
#gffread $dir/${name}-CPC_left.gtf -M -o $dir/${name}_generead.gff3       ######convert the isoform gtf file into loci gff file to do gene expression analysis
#sed 's/\tlocus\t/\tgene\t/g' $dir/${name}_generead.gff3 |sed 's/\ttranscript\t/\tmRNA\t/g' |sed '1d' |sed '1d' > $dir/${name}_genereadnamechange.gff3         ###change the name like cucumber.gff3
#sed -n '/\tgene\t/p' $dir/${name}_genereadnamechange.gff3 |wc
#####cat /psc/bioinformatics/sunyd/genome/Cucumber/Cucumber.gff3 cuff-CPC_leftgenenamechange.gff3 > /psc/bioinformatics/sunyd/genome/Cucumber/CucumberlncRNA.gff3

