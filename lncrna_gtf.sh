dir=$1   ##working dir
name=$2

########get the gtf file for lncRNA##########
grep -Fwf $dir/split_cuff/cuff-CPC_left.name $dir/${name}.combinedhq.gtf > $dir/${name}-CPC_left.gtf   ####generate lncRNA gtf file
cut -d " " -f 4 $dir/${name}-CPC_left.gtf |sort -u |sed 's/;//'> $dir/${name}-CPC_left.gtfname     
diff $dir/${name}-CPC_left.gtfname $dir/split_cuff/cuff-CPC_left.name |grep "<" |sed 's/<\ //g'> $dir/${name}containlncgtf.txt   ####find the contained lncRNA lines, make sure two files are sorted
cut -d " " -f 1,2,3,4 $dir/${name}-CPC_left.gtf > $dir/${name}-CPC_left.gtftmp                
grep -Fwf $dir/${name}containlncgtf.txt $dir/${name}-CPC_left.gtftmp > $dir/${name}-CPC_left.gtftmpdiff      ####make only the transcript ID match
grep -Fvf $dir/${name}-CPC_left.gtftmpdiff $dir/${name}-CPC_left.gtf > $dir/${name}-CPC_left.gtf1     ##delete the contained ones
mv $dir/${name}-CPC_left.gtf1 $dir/${name}-CPC_left.gtf
gffread $dir/${name}-CPC_left.gtf -o $dir/${name}-CPC_left.gff3
sed 's/\ttranscript\t/\tgene\t/g' $dir/${name}-CPC_left.gff3 > $dir/${name}-CPC_leftnamechange.gff3
cut -d " " -f 2,4 $dir/${name}-CPC_left.gtf |sed 's/;//g' |sort -u |wc              ###check if the number is the same
wc $dir/split_cuff/cuff-CPC_left.name
cut -d " " -f 2,4 $dir/${name}-CPC_left.gtf |sed 's/;//g' |sort -u |cut -d " " -f 1 |uniq -c > cuff-CPC_left.gtfisocount   ##count the isoforms
###########extract the gene number for isoform lnRNAs############
#grep -Fwf $dir/split_cuff/cuff-CPC_left.name $dir/${name}.tracking |cut -f 2 |sort -u > $dir/${name}_gene.name     ###extract the gene number
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

