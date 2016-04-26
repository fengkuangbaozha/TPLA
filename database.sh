dir=$1
species=$2
name=$3
go=$4
latin=$5  ###"Solanum_lycopersicum"
genome=$6      ########"ITAG2.4"

mkdir -p $dir/database
###############basic info#############
#sed -n '/\ttranscript\t/p' $dir/${name}-CPC_left.gff3 |cut -d ";" -f 1 |sed 's/ID=//g' |awk -F"\t" '{OFS="\t"; print$9,$1,$4,$5,$7}' > $dir/database/${name}-CPC_left.basic  ###get the transcript start and end info
#cat $dir/split_cuff/cuff_find*/*.cpcsort |sort -u |cut -f1,3 > $dir/database/cpc.score
sort $dir/heatmap/lnc.all |awk -v n=$latin -v g=$genome -F"\t" '{OFS="\t"; print$1,$9,n,g,$10,$2,$3,$4,$5,$6,$7,$8,$11,$12,$13}' |sed 's/\ti\t/\tintron\t/g' |sed 's/\tu\t/\tintergenic\t/g' |sed 's/\to\t/\texon_overlap\t/g' |sed 's/\tx\t/\tantisense\t/g' |sed 's/\t=\t/\texon_overlap\t/g' > $dir/database/${species}.basic.txt 
cat /psc/bioinformatics/sunyd/database/header $dir/database/${species}.basic.txt > $dir/database/${species}.basic.txt1
Rscript /psc/bioinformatics/sunyd/identify/script/strand_change.R $dir/database/${species}.basic.txt1 $dir/database/${species}.basic.txt2
sed 's/X.transcript_id/#transcript_id/g' $dir/database/${species}.basic.txt2 > $dir/database/${species}.basic.txt
rm $dir/database/${species}.basic.txt1
rm $dir/database/${species}.basic.txt2
#perl /psc/home/sunyidan/scripts/extractFromFasta.pl $dir/${name}-CPC_left.lncname.iuox.del.fa list $dir/heatmap/exp.sort >  $dir/database/${species}.lncrna.fa
cut -f1,14 $dir/heatmap/lnc.all > $dir/database/${species}.lncrna.oneline.fa
sed 's/^/>/g' $dir/database/${species}.lncrna.oneline.fa |sed 's/\t/\n/g' |sed '/^$/d' > $dir/database/${species}.lncrna.fa

#######extract gtf file#########
##grep -Fwf $dir/${name}-CPC_left.lncname.iuox.del.name $dir/${name}-CPC_left.gtf > $dir/database/${name}-CPC_left.gtf.iuox.del    ###for the keep gtf file
grep -Fwf $dir/heatmap/exp.sort $dir/${name}-CPC_left.gff3 > $dir/database/${species}.lncrna.gff3

####mirna target and mimic###########
cat $dir/mirna/psrobot.mirna*mimic.name |sed 's/>//g' |cut -f1,2 |sort -u > $dir/mirna/psrobot.mirna.mimic.name
grep -Fwf $dir/heatmap/exp.sort <(cut -f 1,2 $dir/mirna/psrobot.mirna.mimic.name) |sed 's/>//g' > $dir/database/${species}.mimic.txt
#cp $dir/mirna/psrobot.mirna.mimic.name.sort $dir/database/${species}.mimic.lncrna
#awk -F'\t' '{OFS="\t"; print $2,$1}' $dir/mirna/psrobot.mirna.mimic.name |sort |awk -F'\t' 'NF>1{a[$1] = a[$1]","$2}END{for(i in a){print i"\t"a[i]}}'  |sed 's/\t,/\t/g' |sort > $dir/database/${species}.mimic.mirna
grep -Fwf $dir/heatmap/exp.sort <(cut -f 1,2 $dir/mirna/psrobot.mirna.target.name) |sed 's/>//g' > $dir/database/${species}.target.txt
#cp $dir/mirna/psrobot.mirna.target.name.sort $dir/database/${species}.target.lncrna
#awk -F'\t' '{OFS="\t"; print $2,$1}' $dir/mirna/psrobot.mirna.target.name |sort |awk -F'\t' 'NF>1{a[$1] = a[$1]","$2}END{for(i in a){print i"\t"a[i]}}'  |sed 's/\t,/\t/g' |sort > $dir/database/${species}.target.mirna

############## repeat transposon origin#########
grep -Fwf $dir/heatmap/exp.sort $dir/${name}-CPC_left.fa.out |awk -F" " '{OFS="\t"; print$5,$6,$7,$10,$11}' > $dir/database/${species}.repeat.txt

 #############tissue specific ###################
grep -Fwf $dir/heatmap/exp.sort <(awk -F'\t' 'NF>1{a[$1] = a[$1]","$2}END{for(i in a){print i"\t"a[i]}}' $dir/heatmap/normlnc_tissue.oneline.txt) |sed 's/\t,/\t/g' > $dir/database/${species}.tissue.txt

#####expression#############
grep -Fwf $dir/heatmap/exp.sort $dir/heatmap/normlnc_expression_mean.txt > $dir/database/${species}.expression.txt.tmp1
sed -n '1p' $dir/heatmap/normlnc_expression_mean.txt > $dir/database/${species}.expression.txt.tmp2
cat $dir/database/${species}.expression.txt.tmp2 $dir/database/${species}.expression.txt.tmp1 > $dir/database/${species}.expression.txt.tmp
Rscript /psc/bioinformatics/sunyd/identify/script/float_keep.R $dir/database/${species}.expression.txt.tmp $dir/database/${species}.expression.txt 
rm $dir/database/${species}.expression.txt.tmp*

#####coexpression#############
grep -Fwf $dir/heatmap/exp.sort <(awk -F "\t" '{OFS="\t"; print $1,$2,$1,$3,$1,$4,$1,$5,$1,$6,$1,$7,$1,$8,$1,$9,$1,$10,$1,$11,$1,$12,$1,$13,$1,$14,$1,$15,$1,$16,$1,$17,$1,$18,$1,$19,$1,$20,$1,$21}' $dir/coexpression/go.correlation.txt) |sed 's/TCONS/\nTCONS/g' |sed '/^$/d' |sed 's/\t$//g' > $dir/coexpression/go.correlation.2col.txt
grep -Fwf $dir/coexpression/go.correlation.2col.txt $dir/coexpression/go.correlation.value.each.corres.txt |sort -k2 > $dir/database/${species}.each.corres.txt
cut -f 2-4 $go |sort -u > $dir/database/go.tmp
join -a1 -1 2 -2 1 $dir/database/${species}.each.corres.txt $dir/database/go.tmp |awk -F" " '{OFS="\t"; print$2,$1,$3,$4,$5}' |sort -k1,1 -k3,3gr > $dir/database/${species}.each.corres.txt1
Rscript /psc/bioinformatics/sunyd/identify/script/pearson_6power.R $dir/database/${species}.each.corres.txt1 $dir/database/${species}.each.corres.txt
rm $dir/database/${species}.each.corres.txt1
#cp $dir/coexpression/go.correlation.topGO.txt $dir/database/${species}.topgo.txt
#cat $dir/coexpression/normlnc_expression_gene.high.go.txt $dir/coexpression/normlnc_expression_lnc.txt |sort -u > $dir/database/${species}.coexpressed.expression
rm $dir/database/go.tmp

#######put all txt file into a dir########
mkdir $dir/database/txt
mv $dir/database/${species}.* $dir/database/txt

#####RNA structure information##########
#/psc/program/install/ViennaRNA/bin/RNAfold --MEA -d2 -p < $dir/${name}-CPC_left.lncname.iuox.del.fa > ../${species}.fold.txt
#mkdir $dir/database/struc_eps
#mkdir $dir/database/struc_png
#for i in `cut -f1 $dir/database/*.lncrna.oneline.fa |cut -d_ -f2,3`; do cp $dir/stru/${i}_ss.ps $dir/database/struc_eps; cp $dir/stru/${i}_ss-1.png $dir/database/struc_png; done
#rename '-1' '' struc_png/*-1.png
#tar czf $dir/database/struc_eps.tar.gz $dir/database/struc_eps/
#tar czf $dir/database/struc_png.tar.gz $dir/database/struc_png/

##for i in $dir/database/struc_eps/*; do j=`echo $i |rev |cut -d/ -f1 |rev |sed 's/ps/pdf/g'`; z=`echo $i |rev |cut -d/ -f1 |rev |sed 's/\.ps//g'`; ps2pdf $i $dir/database/struc/$j; pdftoppm -png $dir/database/struc/$j $dir/database/struc_png/$z; done

