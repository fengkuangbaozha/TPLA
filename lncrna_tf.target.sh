diff1=$1   #####cuffdiff.028278/ac_rin_pink.isodiff
diff2=$2   ##cuff228-CPC_left.lncname.iuox.del.name
target1=$3  ##cuff85-CPC_left.chip.gff3
target2=$4  ##chip-chip/GSE40257_rin/rin_peaks_mapToFeaturest.bed
dir=$5
name=$6
lncall=$7   ##../heatmap/lnc.all

###########for differentially expressed lncRNA in mutant and CK analysis######
cut -f 1 ${diff1}.txt |sort > ${diff1}.name    ###all the differentially expressed genes and lncRNAs
comm -1 -2 ${diff1}.name $diff2 |sort > ${diff1}.lncrna      ###only differentially expressed  lncRNAs
grep -Fwf ${diff1}.lncrna ${diff1}.txt |cut -f 1,5-7 |sort |awk -F"\t" '{OFS="\t"; print $1,$2"/"$3"/"$4}' > ${diff1}.lncrna.expression ####find the expression up/down regulated the lncRNAs

#############for TF target lncRNAs#########
bedtools intersect -a $target1 -b $target2 |sort -u  > $dir/${name}.lncrna   ####find the chip target location of the lncRNAs  
rev $dir/${name}.lncrna |cut -f 1 |rev |cut -d ";" -f 1 |cut -d "=" -f 2 |sort -u > $dir/${name}.lncrna.name
cut -d ";" -f 1 $dir/${name}.lncrna |rev |cut -d= -f 1 |rev > $dir/${name}.lncrna.targetloc.tmp1
cut -f 3 $dir/${name}.lncrna > $dir/${name}.lncrna.targetloc.tmp2
paste $dir/${name}.lncrna.targetloc.tmp1 $dir/${name}.lncrna.targetloc.tmp2 |sort -u |awk -F'\t' 'NF>1{a[$1] = a[$1]","$2}END{for(i in a){print i"\t"a[i]}}' |sed 's/\t,/\t/g' |sed 's/,transcript//g' |sort > $dir/${name}.lncrna.targetloc
rm $dir/${name}.lncrna.targetloc.tmp*
###########find both diff express and TF target lncRNA#########
comm -1 -2 ${diff1}.lncrna $dir/${name}.lncrna.name |sort > $dir/${name}.target.lncrna
grep -Fwf $dir/${name}.target.lncrna $lncall |sort > $dir/${name}.target.lncrna.info
grep -Fwf $dir/${name}.target.lncrna $dir/${name}.lncrna.targetloc |sort > $dir/${name}.target.lncrna.target
grep -Fwf $dir/${name}.target.lncrna ${diff1}.lncrna.expression |sort > $dir/${name}.target.lncrna.expression
paste $dir/${name}.target.lncrna.expression $dir/${name}.target.lncrna.target $dir/${name}.target.lncrna.info |sort > $dir/${name}.target.lncrna.all

#grep -Fwf ${name}.lncrna $correlation |sed 's/\t/\n/g' |sed '/TCONS/d' |sort -u > ${name}.genename   ##prepare genename for GO annotation
#/psc/program/src/R-3.1.1/bin/Rscript $script/topgo.coexpression.hclus.R $go ${name}.genename $name 0.05 0.05 0.05        ######annotate the genes
#sh ~/sunyd/identify/script/lnc.coexpress.GO.sh cuffdiff.028278/ac_rin_pink.isodiff.lncrna chip-chip/rin.lncrna.name coexpression/floral.correlation0.95.txt /psc/bioinformatics/sunyd/genome/tomato/ITAG2.4.go.merge.csv coexpression/ac_rin_rin
