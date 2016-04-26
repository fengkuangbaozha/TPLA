gene=$1   #####ful1.target.lncrna.all
chipgff=$2  ######../cuff228.expression.gene.chip.gff3
diff=$3  ####cuffdiff.028278.gene/ac_tf_pink.isodiff.txt
chip=$4  #####chip-chip/GSE49125_ful/ful2_peaks_mapToFeaturest.bed
col=$5   ######genename column number
echo "Usage: sh ~/sunyd/identify/script/lncrna_tf.target.nearest.gene.sh ful1.target.lncrna ../cuff228.expression.gene.chip.gff3 cuffdiff.028278.gene/ac_tf_pink.isodiff.txt chip-chip/GSE49125_ful/ful1_peaks_mapToFeaturest.bed"

cut -f $col ${gene}.all |sed '/NONE/d' |sort -u > ${gene}.genename
grep -Fwf ${gene}.genename $chipgff > ${gene}.genename.chip.gff3
grep -Fwf ${gene}.genename $diff |cut -f 1,5-7 |sort |awk -F"\t" '{OFS="\t"; print $1,$2"/"$3"/"$4}' > ${gene}.genename.expression    ########find the differential expression
bedtools intersect -a ${gene}.genename.chip.gff3 -b $chip |cut -d ";" -f 1 |cut -f3,9 |sed 's/ID=//g' |sed 's/Parent=//g' |awk -F"\t" '{OFS="\t"; print $2,$1}' |sort -u|awk -F'\t' 'NF>1{a[$1] = a[$1]","$2}END{for(i in a){print i"\t"a[i]}}' |sed 's/\t,/\t/g' |sort > ${gene}.genename.targetloc  ######find the target location of ful1 

