dir=$1
gtf=$2

###########expression, tissue specific and heatmap analysis############
cut -d_ -f 1 $dir/tissue.merge.cuffnorm/samples.table |sed '1d' > $dir/tissue.merge.cuffnorm/samples.table.group.com
mkdir -p $dir/heatmap
Rscript ~/sunyd/identify/script/fpkm_high.R $dir/tissue.merge.cuffnorm/isoforms.fpkm_table $dir/tissue.merge.cuffnorm/samples.table.group.com ${gtf}.proteincoding.name $dir/cuff*-CPC_left.lncname.iuox.del $dir/heatmap/normlnc   ###get the expression of all genes,lncRNAs
Rscript ~/sunyd/identify/script/tissuespecific.R $dir/heatmap/normlnc_expression_lnc.txt $dir/heatmap/normlnc_tissue.txt
awk -F" " '{OFS="\t"; print $1,$2,$1,$3,$1,$4,$1,$5}' $dir/heatmap/normlnc_tissue.txt |sed 's/TCONS/\nTCONS/g' |awk 'NF > 1' |sort > $dir/heatmap/normlnc_tissue.oneline.txt

############combine the result info together##################
cat $dir/split_cuff/cuff_find*/*-cpc_res_table.txt |sort -u |cut -f1,4 > $dir/heatmap/cpc.score
cat $dir/split_cuff/cuff_find*/*-rfmodel_seqstruall_score.txt |sort -u |cut -f1,2 > $dir/heatmap/rfmodel.score
sed '1d' $dir/coexpression/normlnc_expression_lnc.high.txt |cut -f1 |sort -u > $dir/heatmap/exp.sort     #####check if these two are the same with cuff*-CPC_left.lncname.iuox.del.name

join -a1 -1 5 -2 2 -t $'\t' -e "NONE" <(sort -k5 $dir/cuff*-CPC_left.lncname.iuox.del) <(sort -k2 ${gtf}.ano.all) |awk -F"\t" '{OFS="\t"; print $3,$2,$4,$1,$6,$7,$8}' |sed 's/\t\t\t/\tNONE\tNONE\tNONE/g' |sort -u > $dir/heatmap/normlnc.ano
grep -Fwf $dir/heatmap/exp.sort $dir/cuff*-CPC_left.gff3 |sed -n '/\ttranscript\t/p' |cut -d ";" -f1 |sed 's/ID=//g' |awk -F"\t" '{OFS="\t"; print $9,$1,$4,$5,$7}' |sort > $dir/heatmap/normlnc.location
tr '\n' '\t' < $dir/cuff*-CPC_left.lncname.iuox.del.fa |sed 's/\t>/\n/g' |sed 's/>//g'> $dir/heatmap/normlnc.fa   #####delete the last line tab delimiter
paste $dir/heatmap/normlnc.location <(grep -Fwf $dir/heatmap/exp.sort $dir/cuff*.sequence_info.txt |sort |cut -f2 ) <(grep -Fwf $dir/heatmap/exp.sort $dir/heatmap/rfmodel.score |sort |cut -f2 ) <(grep -Fwf $dir/heatmap/exp.sort $dir/heatmap/cpc.score |sort |cut -f2 ) <(grep -Fwf $dir/heatmap/exp.sort $dir/heatmap/normlnc.ano |sort |cut -f2-4,6-7) <(grep -Fwf $dir/heatmap/exp.sort $dir/heatmap/normlnc.fa |sort |cut -f2) |sed 's/ /_/g' |sort -k 6,6 -k 7,7nr > $dir/heatmap/lnc.all  ###<(sed '1d' ../coexpression/normlnc_expression_summary.txt |sort |cut -f5) 
cat <(sed -n '1p' $dir/coexpression/normlnc_expression_mean.txt) <(grep -Fwf $dir/heatmap/exp.sort $dir/coexpression/normlnc_expression_mean.txt |sort) > $dir/heatmap/lnc.expression

