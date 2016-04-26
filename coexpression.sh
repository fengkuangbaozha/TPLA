dir=$1
ref=$2

#######use fpkm seperately#############
for i in ${dir}/cuffnorm.*/isoforms.fpkm_table; do j=`echo $i |sed 's/table/table.transpose/g'`; (head -n 1 $i && tail -n +2 $i |sort -u) |/psc/home/sunyidan/scripts/transpose.awk > $j; done

cat <(cat ${dir}/cuffnorm.*/isoforms.fpkm_table.transpose |sed -n '/tracking_id/p' |head -n 1) <(cat ${dir}/cuffnorm.*/isoforms.fpkm_table.transpose |sort -u |sed '/tracking_id/d') > ${dir}/tmp  #####vim tmp and make the trackingid the first line
mkdir -p $dir/coexpression
/psc/home/sunyidan/scripts/transpose.awk ${dir}/tmp > ${dir}/coexpression/isoforms.fpkm_table

cut -f1 ${dir}/tmp |sed '1d' |rev |cut -d_ -f 2,3,4 |rev > ${dir}/coexpression/samples.table

Rscript /psc/bioinformatics/sunyd/identify/script/fpkm_high.R ${dir}/coexpression/isoforms.fpkm_table ${dir}/coexpression/samples.table ${ref}.proteincoding.name ${dir}/cuff*-CPC_left.lncname.iuox.del ${dir}/coexpression/normlnc   ###get the high expression of all genes,lncRNAs


######Single lncrna coexpression result############
cat <(sed -n '1p' ${dir}/coexpression/normlnc_expression_gene.high.txt) <(grep -Fwf <(cut -f2 ${ref}.ano.all |sort -u) ${dir}/coexpression/normlnc_expression_gene.high.txt) > ${dir}/coexpression/normlnc_expression_gene.high.go.txt    #####get the GO annotation gene name
Rscript /psc/bioinformatics/sunyd/identify/script/coexpression_power_gene_out.R ${dir}/coexpression/normlnc_expression_gene.high.go.txt ${dir}/coexpression/normlnc_expression_lnc.txt ${dir}/coexpression/go. 6      ###the last one is softpower and correlation genes for each lncrna  and perform it twice because the softpower need to be adjusted

