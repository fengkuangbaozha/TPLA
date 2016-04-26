dir=$1
meth=$2
lncrna=$3
name=$4
location=$5
echo "Usage: sh methylation.intersect.sh ~/sunyd/identify/oryza_rnaseq/SRRpi /psc/bioinformatics/sunyd/identify/oryza_rnaseq/methylation/root/root.gtf /psc/bioinformatics/sunyd/identify/oryza_rnaseq/SRRpi/cuff165.expression.gene.gtf root genebody"

bedtools intersect -a $meth -b $lncrna -wa -wb > $dir/${name}.${location}.gtf
sed -n '/\tCG/p' $dir/${name}.${location}.gtf |cut -f 9,18 |awk '{OFS="\t"; print$4,$2,$3}' |sed 's/;//g' |sed 's/"//g' |sort > $dir/${name}.${location}.cg.meth
Rscript ~/sunyd/identify/script/methylation_level.R $dir/${name}.${location}.cg.meth $dir/${name}.${location}.cg.meth.level $dir/${name}.${location}.cg
sed -n '/\tCHG/p' $dir/${name}.${location}.gtf |cut -f 9,18 |awk '{OFS="\t"; print$4,$2,$3}' |sed 's/;//g' |sed 's/"//g' |sort > $dir/${name}.${location}.chg.meth
Rscript ~/sunyd/identify/script/methylation_level.R $dir/${name}.${location}.chg.meth $dir/${name}.${location}.chg.meth.level $dir/${name}.${location}.chg
sed -n '/\tCHH/p' $dir/${name}.${location}.gtf |cut -f 9,18 |awk '{OFS="\t"; print$4,$2,$3}' |sed 's/;//g' |sed 's/"//g' |sort > $dir/${name}.${location}.chh.meth
Rscript ~/sunyd/identify/script/methylation_level.R $dir/${name}.${location}.chh.meth $dir/${name}.${location}.chh.meth.level $dir/${name}.${location}.chh
