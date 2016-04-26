############## analyze fruit ripen different expression#############
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/arab_rnaseq/SRRback/cuffdiff.cold cold Shoot_cold_control Shoot_cold_treatment
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/arab_rnaseq/SRRback/cuffdiff.heat heat Seedling_heat_wt_control Seedling_heat_wt_treatment
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.drought drought Root_24h_control Root_24h_severe_drought
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress cold Seedling_cold Seedling_control
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress cold Seedling_cold Seedling_control
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress salt Seedling_salt Seedling_control
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress heat Seedling_heat Seedling_control
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress uv Seedling_UV Seedling_control
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/oryza_rnaseq/SRRback/cuffdiff.piroot pi Root_pistar_21d Root_pisuff_21d
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/oryza_rnaseq/SRRback/cuffdiff.salt salt.root Root_control Root_salt
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/oryza_rnaseq/SRRback/cuffdiff.salt salt.shoot Shoot_control Shoot_salt
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.drought drought.6h Shoot_Benning_drought_0h Shoot_Benning_drought_6h
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.drought drought.12h Shoot_Benning_drought_0h Shoot_Benning_drought_12h
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.drought drought.24h Shoot_Benning_drought_0h Shoot_Benning_drought_24h
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.salt salt.1h Root_control Root_salt_1h
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.salt salt.6h Root_control Root_salt_6h
/psc/program/install/R-3.1.1/bin/Rscript ~/sunyd/identify/script/cumme.R ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.salt salt.12h Root_control Root_salt_12h
cat ~/sunyd/identify/arab_rnaseq/SRRback/cuffdiff.cold/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Cold"}'  >> ~/sunyd/identify/arab_rnaseq/SRRback/database/ath.condition.txt
cat ~/sunyd/identify/arab_rnaseq/SRRback/cuffdiff.heat/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Heat"}'  >> ~/sunyd/identify/arab_rnaseq/SRRback/database/ath.condition.txt
cat ~/sunyd/identify/oryza_rnaseq/SRRback/cuffdiff.salt/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Salt"}'  >> ~/sunyd/identify/oryza_rnaseq/SRRback/database/osa.condition.txt
cat ~/sunyd/identify/oryza_rnaseq/SRRback/cuffdiff.piroot/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Pi"}'  >> ~/sunyd/identify/oryza_rnaseq/SRRback/database/osa.condition.txt
cat ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.salt/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Salt"}' >> ~/sunyd/identify/soybean_rnaseq/SRRback/database/gma.condition.txt
cat ~/sunyd/identify/soybean_rnaseq/SRRback/cuffdiff.drought/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Drought"}' >> ~/sunyd/identify/soybean_rnaseq/SRRback/database/gma.condition.txt
cat ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.drought/*isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "Drought"}' >> ~/sunyd/identify/maize_rnaseq/SRRback/database/zma.condition.txt
cat ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress/cold.isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "cold"}' >> ~/sunyd/identify/maize_rnaseq/SRRback/database/zma.condition.txt
cat ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress/salt.isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "salt"}' >> ~/sunyd/identify/maize_rnaseq/SRRback/database/zma.condition.txt
cat ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress/heat.isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "heat"}' >> ~/sunyd/identify/maize_rnaseq/SRRback/database/zma.condition.txt
cat ~/sunyd/identify/maize_rnaseq/SRRback/cuffdiff.stress/uv.isodiff.txt |cut -f1 |sort -u |sed -n '/TCONS/p' |awk '{OFS="\t"; print $1, "UV"}' >> ~/sunyd/identify/maize_rnaseq/SRRback/database/zma.condition.txt
for i in /psc/bioinformatics/sunyd/identify/arab_rnaseq/SRRback/database/*.condition.txt; do j=`echo $i |cut -d "/" -f 1-6`; grep -Fwf $j/SRRback/heatmap/exp.sort $i |sed 's/TCONS_/Ath_TCONS_/g' > $j/SRRback/database/txt/ath.condition.txt; rm ${i}; done
for i in /psc/bioinformatics/sunyd/identify/oryza_rnaseq/SRRback/database/*.condition.txt; do j=`echo $i |cut -d "/" -f 1-6`; grep -Fwf $j/SRRback/heatmap/exp.sort $i |sed 's/TCONS_/Osa_TCONS_/g' > $j/SRRback/database/txt/osa.condition.txt; rm ${i}; done
for i in /psc/bioinformatics/sunyd/identify/soybean_rnaseq/SRRback/database/*.condition.txt; do j=`echo $i |cut -d "/" -f 1-6`; grep -Fwf $j/SRRback/heatmap/exp.sort $i |sed 's/TCONS_/Gma_TCONS_/g' > $j/SRRback/database/txt/gma.condition.txt; rm ${i}; done
for i in /psc/bioinformatics/sunyd/identify/maize_rnaseq/SRRback/database/*.condition.txt; do j=`echo $i |cut -d "/" -f 1-6`; grep -Fwf $j/SRRback/heatmap/exp.sort $i |sed 's/TCONS_/Zma_TCONS_/g' > $j/SRRback/database/txt/zma.condition.txt; rm ${i}; done
