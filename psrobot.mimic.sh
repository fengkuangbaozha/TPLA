#$ -cwd                                                                                                                                                                      
mirna=$1
fasta=$2
dir=$3
######/psc/bioinformatics/sunyd/genome/oryza/oryza_miRbase.fa
psRobot_tar -s $mirna -t $fasta -o $dir/psrobot.mirna.mimic${SGE_TASK_ID} -ts 5 -tp ${SGE_TASK_ID} -gl 8 -p 32 -gn 3
sed -n '/>/p' $dir/psrobot.mirna.mimic${SGE_TASK_ID} > $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp1
sed -n '/Query:/p' $dir/psrobot.mirna.mimic${SGE_TASK_ID} > $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp2
sed -n '/Sbjct/p' $dir/psrobot.mirna.mimic${SGE_TASK_ID} > $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp3
sed -n '/|/p' $dir/psrobot.mirna.mimic${SGE_TASK_ID} > $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp4
paste $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp* |awk -F" " '{OFS="\t"; print $4,$1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' |cut -f1,2,4,6,7,8,10,11,12,13 > $dir/psrobot.mirna.oneline${SGE_TASK_ID}
##paste $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp* |awk -F" " '{OFS="\t"; print $8,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15,$16,$17}' |cut -f1,2,4,6,7,8,10,11,12,13 > $dir/psrobot.mirna.oneline${SGE_TASK_ID}
rm $dir/psrobot.mirna.mimic${SGE_TASK_ID}.tmp*
#cut -f5 $dir/psrobot.mirna.oneline${SGE_TASK_ID} |sed -n '/---/p' > $dir/psrobot.mirna.candi${SGE_TASK_ID}.tmp
#grep -Fwf $dir/psrobot.mirna.candi${SGE_TASK_ID}.tmp $dir/psrobot.mirna.oneline${SGE_TASK_ID} > $dir/psrobot.mirna.candi${SGE_TASK_ID}
#rm $dir/psrobot.mirna.candi${SGE_TASK_ID}.tmp
Rscript ~/sunyd/identify/script/mirna_target_mimic.R $dir/psrobot.mirna.oneline${SGE_TASK_ID} $dir/psrobot.mirna${SGE_TASK_ID}

