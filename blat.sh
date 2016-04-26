dir=$1
number= `echo "${SGE_TASK_ID} - 1" |bc`

blat /psc/bioinformatics/sunyd/identify/script/liftover/ITAG2.4_genomic.2bit $dir/chr${number}.fa $dir/chr${number}.psl
blat /psc/bioinformatics/sunyd/identify/script/liftover.old/ITAG2_genomic.2bit $dir/chr${number}.fa $dir/chr${number}.old.psl -tileSize=12 -minScore=100 -minScore=100 -minIdentity=98
