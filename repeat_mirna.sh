dir=$1
name=$2
species=$3
mirna=$4

 ####psRobot Target##############
mkdir -p $dir/mirna
cut -d " " -f1 $mirna > $dir/mirna/mirna.fa
psRobot_tar -s $dir/mirna/mirna.fa -t $dir/${name}-CPC_left.lncname.iuox.del.fa -o $dir/mirna/psrobot.mirna.target -p 32
paste <(sed -n '/>/p' $dir/mirna/psrobot.mirna.target) <(sed -n '/Query:/p' $dir/mirna/psrobot.mirna.target) <(sed -n '/Sbjct/p' $dir/mirna/psrobot.mirna.target) <(sed -n '/|/p' $dir/mirna/psrobot.mirna.target) |awk -F"\t" '{OFS="\t"; print $3,$1}' |sed 's/>//g' |sort -u  > $dir/mirna/psrobot.mirna.target.name          #####only keep those 
    ####psRobot Target mimics##############
cat $dir/mirna/psrobot.mirna*mimic.name |sed 's/>//g' |cut -f1,2 |sort -u > $dir/mirna/psrobot.mirna.mimic.name

##############identify repeated element RNAs using RepeatMasker############
RepeatMasker -species $species $dir/${name}-CPC_left.fa      ###check the align quality
grep -Fwf <(awk -F " " '{OFS=" ";print $5}' $dir/${name}-CPC_left.fa.out |sort -u |sed '1,3d') $dir/${name}-CPC_left.gtf |awk '{        for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -F" " '{OFS="\t"; print $4,$2,$6,"TE"}' >  $dir/${name}_TE.name


