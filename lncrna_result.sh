dir=$1
script=/psc/home/sunyidan/scripts

cat $dir/cuff_find*/*-protdb_blast_left.fa > $dir/cuff-protdb_blast_left.fa
cat $dir/cuff_find*/*-CPC_left.fa > $dir/cuff-CPC_left.fa
cat $dir/cuff_ml*/*.rfmodel > $dir/cuff.rfmodel
cut -f 1 $dir/cuff.rfmodel |sort -u > $dir/cuff.rfmodelname
sed -n '/>/p' $dir/cuff-CPC_left.fa |sed 's/>//g' |sort -u > $dir/cuff-CPC_left.name       ##the lncRNA name     
cp $script/Library.pm $dir
diff $dir/cuff.rfmodelname $dir/cuff-CPC_left.name |grep "^<" |sed 's/< />/g'> $dir/diff.txt  #########put the diff one into the cpc-left to further analysis
perl $script/convert_extract_sequence.pl $dir/diff.txt $dir/cuff-protdb_blast_left.fa > $dir/cuff-rfmodeldiff.fa   ###get the diff sequence
sed 's/>//g' $dir/diff.txt |cat >> $dir/cuff-CPC_left.name
sort -u $dir/cuff-CPC_left.name > $dir/cuff-CPC_left.name1
mv $dir/cuff-CPC_left.name1 $dir/cuff-CPC_left.name
cat $dir/cuff-rfmodeldiff.fa >> $dir/cuff-CPC_left.fa 
##grep -Fwf cuff-CPC_left.name ../cuffmerge.tracking |cut -f 2 |sort -u |wc      #########check how many genes they belongs

