#$ -cwd
source /psc/program/install/cpc-0.9-r2/profile
#sample=/psc/bioinformatics/sunyd/identify/arab_rnaseq
sample=$1
#name=$2       ####the intial name of file
number=`echo "$SGE_TASK_ID-1" |bc`             #######never foget that, so it can read cuff0
script=/psc/bioinformatics/sunyd/identify/script
#dir=/psc/bioinformatics/sunyd/identify/oryza_rnaseq/SRR1677/
#for i in $sample/${name}*;
#for i in $sample/split_cuff/cuff$number; 
#do 
 #   j=`echo $sample |rev |cut -d/ -f 1 |rev`;	
	cd $TMP/;
	###################run maize lncRNAfinder.pl###################
	cp -r /psc/bioinformatics/sunyd/swiss_prot/ $TMP/
	mkdir $TMP/cuff_find$number
	cd $TMP/cuff_find$number;
	perl $script/database_cpc.pl -scriptonly $script/RF_seqstru.R -scriptstru $script/RF_seqstru_all.R -i $sample/split_cuff/cuff$number -p $TMP/swiss_prot/swiss_prot -trainstru $script/rf_model.top.seqstru.all -trainonly $script/rf_model.top.seqstru -pfam $TMP/swiss_prot/Pfam.hmm -rfam $TMP/swiss_prot/rfam/Rfam.cm.1_1 -o cuff$number -t 30
	cp -r $TMP/cuff_find$number $sample/split_cuff/
	
	#cd $TMP/;
	#mkdir $TMP/allcuff_find$number
	#cd $TMP/allcuff_find$number
	#perl $script/database_cpc.pl -scriptonly $script/RF_seqonly.R -scriptstru $script/RF_seqstru_all.R -i $sample/split_cuff/cuff$number -p $TMP/prot_db/prot_db -trainstru $script/rf_model.top.seqstru.all -trainonly $script/rf_model.top.seqonly -o allcuff$number  -t 30
	#perl $script/lncRNAfinder.pl -i $sample/split_cuff/cuff$number -p /psc/bioinformatics/sunyd/prot_db/prot_db -k /psc/bioinformatics/sunyd/lncrna/arab10.fa -s /psc/bioinformatics/sunyd/identify/arab_rnaseq/small.fa -o cuff$number  -t 30
	#cp -r $TMP/allcuff_find$number $sample/split_cuff/
#	mkdir $j/cuff_find$number;
#done
#sleep 5h
#$ -cwd
#train1=$1
#train2=$2
#sample=$3
#name=$2       ####the intial name of file
#number=`echo "$SGE_TASK_ID-1" |bc`             #######never foget that, so it can read cuff0
#echo "qsub -t 1 /psc/bioinformatics/sunyd/identify/lncrnaml_qsuball_newtrain.sh /psc/bioinformatics/sunyd/lncrna/ml_caret/seqonly_short_top50.feature ~/sunyd/lncrna/ml_caret/seqstru_top50.feature ~/sunyd/identify/arab_rnaseq/SRRall"
#################run my sequence info only rf model with the ORF and blast filter protein#######################
#Rscript /psc/bioinformatics/sunyd/identify/script/RF_seqonly.R $train1 $sample/split_cuff/cuff_find$number/cuff${number}-protdb_blast_left.fa $sample/split_cuff/cuff_seqonly$number

########run my sequence and structure rf model########
#Rscript /psc/bioinformatics/sunyd/identify/script/RF_seqstru.R $train2 $sample/split_cuff/cuff_find$number/cuff${number}-protdb_blast_left.fa $sample/split_cuff/cuff_ml$number/cuff${number}.foldsort $sample/split_cuff/cuff_ml$number/cuff${number}.lfoldsort $sample/split_cuff/cuff_seqstru$number








sleep 3h
