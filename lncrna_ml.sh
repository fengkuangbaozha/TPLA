#$ -cwd


arg_input_seq=$1      #######input sequence , can be absolute path
arg_input_name=$2     #######input sample name
arg_working_dir=$3    #######working directory name

#arg_output_evd_plot_feat_base=$4
echo "$arg_input_name begin at "  
  date

script=/psc/bioinformatics/sunyd/identify/script
#top=
###############caculating the sequence composition and codon bias ####################
mkdir $arg_working_dir
Rscript $script/comp_short.R $arg_input_seq $arg_working_dir/${arg_input_name}.comptxt

###############Structure related features#########################
cd $arg_working_dir
/psc/program/install/ViennaRNA/bin/RNAfold --MEA -d2 -p < $arg_input_seq > $arg_working_dir/${arg_input_name}.foldtxt
/psc/program/install/ViennaRNA/bin/RNALfold -z < $arg_input_seq > $arg_working_dir/${arg_input_name}.lfoldtxt
rm *.ps
#sed -n "/>/p" $arg_working_dir/${arg_input_name}.foldtxt |wc
#sed -n "/>/p" $arg_working_dir/${arg_input_name}.lfoldtxt |wc
Rscript $script/struc_short.R $arg_working_dir/${arg_input_name}.foldtxt $arg_working_dir/${arg_input_name}.lfoldtxt $arg_working_dir/${arg_input_name}.foldformattxt $arg_working_dir/${arg_input_name}.lfoldformattxt

###############cpc related features#################
source /psc/program/install/cpc-0.9-r2/profile
mkdir $arg_working_dir/cpc/
/psc/program/install/cpc-0.9-r2/bin/run_predict.sh $arg_input_seq $arg_working_dir/${arg_input_name}.cpctxt $arg_working_dir/cpc
cp $arg_working_dir/cpc/blastx.feat $arg_working_dir/blastx.txt
cut -f 1,2,4 $arg_working_dir/${arg_input_name}.cpctxt |sed '1i\RNAname\tsequence_length\tcpc_value' > $arg_working_dir/${arg_input_name}.cpctxt1
mv $arg_working_dir/${arg_input_name}.cpctxt1 $arg_working_dir/${arg_input_name}.cpctxt


############manipulate the files and combine them together#################
for i in $arg_working_dir/*; 
do                                         ####sort the file except the first line
	j=`echo $i |sed 's/txt/sort/'`;
	(head -n 1 $i && tail -n +2 $i |sort -u) > $j 
done
#sed -i '$ d' $arg_working_dir/${arg_input_name}.compsort
#sed -i '$ d' $arg_working_dir/${arg_input_name}.foldformatsort
#sed -i '$ d' $arg_working_dir/${arg_input_name}.lfoldformatsort
Rscript $script/format.R $arg_working_dir/${arg_input_name}.cpcsort $arg_working_dir/blastx.sort $arg_working_dir/${arg_input_name}.compsort $arg_working_dir/${arg_input_name}.foldformatsort $arg_working_dir/${arg_input_name}.lfoldformatsort $arg_working_dir/${arg_input_name}.result
Rscript $script/ftop.R $arg_working_dir/${arg_input_name}.result $arg_working_dir/${arg_input_name}.top        #######select the top 50 features
Rscript $script/fcaret_vnew.R $arg_working_dir/${arg_input_name}.top $arg_working_dir/${arg_input_name}.rfmodel

echo "$arg_input_name finished at "
  date

