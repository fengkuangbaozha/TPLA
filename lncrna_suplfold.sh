#$ -cwd
a=$1
b=$2
c=$3
arg_input_seq=${a}${SGE_TASK_ID}      ########the sample path and name 
arg_input_name=${b}${SGE_TASK_ID}     #######input sample name
arg_working_dir=${c}${SGE_TASK_ID}    #######working directory name

#arg_output_evd_plot_feat_base=$4
echo "$arg_input_name begin at "  
  date

script=/psc/bioinformatics/sunyd/identify/script

#############running RNAlfold ############################
#/psc/program/install/ViennaRNA/bin/RNALfold -z < $arg_input_seq > $arg_working_dir/${arg_input_name}.lfoldtxt
Rscript $script/struc_short.R $arg_working_dir/${arg_input_name}.foldtxt $arg_working_dir/${arg_input_name}.lfoldtxt $arg_working_dir/${arg_input_name}.foldformattxt $arg_working_dir/${arg_input_name}.lfoldformattxt
############manipulate the files and combine them together#################

for i in $arg_working_dir/*.lfoldformattxt; 
do                                         ####sort the file except the first line
	j=`echo $i |sed 's/txt/sort/'`;
	z=`echo $i |sed 's/txt/sort1/'`;
	(head -n 1 $i && tail -n +2 $i |sort -u) > $z 
	sed '$d' $z > $j
	rm $z
done

Rscript $script/format.R $arg_working_dir/${arg_input_name}.cpcsort $arg_working_dir/blastx.sort $arg_working_dir/${arg_input_name}.compsort $arg_working_dir/${arg_input_name}.foldformatsort $arg_working_dir/${arg_input_name}.lfoldformatsort $arg_working_dir/${arg_input_name}.result
Rscript $script/ftop.R $arg_working_dir/${arg_input_name}.result $arg_working_dir/${arg_input_name}.top        #######select the top 50 features
Rscript $script/fcaret_vnew.R $arg_working_dir/${arg_input_name}.top $arg_working_dir/${arg_input_name}.rfmodel

echo "$arg_input_name finished at "
  date

