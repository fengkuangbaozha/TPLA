#$ -cwd


arg_input_name=$1     #######input sample name
arg_working_dir=$2    #######working directory name

#arg_output_evd_plot_feat_base=$4
echo "$arg_input_name begin at "  
  date

script=/psc/bioinformatics/sunyd/identify/script


############manipulate the files and combine them together#################
for i in $arg_working_dir/*; 
do                                         ####sort the file except the first line
	j=`echo $i |sed 's/txt/sort/'`;
	(head -n 1 $i && tail -n +2 $i |sort -u) > $j 
done
Rscript $script/format.R $arg_working_dir/${arg_input_name}.cpcsort $arg_working_dir/blastx.sort $arg_working_dir/${arg_input_name}.compsort $arg_working_dir/${arg_input_name}.foldformatsort $arg_working_dir/${arg_input_name}.lfoldformatsort $arg_working_dir/${arg_input_name}.result
Rscript $script/ftop.R $arg_working_dir/${arg_input_name}.result $arg_working_dir/${arg_input_name}.top        #######select the top 50 features
Rscript $script/fcaret_vnew.R $arg_working_dir/${arg_input_name}.top $arg_working_dir/${arg_input_name}.rfmodel

echo "$arg_input_name finished at "
  date

