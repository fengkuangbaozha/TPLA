#$ -cwd
################command lines####################
arg_working_dir=$1
arg_input_name=$2
gff=$3
genome=$4
lncgff=$5
htseq=/psc/program/install/python-2.7.6/bin/htseq-count

echo "$arg_input_name begin at "  
    date

#########################assemble the reads into transcripts####################
export LD_LIBRARY_PATH=/psc/program/install/boost-1.55/lib:$LD_LIBRARY_PATH
tophat -p 32 -G $gff -o $arg_working_dir/${arg_input_name}_topall $genome $arg_working_dir/$arg_input_name/Len_sort_step_4/${arg_input_name}_1.fastq.fq.trimmed.single.cut.single  #$file/SRR331227_1.fastq.fq.trimmed.single.cut.single1 
samtools sort -n $arg_working_dir/${arg_input_name}_topall/accepted_hits.bam  $arg_working_dir/${arg_input_name}_topall/accepted_hits.sorted
samtools view -h $arg_working_dir/${arg_input_name}_topall/accepted_hits.sorted.bam > $arg_working_dir/${arg_input_name}_topall/accepted_hits.sorted.sam 
$htseq -i ID -t gene -m union -s no -q $arg_working_dir/${arg_input_name}_topall/accepted_hits.sorted.sam $lncgff >> $arg_working_dir/${arg_input_name}.htseq
#$tool/rsem-calculate-expression -p 30 --paired-end ${sample}_R1.fq ${sample}_R2.fq $genome $wd3/rep${SGE_TASK_ID}

##$cufflink/cufflinks -o $arg_working_dir/${arg_input_name}_cuff -p 32 -g $gff $arg_working_dir/${arg_input_name}_top/accepted_hits.bam         #####assmeble the reads using cufflinks
#gtf_to_fasta $arg_working_dir/${arg_input_name}_cuff/transcripts.gtf ${genome}.fa $arg_working_dir/${arg_input_name}.fa             #####convert the cufflinks gtf file into fasta
#sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $arg_working_dir/${arg_input_name}.fa > $arg_working_dir/${arg_input_name}_cuff.fa
#mkdir $arg_working_dir/split_cuff
#split -a 3 -d $arg_working_dir/${arg_input_name}_cuff.fa $arg_working_dir/split_cuff/cuff
#for i in $arg_working_dir/split_cuff/cuff*0; do j=`echo $i |sed 's/0/zero/'`; mv $i $j; done
#for i in $arg_working_dir/split_cuff/cuff0*; do j=`echo $i |sed 's/0/zero/'`; mv $i $j; done
#mv $arg_working_dir/split_cuff/cuffzero0 $arg_working_dir/split_cuff/cuffzero1zero
#cd $arg_working_dir/split_cuff/
#rename '0' '' cuff0*
#rename '0' '' cuff0*
#cd $arg_working_dir/
########################another way for doing assemble using trinity###################
#$trinity/Trinity --seqType fq --JM 1500G --min_kmer_cov 2 --output $arg_working_dir/${arg_input_name}_trinity --single $arg_working_dir/$arg_input_name/Len_sort_step_4/${arg_input_name}_1.fastq.fq.trimmed.single.cut.single --CPU 60          #--left $rep1_1L1F  --right $rep1_1L1   --grid_conf $scripts/trinity_conf.txt
#cp $arg_working_dir/${arg_input_name}_trinity/Trinity.fasta $arg_working_dir/
#sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $arg_working_dir/${arg_input_name}_trinity/Trinity.fasta > $arg_working_dir/${arg_input_name}_trin.fa
#mkdir $arg_working_dir/split_trin
#split -d $arg_working_dir/${arg_input_name}_trin.fa $arg_working_dir/split_trin/trin
#cd $arg_working_dir/split_trin
#rename '0' '' $arg_working_dir/split_trin/trin0*
#for i in $arg_working_dir/split_trin/trin*0; do j=`echo $i |sed 's/0/zero/'`; mv $i $j; done
#for i in $arg_working_dir/split_trin/trin0*; do j=`echo $i |sed 's/0/zero/'`; mv $i $j; done
#mv $arg_working_dir/split_trin/trinzero0 $arg_working_dir/split_trin/trinzero1zero
echo "$arg_input_name finished at "
    date

