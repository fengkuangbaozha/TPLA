#!/usr/bin/perl
use warnings; use strict;
use Getopt::Std;
use File::Basename;

my $usage = "$0 
		-r : <out_DIR> running directory default current 
		-n : <sample name> 
		-a : <'read_files can be globed by #'> 
		-s : single?  <T|F>  default F
		-p : <phred> 
		-l : <lenth> <5
		-5 : 'adaptor:-b Index> 
		-d : if run dynamic trim ?  <T|F>  default:T
		-c : if run cutadapter ?    <T|F>  default:T \n";

########
		#  DIR output directory default 'current directory'
		#



###############  get parameters #################
die "$usage" if (@ARGV == 0);
my %opt;
getopts("r:n:a:s:p:l:5:3:d:c:h",\%opt);

die  $usage if ($opt{h});

my $dyna_b = ($opt{d}?$opt{d}:"T");
my $cuta_b = ($opt{c}?$opt{c}:"T");

print STDERR "$dyna_b\t$cuta_b\n";

my($dir,$sample_name,$rds,$sin,$phred,$len,$ad_1) = ($opt{r},$opt{n},$opt{a},$opt{s},$opt{p},$opt{l},$opt{5});
$dir = $dir? $dir :$ENV{PWD};
$sin = $sin? $sin : "F";
my @reads;
if ($rds =~ s/#/\*/g){
	@reads = glob($rds);
}else{
	@reads = split /\s/, $rds;
}
print "READS: $rds\n@reads\n";
############   read adapter     #################

#my %adaptors;
#while (<DATA>){
#	chomp;
#	my ($id,$seq) = split /\t/,$_;
#	$adaptors{$id} = $seq;
#}
#my $uni   = $adaptors{$ad_1};
#$uni =~ tr/ATCGatcg/TAGCtagc/;
#$uni = reverse $uni;

#my $index = $adaptors{$ad_2};
#$index = "A".$index;

#print STDERR "Adaptor:\n$uni\n$index\n";

######## making your main result directory!!! ########


print STDERR "\n############making dir: $dir/$sample_name\n";
my $wd = "$dir/$sample_name";

mkdir $wd;

my $cwd = $ENV{PWD};
my $file_sta = "$wd/$sample_name.sta";


########unzipping the original gz files########
print STDERR "\n############unzipping fastq files:\n\t@reads\n.....\n";
my @raws;

foreach my $r (@reads){
	my $f_n = basename $r;
	my $cat;
	if ($f_n =~ /gz$/){
	   $cat = "zcat";
	}else{
	   $cat ="cat";
	}
	system ("printf '$f_n\tRaw:' >> $file_sta") == 0 or die $!;
	print STDERR "$cat $r |count_line_num.pl $file_sta >$wd/$f_n.fq\n";
	system ("$cat $r |count_line_num.pl $file_sta >$wd/$f_n.fq") == 0  or die $!;
	system ("echo >>$file_sta") ;
	push @raws, "$wd/$f_n.fq";
}



#######    dynamic triming_step_1   ###################
my $len_dir = "$wd/Len_sort_step_2";
my $sta_tri ;  #### determine if count base and reads for cut
my @trims;
if($dyna_b =~ /T/g){
	$sta_tri = $file_sta;
	my $dyn_dir = "$wd/Dynamic_dir_step_1";
	mkdir "$dyn_dir";
	print STDERR "\n############dynamic triming reads1 @raws \n";
	my @Dys;
	foreach my $r (@raws){
		my $n = basename $r;
		print STDERR "DynamicTrim $r -h $phred -d $dyn_dir\n";
		system ("DynamicTrim $r -h $phred -d $dyn_dir") == 0 or die $!;
		my $r_trimmed = "$dyn_dir/$n.trimmed";
		push  @Dys , $r_trimmed;
	}

	############### Length_sorting from dynamic trim ###########
	mkdir "$len_dir";
	if ($sin =~ /T/){
		print STDERR "\n###########lenth sorting step 2############\n";
		foreach (@Dys){
			system ("LengthSort $_ -l $len -d $len_dir") == 0 or die $!;
		}
		@trims = glob("$len_dir/*.single");
	}else{
		print STDERR "\n###########lenth sorting step 2############\n";
		system ("LengthSort @Dys -l $len -d $len_dir") == 0 or die $!;
		@trims = glob ( "$len_dir/*.paired?");

	}
	system ("rm -rf @Dys" ) == 0 or die $!;

}

#system ("rm -f $fq_1") == 0 or die $!;
#system ("rm -f $fq_2") == 0 or die $!;


#######    cut adaptor   !!!#######

my $len_dir_4 = "$wd/Len_sort_step_4";
mkdir "$len_dir_4";
my @cuts;
if ($cuta_b =~ /T/g){
	my $cut_dir = "$wd/Cutadt_dir_step_3";
	mkdir $cut_dir;
#chdir $cut_dir;
	my @cs;
	foreach  (@trims){
		system ("printf 'Trim:' >>$sta_tri") == 0 or die $!;
		my $n = basename $_;
		print STDERR "###########  cut adaptor ###########\n";
		my $arg_adap;
			$arg_adap = " -b $ad_1";
#		if ($_ =~ /paired2*/){
#			$arg_adap = "-a $uni -b $index";
#		}else{
#			$arg_adap = "-a $index -b $uni";
#		}

		print STDERR "cat $_ | count_line_num.pl $sta_tri  | cutadapt $arg_adap  -f fastq - -o  $cut_dir/$n.cut >$cut_dir/$n.sta\n";
		system ("cat $_ | count_line_num.pl $sta_tri  | cutadapt $arg_adap  -f fastq - -o  $cut_dir/$n.cut >$cut_dir/$n.sta") == 0 or die $!;
		system ("echo >> $sta_tri");
		push @cs, "$cut_dir/$n.cut"
	}

	############## Length sort #########
	print STDERR "\n#############length sorting\n";
	if ($sin =~ /T/){
	   foreach (@cs){
		   my $n = basename $_;
		   print STDERR "LengthSor t $_ -l $len -d $len_dir_4\n" or die $!;

		   system ("LengthSort $_ -l $len -d $len_dir_4") == 0 or die $!;
		   system ("printf 'Adaptor:' >>$sta_tri") == 0 or die $!;
		   system ("cat $len_dir_4/$n.single |  count_line_num.pl $file_sta NO ") == 0 or die $!;
		   system ( "echo >>$file_sta");
		   system ("cp $len_dir_4/$n.single $dir") == 0 or die $!;
	   }
	   @cuts = glob("$len_dir_4/*.single");
	   
   }else{
	   print STDERR "LengthSort @cs  -l $len -d $len_dir_4\n" or die $!;
	   system ("LengthSort @cs  -l $len -d $len_dir_4") == 0 or die $!;
	   @cuts = glob("$len_dir_4/*.paired?");
	   foreach my $f (@cuts){
		   system ("printf 'Adaptor:' >>$sta_tri") == 0 or die $!;
		   system ("cat $f |  count_line_num.pl $file_sta NO ") == 0 or die $!;
		   system ( "echo >>$file_sta");
	   }
   }
   system ("rm -rf @cs" ) == 0 or die $!;

}

print "@cuts\n";

system ("echo  >>$file_sta") == 0 or die $!;  #### change to new line

system ( "rm -rf @raws @trims" ) == 0 or die $!;


##########  compress file  ######
__DATA__
Uni	AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Index	GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
Index_1	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
Index_2	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
Index_3	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
Index_4	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
Index_5	GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
Index_6	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
Index_7	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
Index_8	GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
Index_9	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
Index_10	GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
Index_11	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
Index_12	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
Index_13	GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG
Index_14	GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG
Index_15	GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG
Index_16	GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG
Index_18	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG
Index_19	GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG
