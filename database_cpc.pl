#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
#--------------------------configuration of external programs-------------------------------------------------
my $ncbi_blast="/psc/program/install/blast/blast-2.2.26/bin/blastall"; # could be modified according to users' computational environment
my $cpc_home="/psc/program/install/cpc-0.9-r2";#could be modified according to users' computational environment
my $ViennaRNA="/psc/program/install/ViennaRNA/bin"; # could be modified according to users' computational environment
my $hmmerscan="/psc/program/install/hmmer/bin/hmmscan"; # could be modified according to users' computational environment
my $inferscan="/psc/program/install/infernal-1.1.1/cmscan"; # could be modified according to users' computational environment
#----------------------------bowtie1---------------------------
#my $ncbi_blast_formatdb="/psc/program/install/blast/blast-2.2.26/bin/formatdb"; # could be modified according to users' computational environment
#my $bowtie_program="/soft/bowtie/1.0.0/bowtie"; #could be modified according to users' computational environment
#my $bowtie_build_program="/soft/bowtie/1.0.0/bowtie-build"; #could be modified according to users' computational environment
#----------------------------bowtie2---------------------------
#my $bowtie_program="/psc/program/install/bowtie2/bowtie-2.2/bowtie2"; #could be modified according to users' computational environment
#my $bowtie_build_program="/psc/program/install/bowtie2/bowtie-2.2/bowtie2-build"; #could be modified according to users' computational environment
#----------------------------main interface that the program starts from--------------------------------------
&main(); 

#-----------------------------------subroutine begins----------------------------------------------------------
sub main(){
	my $USAGE = qq(
Welcome to use LncRNA_Finder to identify long noncoding RNAs only according to the nucleotide attributes of input sequences!

USAGE:
   LncRNA_Finder.pl -scriptonly <RF_seqonly.R> -scriptstru <RF_seqstru.R> -i <transcript.fasta> -p <protein.fasta> -trainstru <training file> -trainonly <training seqonly file> -pfam <Pfam.hmm> -rfam <rfam.cm> -o <output prefix> [-t <# of thread>] [-r <minimum lncRNA length>] [-f <maximum ORF length>] [-e <E-value of alignment>]
Options:
   -scriptonly <RF_seqonly.R>
   -scriptstru <RF_seqstru.R>
   -i <transcript.fasta>
   -p <protein.fasta>
   -trainstru <training file>
   -trainonly <training seqonly file>
   -pfam <Pfam.hmm>
   -rfam <rfam.cm>
   -o <output prefix>
   -h help
   -t <int> number of thread for the computation || default=4
   -r <int> minimum lncRNA length || default=200
   -f <int> maximum potential ORF length of lncRNAs || default=100
   -m <int> number of mismatch in the alignment with smallRNA || default=0
   -e E-value of the alignment against protein database || default=1.0e-3
   
Requirement:
   Need to install ncbi_blast standalone package, bowtie and cpc in your local environment correctly and set the running path of external programs at the begining of the pipeline
   
Citation:
   Li L, Eichten SR, Shimizu R, Petsch K, Yeh CT, Wu W, Scanlon MJ, Yu JM, Schnable PS, Timmermans MCP, Springer NM, and Muehlbauer GJ: Genome-wide discovery and characterization of maize long non-coding RNAs (lncRNAs). Genome Biology, 2013, revised.
	);
	my ($inputfile,$profile,$outfile,$trainstru,$trainonly,$scriptonly,$scriptstru,$pfam,$rfam); #required parameters
	my ($mintslength,$maxorflength,$evalue,$thread);#optional parameters
	#
	my $parastr = &GetOptions("scriptonly=s{1}"=>\$scriptonly,
							  "scriptstru=s{1}"=>\$scriptstru,
							  "i=s{1}"=>\$inputfile,
							  "p=s{1}"=>\$profile,
							  "trainstru=s{1}"=>\$trainstru,
							  "trainonly=s{1}"=>\$trainonly,
							  "pfam=s{1}"=>\$pfam,
							  "rfam=s{1}"=>\$rfam,
							  "o=s{1}"=>\$outfile,
							  "t=i{0,1}"=>\$thread,
							  "r=i{0,1}"=>\$mintslength,
							  "f=i{0,1}"=>\$maxorflength,
							  "e=s{0,1}"=>\$evalue
	);
	GetOptions ("help|?");
#	system("mkdir $outfile/fasta");
	unless ($parastr && defined($outfile)) {
	   print $USAGE."\n";
	   exit 1;
	}else{
		print "----------------------------------------------------------\n\n";
		unless(-e $scriptstru){
			print "Randomforest predication model file-$scriptstru does not exist\n";
			die("Pipeline halted\n");
		}else{
			print "Randomforest predication model file: $scriptstru\n";
		}
		unless(-e $inputfile){
			print "input transcript sequence file-$inputfile does not exist\n";
			die("Pipeline halted\n");
		}else{
			print "input transcript sequence file: $inputfile\n";
		}
#		unless(-e $profile){
#			print "input protein database file-$profile does not exist\n";
#			die("Pipeline halted\n");
#		}else{
#			print "input protein database file: $profile\n";
#		}		
		#check the optional parameters
		unless(defined($thread)){
			$thread=4;
		}
		unless(defined($mintslength)){
			$mintslength=200;
		}
		unless(defined($maxorflength)){
			$maxorflength=100;
		}
		unless(defined($evalue)){
			$evalue=1.0e-3;
		}
		print "----------------------------------------------------------\n\n";
		print "LncRNA_Finder starts to work with the following parameters:\n";
		print "number of thread for computation:\t$thread\n";
		print "minimum lncRNA length:\t$mintslength\n";
		print "maximum potential ORF length of lncRNAs:\t$maxorflength\n";
		print "E-value of the alignment against protein database:\t$evalue\n";
		print "----------------------------------------------------------\n\n";
	}
		
	my %oriseqhash=getOriSeqinfo($inputfile);
	print "processing transcript length/ORF filter...\n";
	print "----------------------------------------------------------\n\n";
	#transcript length/ORF filter
	my $tslorf=$outfile."-RNAlen_ORFlen_info.txt";
	my $orf=$outfile."-orf.fa";
	my $lenfilterleftseqfile=$outfile."-lenfilterleft.fa";
	open LENOUT,">$tslorf" || die "Cannot create result file:$!";
	open ORF,">$orf" || die "Cannot create result file:$!";
	open LEFTOUT,">$lenfilterleftseqfile" || die "Cannot create result file:$!";
	print LENOUT "seqname\ttranscriptlen\tORF1\tORF2\tORF3\trevORF1\trevORF2\trevORF3\n";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		my $toutstr=$tkey;
		$toutstr.="\t".length($tvalue);
		my $checkstr=checkprotein(translate_frame($tvalue, 1));
		my $checkstr2=checkprotein(translate_frame($tvalue, 2));
		my $checkstr3=checkprotein(translate_frame($tvalue, 3));
		$toutstr.="\t".length($checkstr)."\t".length($checkstr2)."\t".length($checkstr3);
		my $revcom = revcom($tvalue);
		my $checkstr4=checkprotein(translate_frame($revcom, 1));
		my $checkstr5=checkprotein(translate_frame($revcom, 2));
		my $checkstr6=checkprotein(translate_frame($revcom, 3));
		$toutstr.="\t".length($checkstr4)."\t".length($checkstr5)."\t".length($checkstr6);
		print LENOUT $toutstr."\n";
		if(length($checkstr)>$maxorflength || length($checkstr2)>$maxorflength || length($checkstr3)>$maxorflength || length($checkstr4)>$maxorflength || length($checkstr5)>$maxorflength || length($checkstr6)>$maxorflength || length($tvalue)<$mintslength){
			$oriseqhash{$tkey}="Length_ORF";
		}else{
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
			print ORF ">".$tkey."\n".$checkstr."\n";
			print ORF ">".$tkey."\n".$checkstr2."\n";
			print ORF ">".$tkey."\n".$checkstr3."\n";
			print ORF ">".$tkey."\n".$checkstr4."\n";
			print ORF ">".$tkey."\n".$checkstr5."\n";
			print ORF ">".$tkey."\n".$checkstr6."\n";
		}
	}
	close(LEFTOUT);
	close(ORF);
	close(LENOUT);
	print "transcript length/ORF filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing protein db filter...\n";
	print "----------------------------------------------------------\n\n";
	######Hmmer scan #########
	my $prot_hmmer_file=$outfile."-prot_hmmer_res.txt";
	my $hmmercmd=$hmmerscan." --cpu $thread -E $evalue --tblout $prot_hmmer_file $pfam $orf";
	print $hmmercmd."\n";
	if(system($hmmercmd)!=0){
		die("Hmmer package has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	###parse the hmmer result###
	open HMMER,"$prot_hmmer_file" || die "Cannot open $prot_hmmer_file:$!";
	my $hmmerstr;
	while($hmmerstr=<HMMER>){
		if($hmmerstr !~ /^\#/){
		my @barr=split(/\s+/,trim($hmmerstr));
		$oriseqhash{trim($barr[2])}="HMMER";
		}
	}
	close(HMMER);
	my $hmmerleftseqfile=$outfile."-hmmer_left.fa";
	open LEFTOUT,">$hmmerleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne "HMMER" and $tvalue ne "Length_ORF"){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "HMMER filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing BLAST filter...\n";
	print "----------------------------------------------------------\n\n";
	#protein DB filter
#	my $formatdbcmd=$ncbi_blast_formatdb." -i $profile -p T -o T -n prot_db";
#	if(system($formatdbcmd)!=0){
#		die("NCBI BLAST package has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
#	}
#	my $profilenew
#	my $profilenew = $profile =~ s/(.+)\.fa/$1/
	my $prot_blast_file=$outfile."-prot_blast_res.txt";
	my $blastcmd=$ncbi_blast." -d $profile -i $hmmerleftseqfile -a $thread -p blastx -m 8 -e $evalue -o $prot_blast_file"; ####$profile ---- prot_db
	print $blastcmd."\n";
	if(system($blastcmd)!=0){
		die("NCBI BLAST package has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	#parse the protein DB blast result
	open BLASTIN,"$prot_blast_file" || die "Cannot open $prot_blast_file:$!";
	my $blaststr;
	while($blaststr=<BLASTIN>){
		my @barr=split(/\t/,trim($blaststr));
		$oriseqhash{trim($barr[0])}="BLAST";
	}
	close(BLASTIN);
	my $protdbleftseqfile=$outfile."-protdb_blast_left.fa";
	open LEFTOUT,">$protdbleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne "BLAST" and $tvalue ne "HMMER" and $tvalue ne "Length_ORF"){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "protein db filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing INFER filter...\n";
	print "----------------------------------------------------------\n\n";
	######Infernal Rfam scan #########
	my $prot_infer_file=$outfile."-prot_infer_res.txt";
	my $infercmd=$inferscan." --cpu $thread -E $evalue --tblout $prot_infer_file $rfam $protdbleftseqfile";
	print $infercmd."\n";
	if(system($infercmd)!=0){
		die("Infernal package has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	###parse the hmmer result###
	open INFER,"$prot_infer_file" || die "Cannot open $prot_infer_file:$!";
	my $inferstr;
	while($inferstr=<INFER>){
		if($inferstr !~ /^\#/){
		my @barr=split(/\s+/,trim($inferstr));
		$oriseqhash{trim($barr[2])}="INFER";
		}
	}
	close(INFER);
	my $inferleftseqfile=$outfile."-infer_left.fa";
	open LEFTOUT,">$inferleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne "INFER" and $tvalue ne "BLAST" and $tvalue ne "HMMER" and $tvalue ne "Length_ORF"){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "INFER filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing RandomForest filter...\n";
	print "----------------------------------------------------------\n\n";
	#########RNA structure calculation###
	my $fold=$outfile."-foldtxt";
	my $lfold=$outfile."-lfoldtxt";
	system("$ViennaRNA/RNAfold --MEA -d2 -p < $inferleftseqfile > $fold");
	system("$ViennaRNA/RNALfold -z < $inferleftseqfile > $lfold");
	system("rm *.ps");
	#Random forest model potential coding score calculation
	my $seqstru=$outfile."-rfmodel_seqstruall";
#	my $seqonly=$outfile."-rfmodel_seqstrutop";
	system("/psc/program/install/R-3.0.2/bin/Rscript $scriptstru $trainstru $inferleftseqfile $fold $lfold $seqstru");
#	system("/psc/program/install/R-3.0.2/bin/Rscript $scriptonly $trainonly $protdbleftseqfile $fold $lfold $seqonly");
	print "Randomforest model filter completed\n";
	print "----------------------------------------------------------\n\n";
	print "processing CPC filter...\n";
	print "----------------------------------------------------------\n\n";
	#cpc potential protein score calculation
	my $cpc_res_file=$outfile."-cpc_res_table.txt";
	my $cpc_result_evidence_file=$outfile."-cpc_result_evidence.txt";
	my $cpccmd=$cpc_home."/bin/run_predict.sh $inferleftseqfile $cpc_res_file ./ $cpc_result_evidence_file";
	print $cpccmd."\n";
	if(system($cpccmd)!=0){
		die("CPC program has not been installed properly to be invoked by lncRNA_Finder Pipeline\nPipeline\tHalted\n");
	}
	#parse the cpc result
	open CPCIN,"$cpc_res_file" || die "Cannot open $prot_blast_file:$!";
	my $cpcstr;
	while($cpcstr=<CPCIN>){
		my @barr=split(/\t/,trim($cpcstr));
		if(trim($barr[2]) eq "coding"){
			$oriseqhash{trim($barr[0])}="";
		}
	}
	close(CPCIN);
	my $CPCleftseqfile=$outfile."-CPC_left.fa";
	open LEFTOUT,">$CPCleftseqfile" || die "Cannot create result file:$!";
	while (my ($tkey,$tvalue)=each(%oriseqhash)){
		if($tvalue ne "" and $tvalue ne "INFER" and $tvalue ne "BLAST" and $tvalue ne "HMMER" and $tvalue ne "Length_ORF"){
			print LEFTOUT ">".$tkey."\n".$tvalue."\n";
		}
	}
	close(LEFTOUT);
	print "CPC filter completed\n";
	print "----------------------------------------------------------\n\n";
}
sub getOriSeqinfo(){
	my $orifile=shift;
	my %resarr;
	open ORI,$orifile or die "Cannot open $orifile:$!";
	my $tstr;
	while($tstr=<ORI>){
		if(trim($tstr) ne "" && $tstr=~/>/){
			my @tarr=split(/\s+/,trim($tstr));
			my $tseqstr="";
			my $nstr="";
			while ($nstr=<ORI>){
				unless($nstr =~ />/){
					$tseqstr .= trim($nstr);
				}
				if($nstr=~/>/ or eof(ORI) ){
					my $tname=substr(trim($tarr[0]),1);
					$resarr{$tname}=$tseqstr;
					seek(ORI,-length($nstr),1);
					last;
				}
			}
		}
	}
	close(ORI);
	return %resarr;
}
sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }#else{

#            print STDERR "Bad codon \"$codon\"!!\n";
#            exit;
#    }
}

sub dna2peptide {

    my($dna) = @_;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
    return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}

# revcom 
#
# A subroutine to compute the reverse complement of DNA sequence

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

sub checkprotein(){
	my $tstr=shift;
	my $resstr="";
	my @tarr=split(/\*/,trim($tstr));
	my $maxlen=0;
	my $maxstr="NA";
	my $tpro="";
	pop(@tarr);
	foreach $tpro(@tarr){
		if(trim($tpro) ne ""){
			if(index($tpro,"M")!=-1){
				if($maxlen<length(trim($tpro))-index(trim($tpro),"M")){
					$maxlen=length(trim($tpro))-index(trim($tpro),"M");
					$maxstr=substr(trim($tpro),index(trim($tpro),"M"));
				}
			}
		}
	}
	$resstr=$maxstr;
	return $resstr;
}

sub trim(){
	my $string=shift;
	$string=~s/^\s+//;
	$string=~s/\s+$//;
	return $string;
}
		
