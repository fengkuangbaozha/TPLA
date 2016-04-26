#!/usr/bin/perl
use strict;use warnings;
my ($num,$i);
my $a=shift;
my %go;

open (FH,$a) or die $!;
$num=0;
my$c=0;
while (<FH>) {
	
	chomp;
	my @a =split("\t",$_);
    if (!exists($go{$a[0]})){
		$go{$a[0]} =$a[1];
		#	$num++;
		}else{
			#	$go{$a[0]} .="\t$a[1]";
#	$i++;	
	}

}

close FH;
my @rest=@ARGV;
#print STDERR "Total annotated genes\t$num\nTotal Go Ids:\t",($num+$i),"\n";	
 foreach my $id(@rest){
		open (FH,$id) or die $!;
		my @b;
     while (<FH>) {
	if(/GO\.ID/){
		chomp;
		push @b,$_;}
	else{
		chomp;
	my @a =split("\t",$_);
    if (exists$go{$a[0]}){
    $a[1]=$go{$a[0]};
    }
    my $b=join("\t",@a);
    push @b,$b;
    }
} 
close FH;

	open (FH2, ">", "$id\.1") or die $!;
    foreach my $b(@b){
	print FH2 "$b\n"; 	
}

		
}
