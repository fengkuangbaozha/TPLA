use strict;
use Bio::SeqIO;

my $usage = "perl dnamotif.pl <fasta file> <motif>";
my $fasta_filename = shift(@ARGV) or die("Usage: $usage $!");
my $motif = shift(@ARGV) or die("Usage: $usage $!");

my $fasta_parser = Bio::SeqIO->new(-file => $fasta_filename, -format => 'Fasta');
while(my $seq_obj = $fasta_parser->next_seq())
{
#  printf("Searching sequence '%s'...", $seq_obj->id);
  if((my $pos = index($seq_obj->seq(), $motif)) != -1)
  {
    print($seq_obj->id,"\t",$motif,"\t", $pos + 1, "\n");
  }
  else
  {
#    printf("motif not found.\n");
  }
}
