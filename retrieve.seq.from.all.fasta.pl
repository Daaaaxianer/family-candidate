#! /usr/bin/perl

use strict;
use Bio::SeqIO;

print "\n\tparameter assignment\n\tperl $0 idfile infastafile outfastafile\n";

my $idfile = $ARGV[0];
my @id = ();
open(ID, $idfile) or die "cannot open $idfile due to $!\n";
my $i = 0;
while(<ID>)
{
   $_ =~ s/[\n\r]//g;
   $id[$i] = $_;
   $i ++;
}

my $infastafile = $ARGV[1];
my $is = Bio::SeqIO -> new(-format=>'fasta', -file => $infastafile);

my $outfastafile = $ARGV[2];
my $os = Bio::SeqIO -> new(-format=>'fasta', -file => ">".$outfastafile);

while(my $seqobj = $is -> next_seq())
{
   my $id = $seqobj -> display_id();
   for(my $i=0; $i<=$#id; $i++)
   {
      if($id =~ /$id[$i]/)
      {
        $seqobj -> display_id($id[$i]);
        $os -> write_seq($seqobj); 
        last;
      }
   }
}
print "\tID sequence retrieved\n\n";

