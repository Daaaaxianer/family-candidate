#! /usr/bin/perl
use strict;
use Bio::SeqIO;

### use pfam hmm seeds serach gene family candidates
### hmm_seed : hmm seed file of pfam
### protein.fa : fasta file or protein
### e-value1 : threshold of e-value to build specie hmm seed
### e-value1 : threshold of e-value to filter family candidate
### outname : outfile name for id file and fasta file 
### xian 2020-04 ###

die " perl $0 <hmm_seed> <protein.fa> <e-value1> <e-value2> <outname>" unless (@ARGV == 5);


my $seed = $ARGV[0];
my $pro = $ARGV[1];
my $evalue1 = $ARGV[2];
my $evalue2 = $ARGV[3];
my $name = $ARGV[4];

system ("hmmsearch --domtblout seed_domain_search.txt $seed $pro");
system ("perl refine_domain_seq.pl seed_domain_search.txt $pro $name.domain.fasta $evalue1");
system ("muscle -in $name.domain.fasta -phyiout $name.phyi -clwout $name.aln");
system ("hmmbuild $name.hmm $name.phyi");
system ("hmmsearch --domtblout $name.serach.txt $name.hmm $pro");
###
### Avoid conflicts between Perl built-in variables and awk built-in variables
my ($a,$b,$c,$d,$e,$f,$g,$h)= ("$1","$2","$3","$4","$5","$6","$7","$8");
my $awksearch = "awk -F '[ ]+' '{if($a !\~ /#/ && $g < 1e-10 ) {print $a,$g}}' $name.serach.txt > $name.id.txt";
`$awksearch`;
#system ("$awksearch");



