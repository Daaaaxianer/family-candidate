#! /usr/bin/perl
use strict;
use Bio::Seq;
use Bio::SeqIO;
#use Math::BigFloat;

die "perl $0 <hmmout> <protein.fasta> <outfile> <e-value>" unless (@ARGV == 4);

my $hmmout = Bio::SeqIO -> new (-file => "$ARGV[1]", 
	                            -format => "Fasta");

my $domfasta = Bio::SeqIO -> new (-file => ">$ARGV[2]", 
	                            -format => "Fasta");

my %keep=();
open IN, "$ARGV[0]" or die "$!";
while (<IN>){
	chomp;
	next if /^#/;

	my @arr = split(/\s+/,$_);
	next if $arr[6] > $ARGV[3];
	my @brr = ($arr[17],$arr[18]);
	my $target = $arr[0];
	if(!exists $keep{$target}){
		$keep{$target} = \@brr;
	}
}
close (IN);

while (my $seq = $hmmout -> next_seq()){
	my ($id,$sequence,$desc) = ($seq -> id,$seq -> seq,$seq -> desc);
	if(exists $keep{$id}){
		print "$id $keep{$id}[0] $keep{$id}[1]\n";
		### use trunc
        my $newseqobj = $seq -> trunc($keep{$id}[0],$keep{$id}[1]);
        ### use subseq
        # my $newseq = $seq -> subseq($keep{$id}[0],$keep{$id}[1]);

        # my $newseqobj = Bio::Seq ->new(-seq => $newseq,
        # 	                          -desc => 'subseq',
        # 	                          -id => $id);
		$domfasta -> write_seq($newseqobj);
	}
}
$hmmout -> close();
$domfasta -> close();
