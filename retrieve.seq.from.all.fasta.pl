#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
my $python = $ENV{PYTHON} || 'python3';
exec {$python} $python, "$Bin/retrieve.seq.from.all.fasta.py", @ARGV;
die "Cannot execute $python: $!\n";
