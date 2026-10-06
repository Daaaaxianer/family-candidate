#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
my $python = $ENV{PYTHON} || 'python3';
exec {$python} $python, "$Bin/refine_domain_seq.py", @ARGV;
die "Cannot execute $python: $!\n";
