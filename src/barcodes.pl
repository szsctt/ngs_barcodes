# script to count barcodes in NGS data
# each read may contain up to two barcodes at defined locations within the read
# each set of barcodes is specified by a yaml file, and its start coordinate within the read must be given
# output file with counts for each barcode (or combination of barcodes)
# allow for mismatches and variability in barcode location

#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use Getopt::Long;
use YAML::Tiny 'LoadFile';
use Data::Dumper;

# allocate
my $barcodes; # barcodes (yml file)
my $reads; # reads
my $out = "counts.txt"; #output file for counts
my $start; #1-based start coordinate of barcode (if only one)
my $start1; #1-based start coordinate of second barcode
my $start2; #1-based start coordinate of second barcode
my $mismatches; #number of allowed mismatches
my $search; # search this number of bases either side from the expected position of the barcode
my $toprint = 500000; #for updates every $toPrint reads

my $help;

# get inputs from command line
GetOptions('barcodes=s' 		=> \$barcodes,
		   'reads=s' 		=> \$reads,
		   'out=s'  		=> \$out,
		   'start=i'  		=> \$start,
		   'start1=i'  		=> \$start1,
		   'start2=i'  		=> \$start2,
		   'mismatches=i'	=> \$mismatches,
		   'search=i'		=> \$search,
	       'help',    		=> \$help);

if ($help) { printHelp(); }

# check for required inputs
unless ($barcodes and $reads and ($start1 or $start)) { printHelp(); }

if ($start1) { unless ($start2) { printHelp(); } } 

#----------------------------------------------------------------------------------------------------------------------------------------
### Read barcode yaml(s) ###
#----------------------------------------------------------------------------------------------------------------------------------------

## barcode yaml(s) should contain a 'dictionary' of barcode names and sequences for each set of barcodes
##
## example yaml
##
##  barcodes1:
## 		bc1: ATCTAT
## 		bc2: TCTGAA
## 		bc3: TTAGGA
##  barcodes2:
##		bc4: GTACA
##		bc5: GTGAG
## 	

## read yaml files 
#load data from yaml into hashes (one for each set of barcodes)
my $barcs_yaml;
my %barcs1;
my %barcs2;
my $barcodes_name;

# load yaml
$barcs_yaml = LoadFile($barcodes);

#check number of sets of barcodes in yaml
if (scalar keys %{$barcs_yaml} > 2) { die "found too many sets of barcodes in yaml"; }
if (scalar keys %{$barcs_yaml} == 0) { die "didn't find any barcodes in yaml"; }

# get first set of barcodes
$barcodes_name = (keys %{$barcs_yaml})[0];
print "loading barcode set 1: $barcodes_name\n";
%barcs1 = %{ $barcs_yaml->{$barcodes_name} };

# if second set of barcodes, get these too
if (scalar keys %{$barcs_yaml} > 1) { $barcodes_name = (keys %{$barcs_yaml})[1]; }
print "loading barcode set 2: $barcodes_name\n";
%barcs2 = %{ $barcs_yaml->{$barcodes_name} };

# check that number of starts specified corresponds with number of sets of barcodes in input file
if ($barcs_yaml->{barcodes2}) { 
	unless ($start1 and $start2) { 
		print "Specify --start1 and --start2 for two sets of barcodes\n"; 
	} 	
} 
if ($barcs_yaml->{barcodes1}) { 
	unless ($start) { 
		die "Specify --start for one set of barcodes\n "; 
	} 
}

# construct data structure to store counts of barcodes
# if only one set of barcodes, this is a hash with barcode sequences as keys and count as values
# if two sets of barcodes, this is a 



#status update
print "finding barcodes in reads...\n";
my $total = `wc -l < $reads` / 4;
print "$total reads to check\n";

#look through fastq file for reads where barcodes match a specified part of sequence
my $n_lines = 0;
my $query;
open(READS, $reads) || die "Could not open read file: $reads\n";

while (my $line = <READS>) {

	$n_lines++; #increment line counter
	
	#is this a line with a read sequence?
	unless ((($n_lines - 2) % 4) == 0) { next; } 
	
	#do we need to reverse-complement?
	#unless ($line =~ /$fwdPrimer/i) { $line = reverseComp($line); }
	
	#get sub-string from start to stop
	#$query = substr($line, $start, $length);  #get barcode from read
	
	#increment barcode counter
	if (exists $counts{$query}) { $counts{$query} += 1; } #if matched (assumes one unique match)
	else { $counts{"none"} += 1; } #if not matched
	
	#print out progress
	my $n_reads = ($n_lines - 2) / 4;
	if (($n_reads % $toprint) == 0) { print "seached $n_reads reads\n" };

}
close(READS);

#print out results
print "\n counts for barcodes: \n";

#print out counts to terminal and file
open(OUT, ">", $out) || die "Could not open output file: $out\n";
print OUT "name\tseq\tcount\n";
foreach my $barcode (keys %counts) { 

	my $count = $counts{$barcode};
	my $name = $lookup{$barcode}; 
	
	#deal with "none"
	if ($barcode =~ /none/) {$name = "umatched"; }
	
	print "$name ($barcode) : $count \n"; 
	print OUT "${name}\t${barcode}\t${count}\n";

}
close(OUT);

exit;

#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "\nScript for counting barcodes in reads\n\n";
	
	print "Usage:\n";
	print "\tperl barcodes.pl --barcodes <txt> --reads <fastq> --out <txt> --fwdPrimer <seq> --start <num> --length <num> --help\n\n";
	print "Arguments:\n";
	print "\t--barcodes:   Yaml file containing one or two sets of barcodes\n";
	print "\t--reads:  Fastq file containing reads\n";
	print "\t--out:  Output files for counts of each barcode (txt) [counts.txt]\n";
	print "\t--start:  1-based start coordinate of barcode [1]\n";
	print "\t--length:  length of barcode [length of read]\n";
	print "\t--help:    Print this help message and quit\n";

	print "\n\t Required arguments: barcodes, reads, fwdPrimer, start, length\n";
	exit;
}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

