#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use Getopt::Long;

my $barcodes; #input barcodes
my $reads; # reads
my $out = "counts.txt"; #output file for counts
my $start; #1-based start coordinate of barcode
my $length; #1-based stop coordinate of barcode
my $fwdPrimer; #forward primer used for amplicon sequencing
my $toprint = 500000; #for updates every $toPrint reads

my $help;

GetOptions('barcodes=s' 	=> \$barcodes,
		   'reads=s' 		=> \$reads,
		   'out=s'  		=> \$out,
		   'start=i'  		=> \$start,
		   'length=i' 		=> \$length,
		   'fwdPrimer=s'  	=> \$fwdPrimer,
	       'help',    		=> \$help);

if ($help) { printHelp(); }

unless ($barcodes and $reads and $fwdPrimer and $start and $length) { printHelp(); }

## first read barcode file
## create an array of barcode sequences
## also two hashes: one with sequences as keys and names as values, and one with name as key and count as value

my %counts; #to store sequence-count pairs
my %lookup; #to store sequence-name pairs

#add "none" to %counts to keep track of reads where there was no match
$counts{"none"} = 0;

print "reading barcode file...\n";

open(BARC, $barcodes) || die "Could not open barcode file: $barcodes\n";
while (my $line = <BARC>) {

	my @parts = split(" ", $line); #get fields from line (ID, seq, name)
	
	if ($parts[2]) { #if barcode has a name
	
		$counts{$parts[1]} = 0; #start count at 0 (seq : count)
		$lookup{$parts[1]} = $parts[2]; #add to lookup table (seq : name)
	}
}
close(BARC);

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
	unless ($line =~ /$fwdPrimer/i) { $line = reverseComp($line); }
	
	#get sub-string from start to stop
	$query = substr($line, $start, $length);  #get barcode from read
	
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
	print "\t--barcodes:   File with barcodes (txt). First column contains barcode ID, second contains sequence, third contains name\n";
	print "\t--reads:  Fastq file containing reads\n";
	print "\t--out:  Output files for counts of each barcode (txt) [counts.txt]\n";
	print "\t--fwdPrimer:  1-based stop coordinate of barcode [end of read]\n";
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

