#!/usr/bin/perl
# It's not a bug, it's a feature

use strict;
use warnings;

use Getopt::Long;

my $barcode_file; #input sample barcodes to add
my $reads; # reads
my $out = "counts.txt"; #output file for counts
my $fwdPrimer; #forward primer used for amplicon sequencing
my $toprint = 500000; #for updates every $toPrint reads

my $help;

GetOptions('barcodes=s' 	=> \$barcode_file,
		   'reads=s' 		=> \$reads,
		   'out=s'  		=> \$out,
		   'fwdPrimer=s'  	=> \$fwdPrimer,
	       'help',    		=> \$help);

if ($help) { printHelp(); }

unless ($barcode_file and $reads and $fwdPrimer) { printHelp(); }

## first read barcode file
## create array of barcodes

my @barcodes; #to store sequence-name pairs

print "reading barcode file...\n";

open(BARC, $barcode_file) || die "Could not open barcode file: $barcode_file\n";
while (my $line = <BARC>) {

	chomp $line;
	
	my @parts = split(" ", $line); #get fields from line (ID, seq, name)
	
	push(@barcodes, $parts[1]);
}

close(BARC);

#if we didn't get any barcodes from file, die
if ($#barcodes < 0) { die "Couldn't find any barcodes in file $barcode_file\n"; }

#use length of barcode to construct string to add to qualities to account for extra bases
my $new_qual = "J";
$new_qual x= (length $barcodes[0]);

#status update
print "adding barcodes to reads...\n";
my $total = `wc -l < $reads` / 4;
print "total reads: $total\n";

#look through fastq file for reads where barcodes match a specified part of sequence
my $n_lines = 0;
my $query;
open(READS, $reads) || die "Could not open read file: $reads\n";
open(OUT, ">", $out) || die "Could not open output file: $out\n";

while (my $line = <READS>) {

	$n_lines++; #increment line counter
	
	chomp $line;
	
	#is this a line to skip?
	unless ((($n_lines - 2) % 2) == 0) { print OUT "$line\n"; next; } 
	
	#is this a line with qualities?
	unless ((($n_lines - 2) % 4) == 0) { print OUT "$new_qual$line\n"; next; } 
	
	#do we need to reverse-complement?
	unless ($line =~ /$fwdPrimer/i) { $line = reverseComp($line); }
	
	#add barcode to start of line and print
	print OUT "$barcodes[rand @barcodes]$line\n";
	
	#print out progress
	my $n_reads = ($n_lines - 2) / 4;
	if (($n_reads % $toprint) == 0) { print "processed $n_reads reads\n" };

}
close(READS);



#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "\nScript to add barcodes to the start of reads\n\n";
	print "Usage:\n";
	print "\tperl add_sample_barcodes.pl --barcodes <txt> --reads <fastq> --out <txt> --fwdPrimer <seq> --help\n\n";
	print "Arguments:\n";
	print "\t--barcodes:   File with barcodes (txt). First column contains barcode ID, second contains sequence, third contains name\n";
	print "\t--reads:  Fastq file containing reads\n";
	print "\t--out:  Output files for counts of each barcode (txt) [counts.txt]\n";
	print "\t--fwdPrimer:  1-based stop coordinate of barcode [end of read]\n";
	print "\t--help:    Print this help message and quit\n";

	print "\n\t Required arguments: barcodes, reads, fwdPrimer\n";
	exit;
}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}
