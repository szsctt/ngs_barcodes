# script to count barcodes in NGS data
# each read contains a barcode at a location defined by start
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
my $prefix; # prefix for output file names [barcode names]
my $outpath = "."; #output folder for counts and read files
my $start; #1-based start coordinate of barcode
my $mismatches; #number of allowed mismatches
my $toprint = 500000; #for updates every $toPrint reads

my $help;

# get inputs from command line
GetOptions('barcodes=s' 	=> \$barcodes,
		   'reads=s' 		=> \$reads,
		   'outpath=s'  	=> \$outpath,
		   'prefix=s'  		=> \$prefix,
		   'start=i'  		=> \$start,
		   'mismatches=i'	=> \$mismatches,
	       'help',    		=> \$help);

if ($help) { printHelp(); }

# check for required inputs
unless ($barcodes and $reads and $start and $outpath) { printHelp(); }

#check outpath ends wihth "/"
unless ($outpath =~ /$\//) { $outpath = $outpath."/"; }


#----------------------------------------------------------------------------------------------------------------------------------------
### Read barcode yaml ###
#----------------------------------------------------------------------------------------------------------------------------------------

## barcode yaml should contain a hash of barcode sets with names and sequences for each set of barcodes
##
## example yaml
##
## 	bc1: ATCTAT
## 	bc2: TCTGAA
## 	bc3: TTAGGA


## read yaml files 

my $barcs_yaml;
my %name2seq;
my $barcodes_name;
my $barcs_len;

# load yaml
$barcs_yaml = LoadFile($barcodes);

#check number of sets of barcodes in yaml
if (scalar keys %{$barcs_yaml} > 1) { die "Error: found too many barcodes in yaml\n"; }
if (scalar keys %{$barcs_yaml} == 0) { die "Error: didn't find any barcodes in yaml\n"; }

# extract useful parts of yaml object
$barcodes_name = (keys %{$barcs_yaml})[0];
print "loading barcodes: $barcodes_name\n";
%name2seq = %{ $barcs_yaml->{$barcodes_name} }; # hash of barcode names and sequences
$barcs_len = length((values %name2seq)[0]); # barcode length


#----------------------------------------------------------------------------------------------------------------------------------------
### Count barcodes ###
#----------------------------------------------------------------------------------------------------------------------------------------

# construst hash with barcode names as keys and count as values to store results
my %counts;
foreach my $name (keys %name2seq) { $counts{$name} = 0; } # initialise all counts to 0
$counts{"none"} = 0; # to count reads with no matches
$counts{"ambiguous"} = 0; # to count reads with matches to more than one barcode or where forward and reverse reads match different barcodes

# construct hash to store filenames of output read files
my %outnames;
if ($prefix) {
	# construct filenames with prefix
	foreach my $name (keys %counts) { $outnames{$name} = join("", $outpath, $prefix, "_", $name, ".fastq"); } 
}
else {
	# construct filenames without prefix
	foreach my $name (keys %counts) { $outnames{$name} = join("", $outpath, $name, ".fastq"); } 
}

# search strategy for barcodes
my %barcodes;

## if we are allowing mismatches, construct regexes to use to find barcodes in reads
## create hash with regex as key and name of barcode as value
if ($mismatches) { 

	#first check that we aren't asking for more mismatches than the length of the barcode
	if ($mismatches >= $barcs_len) { die "Number of mismatches cannot be greater than or equal to barcode length\n"; }
	
	#for each barcode, construct regex
	while (my ($name, $seq) = each (%name2seq)) {
		#construct array which will be collapsed into regex
		my @regex;
		 push(@regex, $seq, reverseComp($seq));
		
		#add '.' to every position in strings in regex array $mismatch number of times
		my $i = 0;
		while ($i < $mismatches) { @regex = addAny(@regex); $i++; }
		
		#collapse array to create regex
		$barcodes{(join("|", @regex))} = $name;
		
	}
}
# if we aren't allowing any mismatches, it's much faster to use hash lookup to see if a sequence is one of the barcodes
## in this case, construct a hash with the barcode sequences as keys and the name of the barcode as values
## to allow for reverse reads, also add the reverse complements of the barcodes as keys
else {
	%barcodes = reverse %name2seq; 
	#also need to add the reverse complements of each sequence, in case reads are reversed
	while (my ($name, $seq) = each (%name2seq)) { $barcodes{reverseComp($seq)} = $name; }
}


#status update
print "finding barcodes in reads...\n";
my $total = `wc -l < $reads` / 4;
print "$total reads to check\n";

#look through fastq file for reads where barcodes match a specified part of sequence
my $n_lines = 0;
my @lines;
my $query_fwd;
my $query_rev;

open(READS, $reads) || die "Could not open read file: $reads\n";
while (my $line = <READS>) {
	
	#get fastq four lines at a time
	push @lines, $line;
    if (scalar(@lines) == 4) {

		#get read sequence from @lines
		my $seq = $lines[1];
		chomp $seq;
	
		#get sub-string from start to stop for forward and reverse case
		my $query_fwd = substr($seq, $start, $barcs_len);
		my $query_rev = substr($seq, (length($line) - $start), $barcs_len);
	
		my @matches = ();
	
		#check for barcode (mismatches allowed)
		if ($mismatches) {
	
			#find regex matches
			@matches = grep { $query_fwd =~ /$_/ } keys %barcodes;
			push @matches, grep { $query_rev =~ /$_/ } keys %barcodes;

			#check number of matches, write to outfile and increment counter
			checkMatches(\@matches, \%counts, \%outnames, \@lines, \%barcodes);
		
		}
	
		#check for barcode (no mismatches allowed)
		else {
		
			#check if query is key in %barcodes
			if (exists $barcodes{$query_fwd}) { push @matches, $query_fwd; }
			if (exists $barcodes{$query_rev}) { push @matches, $query_rev; }
		
			#check number of matches, write to outfile and increment counter
			checkMatches(\@matches, \%counts, \%outnames, \@lines, \%barcodes);
		}
	
		#print out progress
		$n_lines++;
		if (($n_lines % $toprint) == 0) { print "seached $n_lines reads\n" };
		
		#clear line buffer
		@lines = ();
		
    }

}
close(READS);

#print out results
print "\n counts for barcodes: \n";

#print out counts to terminal and file
my $countfile = join("", $outpath, $prefix, "_counts.txt");
open(OUT, ">", $countfile) || die "Could not open output file: $countfile\n";
print OUT "name\tseq\tcount\n";
foreach my $name (keys %counts) { 

	my $count = $counts{$name};
	my $seq= $name2seq{$name}; 
	
	#deal with "none"
	if ($name =~ /none/) {$seq = "NA"; }
	if ($name =~ /ambiguous/) {$seq = "NA"; }
	
	print "$name ($seq) : $count \n"; 
	print OUT "${name}\t${name}\t${count}\n";

}
close(OUT);

exit;

#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "\nScript for counting barcodes in reads\n\n";
	
	print "Usage:\n";
	print "\tperl barcodes.pl --barcodes <yml> --reads <fastq> --outpath <path> --start <num> --help\n\n";
	print "Arguments:\n";
	print "\t--barcodes:\tYaml file containing one or two sets of barcodes\n";
	print "\t--reads:\tFastq file containing reads\n";
	print "\t--outpath:\tOutput folder for counts of each barcode [.]\n";
	print "\t--prefix:\tPrefix for output file names [barcode names]\n";
	print "\t--start:\t1-based start coordinate of barcode\n";
	print "\t--mismatches:\tAllow this many mismatches in barcodes [0]\n";
	print "\t--help:\t\tPrint this help message and quit\n";

	print "\n\tRequired arguments: barcodes, reads, start\n";
	exit;
}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);
}

sub addAny {
### take an array of strings and add new elements in which '.' replaces every possible position in each string
### eg if input is ['ATC'] function returns ['.TC', 'A.C' and 'AT.'] 

	my @array = @_;
	
	my @new_array;
	
	foreach my $element (@array) { #iterate over elements of array
		for (my $j=0; $j<(length $element); $j++) { #iterate over letters in string
			my $temp = $element;
			substr($temp, $j, 1) = '.'; #replace ith letter with '.'
			push(@new_array, $temp);
		}
	}
	
	return @new_array;
	
}

sub incrementAndWrite {
#increment counter in hash and write lines in array to outfile

	my ($countRef,  $outfileRef, $key, $arrayRef) = @_;
	
	#increment counter
	$countRef->{$key} += 1;
	
	#open file for writing or appending
	if ($countRef->{$key} == 0) {
		open(OUTFILE, '>', $outfileRef->{$key}) || die "Could not open file $outfileRef->{$key})";
	}
	else {
		open(OUTFILE, '>>', $outfileRef->{$key}) || die "Could not open file $outfileRef->{$key})";
	}
	#write lines
	print OUTFILE @$arrayRef;
	
	#close file
	close(OUTFILE);

}

sub checkMatches {
#check the number of matches, increment and write to the appropriate counters and files

	my ($matchArrayRef, $countHashRef, $outfileHashRef, $writeArrayRef, $barcodesHashRef) = @_;
	
	#check number of matches, increment counter and write to output file
	#one match
	if ( (scalar(@$matchArrayRef)) == 1) {
		 incrementAndWrite($countHashRef, $outfileHashRef, $barcodesHashRef->{($$matchArrayRef[0])}, $writeArrayRef);
	}
	#more than one match
	elsif ((scalar(@$matchArrayRef)) > 1) { 
		incrementAndWrite($countHashRef, $outfileHashRef, "ambiguous", $writeArrayRef);
	}
	#no matches
	else { 
		incrementAndWrite($countHashRef, $outfileHashRef, "none", $writeArrayRef);
	}


}