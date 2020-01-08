# script to count barcodes in NGS data
# each read contains a barcode at a location defined by start
# each set of barcodes is specified by a yaml file, and its start coordinate within the read must be given
# output file with counts for each barcode by set and combination of barcodes
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
my $extend_search = 0; #extend search this many bases either side of the expected position of the barcode
my $prefix; # prefix for output file names [""]
my $outpath = "."; #output folder for counts and read files
my $mismatches = 0; #number of allowed mismatches
my $fwdPrimer; #forward primer to check if read is forward or reverse
my $n_to_process;
my $toprint = 10000; #for updates every $toPrint reads
my $help;


# get inputs from command line
GetOptions('barcodes=s' 		=> \$barcodes,
		   'reads=s' 			=> \$reads,
		   'outpath=s'  		=> \$outpath,
		   'prefix=s'  			=> \$prefix,
		   'extend_search=i'  	=> \$extend_search,
		   'mismatches=i'		=> \$mismatches,
		   'fwdPrimer=s'		=> \$fwdPrimer,	   
		   'n_reads=i'			=> \$n_to_process,	  
	       'help',    			=> \$help);

if ($help) { printHelp(); }

# check for required inputs
unless ($barcodes and $reads and $outpath) { printHelp(); }

#check outpath ends wihth "/"
unless ($outpath =~ /\/$/) { $outpath = $outpath."/"; }

#if we want to check for mismatches or extend the search use a regex 
#otherwise use hash lookup
my $regex;
if ($mismatches || $extend_search) { $regex = 1; }

#keep these while debugging
my %barcodes;
my %outnames;
#delete after rewrite
my $start;
my %name2seq;
my $barcs_len;
my $printFastq;

#----------------------------------------------------------------------------------------------------------------------------------------
### Read barcode yaml ###
#----------------------------------------------------------------------------------------------------------------------------------------

## barcode yaml should contain a hash of barcode sets with names and sequences for each set of barcodes
##
## example yaml
##
##	barcodes1:
##		start: 72
##		barcodes:
## 			bc1: ATCTAT
## 			bc2: TCTGAA
## 			bc3: TTAGGA
##	barcodes2:
##		start: 0
##		barcodes:
## 			bc3: ACCTCT
## 			bc4: TCTGCC
## 			bc5: GTAGGG


## read yaml files 

# load yaml
my $barcs_yaml;
$barcs_yaml = LoadFile($barcodes);

my @allBarcs; #need this for constructing barcode combination counter
my @setOrder; #keep track of order because lookup is random

#add extra information to yaml needed for search
foreach my $set (keys %{ $barcs_yaml}) {

	#because hash lookup is random, need to keep track of order for alter
	push @setOrder, $set;

	#get barcode length
	$barcs_yaml->{$set}->{"length"} = length((values %{ $barcs_yaml->{$set}->{"barcodes"} })[0]);
	
	#get array of barcodes to create counter later on
	my @barcsArray = keys %{$barcs_yaml->{$set}->{"barcodes"}};
	push @barcsArray, qw(none ambiguous);
	push @allBarcs, \@barcsArray;
	
	#construct counter to keep track of counts of barcodes
	$barcs_yaml->{$set}->{"counts"} = {};
	foreach my $barc (keys %{$barcs_yaml->{$set}->{"barcodes"}} ) {
		$barcs_yaml->{$set}->{"counts"}->{$barc} = 0;
	}
	#add 'none' and 'ambiguous'
	$barcs_yaml->{$set}->{"counts"}->{"none"} = 0;
	$barcs_yaml->{$set}->{"counts"}->{"ambiguous"} = 0;
	
	#construct regexes if we are doing regex search
	if ($regex) {
	
		foreach my $barc (keys %{$barcs_yaml->{$set}->{"barcodes"}} ) {
			my @regex;
			
			#if we know the forward primer, only need to search for barcode
			#otherwise need to search for reverse complement as well
			my $seq = $barcs_yaml->{$set}->{"barcodes"}->{$barc};
			if ($fwdPrimer) { push(@regex, $seq);}
			else { push(@regex, $seq, reverseComp($seq));}
			
			#if we are allowing mismatches
			if ($mismatches) {
			
				#add '.' to every position in strings in regex array $mismatch number of times
				my $i = 0;
				while ($i < $mismatches) { @regex = addAny(@regex); $i++; }
			}
			
			#collapse array to create regex
			$barcs_yaml->{$set}->{"search"}->{(join("|", @regex))} = $barc;
		}
	}
	#for 'search' hash (for hash lookup), we need barcodes as keys and names as values (ie reverse of "barcodes" has)
	else {
		foreach my $barc (keys %{$barcs_yaml->{$set}->{"barcodes"}} ) {
			#add sequence as key to hash with barcode name as value
			my $seq = $barcs_yaml->{$set}{"barcodes"}{$barc};
			$barcs_yaml->{$set}{"search"}{$seq} = $barc;
			#if we don't know the forward primer, we also need to add reverse complements of barcodes as keys
			unless ($fwdPrimer) {
				$barcs_yaml->{$set}{"search"}{reverseComp($seq)} = $barc;
			}
		}
	}	
}

#create data structure to store counts of combinations of barcodes
#this is a nested hash, where the number of levels corresponds to the number of sets of barcodes
#each path through the nested hash corresponds to a unique combination of barcodes 
my $counts = {};

#first get all possible combinations of barcodes
my $combinations = cartProduct(\@allBarcs);

#then create nested hash with as many levels as sets of barcodes
#leaf values are initialised to zero
my %counts;
foreach my $combination (@$combinations) {
	setValue(\%counts, $combination, 0);
}
	

#----------------------------------------------------------------------------------------------------------------------------------------
### Count barcodes ###
#----------------------------------------------------------------------------------------------------------------------------------------

#status update
my $total = `wc -l < $reads`;
$total /= 4;
print "finding barcodes in $total reads\n";

#first de-replicate fastq file
#for each read, add to hash with sequence of read as key and count of read as value
my %reads;
my @lines;
my $i = 0;
my $dropped = 0;

open(READS, $reads) || die "Could not open read file: $reads\n";
while (my $line = <READS>) {
	
	$i++;
	
	if ($i % 4 == 2) {

		chomp $line;
		
		#if we know forward primer, check if read is forward or reverse
		if ($fwdPrimer) {
			if ($line !~ /$fwdPrimer/i) { 
				$line = reverseComp($line);
				#if still doesn't match, drop it
				unless ($line =~ /$fwdPrimer/i) { 
					$dropped += 1;
					next;
				 }
			}
		}
		
		#add sequence to hash and increment counter
		if (exists $reads{$line}) { $reads{$line} += 1; }
		else { $reads{$line} = 0; }
    }
    
    if ($n_to_process) {
		if ($i/4 > $n_to_process) { last; }
	}
}	


close(READS);

print "dropped $dropped reads that don't match the forward primer after reversing\n";

my $n_uniq = scalar keys %reads;
print "after de-duplicaton, $n_uniq unique reads to process\n";

#look through %reads and count reads where barcodes match a specified part of sequence
my $n_lines = 0;
my @matches;
my $line_len;
$i = 0;

foreach my $seq (keys %reads) {
		
	#to store references to arrays containing names of matched barcodes for each set
	@matches = ();
	
	#get length of line
	$line_len = length($seq);
	
	#get queries from each set and do search
	foreach my $set (@setOrder) {
	
		if ($n_lines == 0) {
		#check that length of line is more than $start + barcode length
			unless (($barcs_yaml->{$set}->{"start"} + $barcs_yaml->{$set}->{"length"}) <= $line_len) {
				print "Error: line is shorter than expected end of barcode\n";
				die "Is the start position correct and the read file unzipped?\n"
			}
		}
	
		#get start and stop of part of read to search
		my $start = $barcs_yaml->{$set}->{"start"} - $extend_search;
		if ($start < 0) { $start = 0; }
		if ($start > $line_len) { $start = $line_len; }
		my $length = $barcs_yaml->{$set}->{"length"} + ($extend_search * 2);
		if ($line_len < ($length + $start)) { $length = $line_len - $start; }
		
		#get relevant part of read
		my @setQueries;
		push @setQueries, uc(substr($seq, $start, $length));
		
		#do the same if we need to check the reverse reads
		unless ($fwdPrimer) {
			
			push @setQueries, reverseComp(uc(substr($seq, -1*($line_len - $start), $length)));
		}
		
		#if string will be zero length, don't bother
		if ($start == $line_len) { 
			push @matches, [];
			next; 
		} #if string will be zero length, don't bother
		
		
		#get matches
		my @setMatches;
		if ($regex) {
			foreach my $query (@setQueries) {
				#do regex matching
				push @setMatches, grep { $query =~ /$_/i } keys %{$barcs_yaml->{$set}->{"search"}};
			}
		}
		else {
			#check for barcode by hash lookup
			foreach my $query (@setQueries) {
				if (exists $barcs_yaml->{$set}->{"search"}->{$query}) {
					push @setMatches, $query;
				}
			}
		}
		#get unique names from query matches
		foreach my $match (@setMatches) {
			$match = $barcs_yaml->{$set}->{"search"}->{$match};
		}
		my %seen = ();
		my @uniqSetMatches= grep { ! $seen{$_} ++ } @setMatches;
		
		#push set matches to @matches
		push @matches, \@uniqSetMatches;	
	}
	
	#check matches and increment counter
	checkAllMatches(\@matches, $barcs_yaml, \%counts, \@setOrder, $reads{$seq});
	
	$i++;
	if ($i % $toprint == 0) { print "processed $i reads\n"; }	
}




#----------------------------------------------------------------------------------------------------------------------------------------
### Write counts to file  ###
#----------------------------------------------------------------------------------------------------------------------------------------

#print out results
@allBarcs = (); #need this for writing combination file
my $countfile;

#print out counts to terminal and file
foreach my $set (@setOrder) {

	#push set barcode names to @allbarcs for combination file
	my @barcsArray = keys %{$barcs_yaml->{$set}->{"barcodes"}};
	push @barcsArray, qw(none ambiguous);
	push @allBarcs, \@barcsArray;
	
	print "\ncounts for barcodes: $set \n";
	
	#make file name
	if ($prefix) { $countfile = "${outpath}${prefix}_${set}_counts.txt"; }
	else { $countfile = "${outpath}${set}_counts.txt"; }
	
	#open output file and write header
	open(OUT, ">", $countfile) || die "Could not open output file: $countfile\n";
	print OUT "name\tseq\tcount\n";
	
	#get name, sequence and count and write to file
	foreach my $name (@barcsArray) {
		my $seq = $barcs_yaml->{$set}->{"barcodes"}->{$name};
		my $count = $barcs_yaml->{$set}->{"counts"}->{$name};
		
		#deal with "none"
		if ($name =~ /none/) {$seq = "NA"; }
		if ($name =~ /ambiguous/) {$seq = "NA"; }
		
		#print to terminal and file
		print "\t$name ($seq) : $count\n"; 
		print OUT "${name}\t${seq}\t${count}\n";
	}
	close(OUT);
}

#get all possible combinations of barcodes which is the complete set of keys for %counts
@$combinations = ();
$combinations = cartProduct(\@allBarcs);

#make file name
#make file name
if ($prefix) { $countfile = "${outpath}${prefix}_counts.txt"; }
else { $countfile = "${outpath}counts.txt"; }

#open file for combined counts
open(OUT, ">", $countfile) || die "Could not open output file: $countfile\n";
my $header = ("name\t" x scalar(@setOrder))."count\n";
print OUT $header;

#then get every value from nested hash and write to file
#leaf values are initialised to zero
foreach my $combination (@$combinations) {

	#join names
	my $toprint = join("\t", @$combination);
	
	#write value
	my $count = getValue(\%counts, $combination);

	#print to outfile
	print OUT "${toprint}\t${count}\n";
}
		
close(OUT);
		
exit;

#----------------------------------------------------------------------------------------------------------------------------------------
### Subroutines ###
#----------------------------------------------------------------------------------------------------------------------------------------

sub printHelp {
	print "\nScript for counting barcodes in reads\n\n";
	
	print "Usage:\n";
	print "\tperl barcodes.pl --barcodes <yml> --reads <fastq> --outpath <path> --help\n\n";
	print "Arguments:\n";
	print "\t--barcodes:\tYaml file containing one or two sets of barcodes\n";
	print "\t--reads:\tFastq file containing reads\n";
	print "\t--outpath:\tOutput folder for counts of each barcode [.]\n";
	print "\t--prefix:\tPrefix for output file names [barcode names]\n";
	print "\t--start:\t1-based start coordinate of barcode\n";
	print "\t--mismatches:\tAllow this many mismatches in barcodes [0]\n";
	print "\t--extend_search:\tSearch this many bases each side of expected barcode positiohn [0]\n";
	print "\t--fwdPrimer:\tOptionally provide sequence of foward primer common to all reads to improve search \n";
	print "\t--n_reads:\tLimit analysis to the first n reads\n";
	print "\t--help:\t\tPrint this help message and quit\n";

	print "\n\tRequired arguments: barcodes, reads\n";
	exit;
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

sub cartProduct {
#take an array of array references, and generate the cartesian product of all the referenced arrays
#note that this destroys the input array

	my $inArrayRef = pop @_;
	
	#create a reference to an array with one empty array reference for output
	my $outArrayRef = [];
	push @$outArrayRef, [];
	
	while (@{$inArrayRef}) {
	
		#get input array for adding
		my $currInArrayRef = pop @{$inArrayRef};
		
		#for each array in @$outArrayRef, remove it and add new arrays 
		#each consisting of that array with every element from @$currInArrayRef added
		my $i = 0;
		my $outArrayLen = scalar @{$outArrayRef};
		while ($i < $outArrayLen) {
			my $currOutArrayRef = pop @$outArrayRef;
			
			foreach my $new (@$currInArrayRef) {
				my $newArrayRef = [];
				#add new element to old out array
				push @$newArrayRef, @$currOutArrayRef, $new;
				
				#add old out array to out array of arrays
				unshift @$outArrayRef, $newArrayRef;
				
			}
			$i++;
		}
	}
	
	return $outArrayRef;
}

sub checkAllMatches {
#check the number of matches, increment and write to the appropriate counters and files

	my ($matches, $barcs_yaml, $counts, $setOrder, $count) = @_;
	
	#array containing unique entry to increment in combined counter
	my @combinedMatches = ();
	
	#first increment individual counters in $barcs_yaml
	foreach my $i (0..$#{$setOrder} ) {
		my $set = $setOrder->[$i];
		my $matchRef = $matches->[$i];
		incrementSetMatches($matchRef, $barcs_yaml, $set, \@combinedMatches, $count);
	}
	
	@combinedMatches = reverse @combinedMatches;
	
	#then increment combined counter
	incrementValue($counts, \@combinedMatches, $count);

}

sub getValue {
#create or update one leaf in recursive hash
#pass in reference to hash, array with each element corresponding to a key
#and value to set

	my ($ref, $arrayRef) = @_;
	
	#tail will be leaf key
	my $tail = pop @{$arrayRef};
	#cursor is current location in hash
	my $cursor = $ref;
	
	#iterate over levels
	foreach my $level (@{$arrayRef}) {

		#traverse down
		$cursor = $cursor->{$level};
	}
	
	#return value
	return $cursor->{$tail};

}

sub incrementSetMatches {
#increment counter in hash and write lines in array to outfile

	my ($matchArrayRef, $barcs_yaml, $set, $combinedMatches, $count) = @_;
	
	#check number of matches, increment counter and write to output file
	#one match
	if ( (scalar(@$matchArrayRef)) == 1) {

		#increment counter
		$barcs_yaml->{$set}->{"counts"}->{@{$matchArrayRef}[0]} += $count;
		#push to combinedMatches
		push @$combinedMatches, @{$matchArrayRef}[0];
	}
	#more than one match
	elsif ((scalar(@$matchArrayRef)) > 1) { 
		#increment counter
		$barcs_yaml->{$set}->{"counts"}->{"ambiguous"} += $count;
		#push to combinedMatches
		push @$combinedMatches, "ambiguous";
	}
	#no matches
	else { 
		#increment counter
		$barcs_yaml->{$set}->{"counts"}->{"none"} += $count;
		#push to combinedMatches
		push @$combinedMatches, "none";
	}

}

sub incrementValue {
#create or update one leaf in recursive hash
#pass in reference to hash, array with each element corresponding to a key
#and value to set

	my ($ref, $arrayRef, $count) = @_;
	
	#if count not specified, add 1
	unless ($count) { $count = 1; }
	
	#tail will be leaf key
	my $tail = pop @{$arrayRef};
	#cursor is current location in hash
	my $cursor = $ref;
	
	#iterate over levels
	foreach my $level (@{$arrayRef}) {

		#traverse down
		$cursor = $cursor->{$level};
	}
	
	#set value
	$cursor->{$tail} += $count;

}

sub setValue {
#create or update one leaf in recursive hash
#pass in reference to hash, array with each element corresponding to a key
#and value to set

	my ($ref, $arrayRef, $value) = @_;
	
	my $tail = pop @{$arrayRef};
	my $cursor = $ref;
	
	#iterate over levels
	foreach my $level (@{$arrayRef}) {
		#create a new hash if there isn't one
		$cursor->{$level} ||= {};
		#traverse down
		$cursor = $cursor->{$level};
	}
	
	#set value
	$cursor->{$tail} = $value;

}

sub reverseComp {
### return reverse compliment of sequence
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/NATGCatgc/NTACGtacg/;
	return($rev);
}
