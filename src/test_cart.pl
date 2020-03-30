use strict;
use warnings;

use Getopt::Long;
use Data::Dumper;


my @one = qw(one two three);
my @two = qw(four five six);
my @three = qw(seven eight nine);
my @comb;
$comb[0] = \@one;
$comb[1] = \@two;
$comb[2] = \@three;

#print Dumper(\@comb);
my @product;
cartProduct(\@comb, \@product);


sub cartProduct {
# take an array of array references, 
#and generate the cartesian product of all the referenced arrays (in an output array references)
# note that this destroys the input array

	my ($inArrayRef, $outArrayRef) = @_;
	my $len = $#$inArrayRef;
	
	#if $inArrayRef is empty, return
	if ($len == 0) {
		return;
	}
	
	my @newOutArray = ();
	foreach my $array (@{$inArrayRef}) {
		foreach my $element ( @{${$outArrayRef}[0]} )
		 {
			my @tmp = push @{$array}, $element;
			push @newOutArray, \@tmp;
		}
	}
	
	cartProduct(\{@{${$inArrayRef}[1..$len]}}, \@newOutArray);
}