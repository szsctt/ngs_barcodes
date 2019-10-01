#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $barcodes;
my $reads;
my $outP;
my $outS;
my $unmatched;
my $mismatches;

GetOptions('barcodes=s'  => \$barcodes,
		   'reads=s'     => \$reads,
		   'outP=s'      => \$outP,
		   'outS=s'      => \$outS,
		   'unmatched=s' => \$unmatched,
		   'mismatches', => \$mismatches);

my %barcodes2id;
my %barcodes2reads;
my @barcodes;

open(BARCODES, $barcodes) || die "Could not open file: $barcodes\n";
while (my $line = <BARCODES>) {
	chomp $line;
	my @parts = split("\t",$line);
		
	my $forward = uc($parts[1]);
	my $reverse = uc($parts[2]);

	$forward =~ tr/Nn/../;
	$reverse =~ tr/Nn/../;

	my $forwardRev = reverseComp($forward);
	my $reverseRev = reverseComp($reverse);

	my $barcode = join("_",($forward,$reverse,$forwardRev,$reverseRev));

	$barcodes2id{uc($barcode)} = $parts[0];
	push(@barcodes, uc($barcode));
}
close BARCODES;

my $header;
my $seq;
my $split;
my $qual;
open (READS,     $reads)        || die "Could not open file: $reads\n";
open (UNMATCHED, ">$unmatched") || die "Could not open file: $unmatched\n";
my $i = 0;

while (my $line = <READS>) {
	$i++;
	chomp $line;
	if ($split) {
		$qual = $line;

		my $matched = 0;
		do {
			foreach my $barcode (@barcodes) {
				my ($fwd, $rev, $fwdR, $revR) = split("_",$barcode);

				my @matches = ($fwd,$rev,$fwdR,$revR);
				if ($mismatches) {
					for (my $i = 0; $i<4; $i++) {
						my $reg = $matches[$i];
						for (my $p = 0; $p < length($reg); $p++) {
							my $temp = $reg;
							substr($temp, $p, 1) = '.';
							push (@matches, $temp);
						}
					}
				}

				my $regex = join("|",@matches);

				if ($seq =~ /^($regex)/i or $seq =~ /($regex)$/i) { # Barcode must be present at start or end of read
					$matched = 1;

					my $output = join('', $outP, $barcodes2id{$barcode}, $outS);
					open (OUTPUT, ">>$output") || die "Could not open file: $output\n";
					print OUTPUT "$header\n$seq\n$split\n$qual\n";
					close OUTPUT;
				}
			}
			unless ($matched == 1) { 
				print UNMATCHED "$header\n$seq\n$split\n$qual\n";
				$matched = 1; 
			}	
		} until($matched == 1);
		
		$header = '';
		$seq    = '';
		$split  = '';
		$qual   = '';
	}

	elsif ($seq)    { $split  = $line; }
	elsif ($header) { $seq    = $line; }
	else            { $header = $line;}
}
close READS;
close UNMATCHED;

exit;

sub reverseComp {
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;
	return($rev);

}