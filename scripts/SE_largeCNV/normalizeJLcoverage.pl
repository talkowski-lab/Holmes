#!/usr/bin/perl

#### This script is written by Serkan Erdin in July 2015 to generate log2 ratio profile for 
#### a given sample

use strict;
use warnings;

my $querylist = shift;
my $reflist = shift;
my $sample = shift;

my %refcov=();
my %refvar=();
open(INFO,$reflist) or die;
while(<INFO>){
	chomp;
	if(!/^data/){
		my @array = split /\t/;
		$refcov{$array[0]}{$array[1]} = $array[5];
		$refvar{$array[0]}{$array[1]} = $array[6];
	}
}
close(INFO);

print "Chr\tStart\tEnd\tRefCov\tSampleCov\tLog2\tSignificance\n";

my $selectedcol=();
open(INFO,$querylist) or die;
while(<INFO>){
	chomp;
        if(/^data/){
		my @header = split /\t/;
		for(my $i=0;$i < scalar(@header); $i++){
			if($sample eq $header[$i]){
#				print "$sample\t$header[$i]\n";
				$selectedcol = $i;
			}
		}
	}else{
		my @array = split /\t/;
		if(defined($refcov{$array[0]}{$array[1]})){
			if($refcov{$array[0]}{$array[1]} == 0){
				print "$array[0]\t$array[1]\t$array[2]\t$refcov{$array[0]}{$array[1]}\t$array[$selectedcol]\tNA\tnoref\n";
			}elsif($refcov{$array[0]}{$array[1]} > 0){
				if($array[$selectedcol] == 0){
					 print "$array[0]\t$array[1]\t$array[2]\t$refcov{$array[0]}{$array[1]}\t$array[$selectedcol]\tNA\tnodata\n";
				}elsif($array[$selectedcol] > 0){
					my $log2 = log($array[$selectedcol]/$refcov{$array[0]}{$array[1]})/log(2);
					my $diff = abs($array[$selectedcol] - $refcov{$array[0]}{$array[1]});
					my $significance = $diff/$refvar{$array[0]}{$array[1]};
					print "$array[0]\t$array[1]\t$array[2]\t$refcov{$array[0]}{$array[1]}\t$array[$selectedcol]\t$log2\t$significance\n";
				}
			}
		}
	}
}
close(INFO);
