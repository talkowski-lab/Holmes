#!/usr/bin/perl

use strict;
use warnings;

my $list = shift;

my %countzero=();
my %countall =();	
open(INFO,$list) or die;
while(<INFO>){
	chomp;
	if(!/^data/){
	my @array = split /\t/;
	$countall{$array[0]}++;
	if($array[3] == 0){
		$countzero{$array[0]}++;
	}
	}
}
close(INFO);

foreach my $key (keys %countall){
	if(!defined($countzero{$key})){
		$countzero{$key}=0;
	}
	print "$key\t$countall{$key}\t$countzero{$key}\n";
}	
