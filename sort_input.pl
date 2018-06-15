#!/usr/bin/perl

use strict;

my $file = "/home/pablo/genetica/NeuralNetNeumo/data/repucri16052018_filtro.txt";
my @mortality1;
my @mortality0;
my $header;

open(READ, $file) or die "Cannot open file ERROR:$!";

while(<READ>){
    chomp;
    my $line = $_;
    my @bits = split/\t/;
    shift(@bits);
    if($bits[12] eq "MORTALITY"){
	$header = $line;
    }   
    elsif($bits[12] == 1){

	my $l = "";
	foreach my $b (@bits){
	    $l .= $b ."\t";
	}
	chop($l);
	
	push @mortality1, $l;
    }
    elsif($bits[12] ne ""){

	my $l = "";
	foreach my $b (@bits){
	    $l .= $b ."\t";
	}
	chop($l);
	
	push @mortality0, $l;
    }
}
close(READ);

my $outfile = "/home/pablo/genetica/NeuralNetNeumo/data/repucri16052018_filtro_sorted.txt";
open(OUT, ">$outfile") or die "Cannot open file $outfile ERROR:$!";

print OUT $header ."\n";

foreach my $line (@mortality1){
    print OUT $line ."\n";
}

foreach my $line (@mortality0){
    print OUT $line ."\n";
}
close(OUT);
