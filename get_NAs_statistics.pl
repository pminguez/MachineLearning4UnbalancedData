#!/usr/bin/perl

use strict;

my $file = "/home/pablo/genetica/NeuralNetNeumo/data/repucri16052018_filtro.NAS.t.txt";

open(READ, $file) or die "Cannot open file $file ERROR:$!";

while(<READ>){
    chomp;
    my @cols = split/\t/;
    
    my $mirna = shift @cols;
    
    my $nas = 0;
    
    foreach my $val (@cols){
	if($val eq "NA"){
	    $nas++;
	}
    }
    
    print $mirna ."\t". $nas ."\t". ($nas * 100) / scalar(@cols) ."\n";
}
close(READ);
