#!/usr/bin/perl

use strict;

my $datadir = "/home/pablo/genetica/miRNAs_Vito/data/second_round/";

my @headers = ("miR-144-3p","miR-144-5p","miR-1246","miR-486-5p","miR-320a","miR-320b","miR-185-5p","miR-21-5p","miR-629-5p", "periostina");
my %samples;
my @samples_array;

## Read healthy
my %values;

my $healthy_file = $datadir ."healthy.matrix";
open(HEALTH, $healthy_file) or die "Cannot open file $healthy_file ERROR:$!";

while(<HEALTH>){
    chomp;
    my @cols = split/\t/;
    my $sample = shift @cols;

    if(scalar(@cols) == 8){
	push @cols, "NA";
    }

    unless($sample eq "#NAMES"){

	$samples{$sample} = "healthy";
	push @samples_array, $sample;

	for(my $i=0;$i<scalar(@cols);$i++){

	    my $value = $cols[$i];  
	    if($value eq ""){
		$value = "NA";
	    }
	    
	    $values{$headers[$i]}{$sample} = $value;	    
	}
    }
}
close(HEALTH);

## Read asthma sample
my @sets = ("gema1", "gema2", "gema3", "gema4");

foreach my $set (@sets){
    my $file = $datadir . $set .".matrix";
    open(GEMA, $file) or die "Cannot open file $file ERROR:$!";
    while(<GEMA>){
	chomp;
	my @cols = split/\t/;
	my $sample = shift @cols;

	if(scalar(@cols) == 8){
	    push @cols, "NA";
	}
	
	if($samples{$sample}){
	    print "Sample: ". $sample ." is repeated in ". $set ." and ". $samples{$sample} ."\n";
	    exit;
	}
	else{	    
	    unless($sample eq "#NAMES"){

		$samples{$sample} = $set;
		push @samples_array, $sample;
		
		for(my $i=0;$i<scalar(@cols);$i++){
		    
		    my $value = $cols[$i];    
		    if(!$value){
			$value = "NA";
		    }
	    
		    $values{$headers[$i]}{$sample} = $value;
		}
	    }
	}
    }
    close(GEMA);
}

## Write the matrix
my $results_file = $datadir ."expression_matrix.txt";
open(RES, ">$results_file") or die "Cannot open file $results_file ERROR:$!";

print RES "#NUMBER_FEATURES\t10\n";
print RES "#NUMBER_SAMPLES\t". scalar(keys %samples) ."\n";
print RES "#VARIABLE\tclass\tCATEGORICAL{healthy,gema1,gema2,gema3,gema4}\tVALUES{";

my $sample_string;

foreach my $s (@samples_array){
    $sample_string .= $samples{$s} .",";
}
chop($sample_string);

print RES $sample_string ."}\tDESCRIPTION{}\n";

print RES "#NAMES";

foreach my $s (@samples_array){
    print RES "\t". $s;
}
print RES "\n";

foreach my $miRNA (@headers){
    print RES $miRNA;
    foreach my $sample (@samples_array){
	#print $miRNA ."\t". $sample ."\t". $values{$miRNA}{$sample} ."\n";
	print RES "\t". $values{$miRNA}{$sample};
    }
    print RES "\n";
}
close(RES);
