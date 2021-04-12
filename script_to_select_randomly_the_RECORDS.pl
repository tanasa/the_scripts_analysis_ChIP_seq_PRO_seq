#!/usr/bin/perl

###################################################################
###################################################################
# The scripts reads a set of genes and selects 837 random entries

use warnings;
use POSIX;
use strict;

my $inputfile = '';
my $outputfile = '';

my $line = '';
my $range = 0;

my @v = ();
my $i = 1; 
my @random_numbers = ();
my @lines = ();

$inputfile = $ARGV[0];

$outputfile = $inputfile.".randomly.837.genes.txt";

# Read the INPUT file 

open(INPUT, "< $inputfile") or die "cannot open $inputfile : $!\n";
open(OUTPUT, "> $outputfile");

my $j = 1;
my @taken = (); 

###################################################################
###################################################################
###################################################################

while(defined($line = <INPUT>)) 
{
	chomp $line;
	$lines[$j] = $line;
        $taken[$j] = 0; 
	$j = $j + 1;
}

$range = $#lines ;

###################################################################
###################################################################
###################################################################

my $k = 0;
my $random_number = 0;

while ($k < 837)  ### select the number of 837 RANDOM ENTRIES ....
{
    $random_number = int(rand($range));
    
    if ($taken[$random_number] == 0)
    {    
        print OUTPUT "$lines[$random_number]\n";
        $k = $k + 1;
        $taken[$random_number] = 1;
    }
}    
    
close (INPUT);
close (OUTPUT);
exit;

###################################################################
###################################################################
###################################################################
###################################################################
