#!/usr/bin/perl -w

use strict;
use warnings;


#-------------------------------------------------------------------------
# gffexpander.pl
#-------------------------------------------------------------------------
# This script will accept features in a gff file and expand the region
# of the gff features by Xbp either upstream, downstream of both.
# This can be used with bedtools intersect to identify whether the newly 
# expanded features eg. promotors overlap nearby features of interest 
# eg. genes

my $usage = "gffexpander.pl <+,- or +-> <distanceto_expand(bp)> <feature_file.gff>";
#-------------------------------------------------------------------------
# 1. Initiate variables
#-------------------------------------------------------------------------

# Set whether extending upstream of the feature, downstream of the feature or both
my $UpDown = shift;
unless ($UpDown =~ m/\+|\-/) {die "$usage"};
# Set the no. bp to extend the feature by	
my $distance = shift;
unless ($distance =~ m/\d+/) {die "$usage"};
# Open input file
my $gffFile = shift;

#-------------------------------------------------------------------------
# 2. Open input &  collect sequences in two parallel arrays:
#	one for names and one for aa sequences.
#-------------------------------------------------------------------------

open (GFF_FILE, "$gffFile ") || die "Cannot open file \"$gffFile\"\n\n";

while (<GFF_FILE>) {
   	chomp;
   	my $thisLine = $_;
   	my @gffFeatureIn = split("\t", $thisLine);
   
   	my $col1 = shift @gffFeatureIn;			# Sequence id
	my $col2 = shift @gffFeatureIn;			# Source
	my $col3 = shift @gffFeatureIn;			# Type
	my $col4 = shift @gffFeatureIn;			# Start
	my $col5 = shift @gffFeatureIn;			# End
	my $col6 = shift @gffFeatureIn;			# Score
	my $col7 = shift @gffFeatureIn;			# Strand
	my $col8 = shift @gffFeatureIn;			# Phase
	my $col9 = "@gffFeatureIn";				# Attribute
	
    ($col4, $col5) = mod_gff($col4, $col5, $col7, $UpDown, $distance);	# $mimpStart, $mimpEnd, $strand
	print join ("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9) . "\n";
}		

close (GFF_FILE);
    
exit;
#-------------------------------------------------------------------------
# Sub 1. Modify .gff features.
#-------------------------------------------------------------------------

sub mod_gff {
	my ($mimpStart, $mimpEnd, $strand, $UpDown, $distance) = @_;
    my $upstream = 'no'; 
    my $downstream = 'no';
    if ($UpDown =~ m/\-/) { $upstream = 'yes'; }
    if ($UpDown =~ m/\+/) { $downstream = 'yes'; }
    if ($strand eq '+') {
    	if ($downstream eq 'yes') {
    		$mimpEnd += $distance;
    	}
    	if ($upstream eq 'yes') {
    		$mimpStart -= $distance;
    	}
    } elsif ($strand eq '-') {
    	if ($downstream eq 'yes') {
    		$mimpStart -= $distance;
    	}
    	if ($upstream eq 'yes') {
    		$mimpEnd += $distance;
    	}
    }
    if ($mimpStart <= 1) {
    	$mimpStart = 1;
    }
	return ($mimpStart, $mimpEnd);
}
  