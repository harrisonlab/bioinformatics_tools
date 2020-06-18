#!/usr/bin/perl
use strict;
use warnings;

# extract_gff_for_sigP_hits.pl
# Written by A.D Armitage, East Malling Research.
# This script will accept a file of gene names and additional information from
# fasta file headers that relate to CDS features in an augustus .gff3 file
# and will output a new feature track for the hit transcript, using the features supplied.
#


# build a hash of gene names present in the gene_list file.
# Associated information detailing the sigP Hmm score and cleavage site location
# is stored in the hash table.

# read each line of gff file
#	split the gff feature line
#	if the third line contatins 'gene'
#		then store this as the "previous_gene"
# 	if the third column contains 'transcript'
# 		then split the 9th column to get the transcript ID
#		if the transcript ID matches an entry in the hash then
#			print the previous gene as the parent feature (in gff3 format)
# 			print array of line
# 			...and column 2 with the specifed program...
# 			...and column 9 containing 'transcript_id "gene_name+transcript_no"; gene_id "gene_name"; parent "previous_gene_name";'
# 	repeat for next gff line

my $usage = "gene_list_to_gff.pl <gene_names.txt> <gene_models.gff> [prediction_source] [transcript_id_to_search_for_(ie.ID/Name)]\n\n";

my $infile = shift or die $usage;
my $gene_models = shift or die $usage;
my $prediction_source = shift or die $usage;
my $feature = shift or die $usage;
my %gene_hash;
my $transcript_ID;
my $gene_name;
my $gff_note;
my $gff_line;

my $gene_num;
my $Prev_num;
my $PrevCol1;
my $PrevCol2;
my $PrevCol3;
my $PrevCol4;
my $PrevCol5;
my $PrevCol6;
my $PrevCol7;
my $PrevCol8;
my $PrevCol9;

# This is a switch to turn on printing subfeatures with parents
# that are in the list of desired transcripts.
my $print_subfeatures = 0;

# read each line of gff file
#	split the gff feature line
#	if the third line contatins 'gene'
#		then store this as the "previous_gene"
# 	if the third column contains 'transcript'
# 		then split the 9th column to get the transcript ID
#		if the transcript ID matches an entry in the hash then
#			print the previous gene as the parent feature (in gff3 format)
# 			print array of line
# 			...and column 2 with the specifed program...
# 			...and column 9 containing 'transcript_id "gene_name+transcript_no"; gene_id "gene_name"; parent "previous_gene_name";'
# 	repeat for next gff line



print "##gff-version 3\n";

open (INFILE, $infile) or die "can't open: $infile\n$usage\n\n";
open (GENE_MODELS, $gene_models) or die "can't open: $gene_models\n$usage\n\n";

while (my $header_line = <INFILE>) {
  chomp $header_line;
  my @split_line = split ('\t', $header_line);
  # print "$split_line[0]\n";
  $gene_name = shift @split_line;
  chomp $gene_name;
  $gene_name =~ s/ $//g;
  # print "$gene_name\n";
  # $gene_name .= '.t1';
  my $gffnote = join(' ', @split_line);
  $gffnote =~ s/=/:/g;
  $gffnote =~ s/ //g;
  $gffnote =~ s/\t//g;
  if ($gene_hash{$gene_name}) {
    $gene_hash{$gene_name} .= ",$gffnote";
  } else {
    $gene_hash{$gene_name} = "Note=$gffnote";
    # print "$gene_name\t$gene_hash{$gene_name}\n";
  }
}



while ($gff_line = <GENE_MODELS>) {
	chomp $gff_line;
  # print "$gff_line\n";
	my @split_line = split ('\t', $gff_line);
	if ($split_line[0] =~ m/gff-version 3/) { next;
	} elsif ($split_line[2] eq 'gene') {
      # print "monkeys\n";
			$PrevCol1 = $split_line[0];
			$PrevCol2 = $prediction_source;
			$PrevCol3 = "gene";
			$PrevCol4 = $split_line[3];
			$PrevCol5 = $split_line[4];
			$PrevCol6 = $split_line[5];
			$PrevCol7 = $split_line[6];
			$PrevCol8 = $split_line[7];
			$PrevCol9 = $split_line[8];
			$print_subfeatures = 0;
	} elsif ($split_line[2] eq 'transcript' || $split_line[2] eq 'mRNA') {
    # print "badgers\n";
		$split_line[8] =~ m/$feature=.*/;
		$transcript_ID = $&;
		$transcript_ID =~ s/^.*?=//;
		$transcript_ID =~ s/;.*//;
		if ($gene_hash{$transcript_ID} && ($split_line[2] eq 'transcript' || $split_line[2] eq 'mRNA')) {
 			$transcript_ID =~ m/\.t[\d]\b/;
			my $gene_num = $`;
			my $col1 = $split_line[0];
			my $col2 = $prediction_source;
			my $col3 = $split_line[2];
			my $col4 = 	$split_line[3];
			my $col5 = 	$split_line[4];
			my $col6 = 	$split_line[5];
			my $col7 = 	$split_line[6];
			my $col8 = 	$split_line[7];
			my $col9 = 	"$split_line[8]" . ";$gene_hash{$transcript_ID}";
			print join("\t", $PrevCol1, $PrevCol2, $PrevCol3, $PrevCol4, $PrevCol5, $PrevCol6, $PrevCol7, $PrevCol8, $PrevCol9), "\n";
			print join("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9), "\n";
      # print "monkeys\n";
      # print "$gene_hash{$transcript_ID}";
			$print_subfeatures = 1;
		}
	} elsif ($print_subfeatures == 1 && index $split_line[8], "$transcript_ID") {
			my $gene_num = $&;
			my $col1 = $split_line[0];
			my $col2 = $prediction_source;
			my $col3 = $split_line[2];
			my $col4 = 	$split_line[3];
			my $col5 = 	$split_line[4];
			my $col6 = 	$split_line[5];
			my $col7 = 	$split_line[6];
			my $col8 = 	$split_line[7];
			my $col9 = 	$split_line[8];
			print join("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9), "\n";
	}
}

# while (my $gff_line = <GENE_MODELS>) {
# 	chomp $gff_line;
#   my @split_line = split ('\t', $gff_line);
#   if ($split_line[0] =~ m/gff-version 3/) { next;
#   } elsif ($split_line[2] eq 'CDS') {
#     $split_line[8] =~ m/$feature = \w+/;
#     $transcript_ID = $&;
#     $transcript_ID =~ s/^.*?= //;
#     if ($gene_hash{$transcript_ID} && $split_line[2] eq 'CDS') {
#       my $gene_num = $`;
#       my $col1 = $split_line[0];
#       my $col2 = $prediction_source;
#       my $col3 = $split_line[2];
#       my $col4 = $split_line[3];
#       my $col5 = $split_line[4];
#       my $col6 = $split_line[5];
#       my $col7 = $split_line[6];
#       my $col8 = $split_line[7];
#       my $col9 = "$split_line[8]" . "; \"$gene_hash{$transcript_ID}\"";
#       print join("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9), "\n";
#     }
#   }
# }

# print %gene_hash;

exit;
