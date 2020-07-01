#!/home/gomeza/miniconda3/envs/perly_env/bin/perl
use strict;
use Cwd;


my $usage = "blast_parse.pl <homology_table.csv> > <outfile.csv>";
my $query_file = shift or die $usage;
my @ao_column;
my $group;
my %hash_tab;
my %out_hash;
my $column_no = 0;
my $group_no = 1;
my %seen;
my $headers_col;
my $itterations = 1;


#-------------------------------------------------------
# 		Step 1.		Collect inputs into hash
#-------------------------------------------------------

# Build a hash table from an array of inputs

open (INFILE, "$query_file") or die $usage;
#my @ao_inlines = ('header	gene1	gene2	gene3	gene4', 'gene1	1	-	-	-', 'gene2	-	1	-	-',	'gene3	-	-	1	1', 'gene4	-	-	1	1');
my @ao_inlines = <INFILE>;
foreach (@ao_inlines) {	
	my @inline = split ('\t', $_);
	$headers_col = shift @inline;
	$hash_tab{"$headers_col"} = "@inline";
}



#-------------------------------------------------------
# 		Step 2.		Find homolog groups
#-------------------------------------------------------
#	a-	Collect the hash table header and set as array.
# 	b-	Build an array of all the elements in the first array
#		in the hash table.
#	c-	If the column is identicle to a previous column, skip
#		it.
#	d-	Substitute the "hits" in the array to a homolog group.
#	e-	Create an array of modified elements.
#	f-	Append each element in the array to the next line
#		in a new hash table.

my @keys = split (' ', $hash_tab{'header'});	

foreach (@keys) {
	my $cur_column;
	foreach (@keys) {
		my @cur_line = split (' ', $hash_tab{"$_"});
		$cur_column .= "$cur_line[$column_no];";
#		print "$cur_column\n";
#		print "@cur_line\n";
	}
	$column_no ++;
	next if $seen{$cur_column}++;	# Skips if previously seen column
	next if (($cur_column =~ tr/1//) == 1);	# Skips if column contains a single element.
	$cur_column =~ s/1/$group_no/g;
	my @ao_elements = (split (';', "$cur_column"));
	$group_no ++;
	foreach (@keys) {
		$out_hash{"$_"} .= shift @ao_elements;
		$out_hash{"$_"} .= ";";
	}
}

#-------------------------------------------------------
# 		Step 3.		Print the final hash table
#-------------------------------------------------------


print "ID";
while ($itterations != $group_no){ print "\tGrp$itterations"; $itterations ++ ;}
print "\n";

foreach (@keys) {
	if ($out_hash{$_}) {
	my $outline = (join ("\t", (split (';', $out_hash{$_}))));
	print "$_\t$outline\n";
	} else {print "$_\n"}
}

exit;