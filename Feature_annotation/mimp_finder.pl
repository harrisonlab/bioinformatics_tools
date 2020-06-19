#!/usr/bin/perl -w
#
use strict;
use warnings;
#
#
#-------------------------------------------------------------------------
# find_rxlr.pl
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# 1. Initiate variables
#-------------------------------------------------------------------------
my $seqName = "";
my @seqNames = ();
my $dnaSeq = "";
my @dnaSeqs = ();
my $total = 0;
my $tempSeq=();
my $thisLine=();
my $inSequence = 0;
my $seqComplete = 0;

# Input file
my $fastaFile = shift;
chomp $fastaFile;
# Fasta Output
my $outFile = shift;
chomp $outFile;
# .gff Output
my $outGff = shift;
chomp $outGff;


#-------------------------------------------------------------------------
# 2. Open input &  collect sequences in two parallel arrays:
#	one for names and one for aa sequences.
#-------------------------------------------------------------------------

open (FASTA_FILE, "$fastaFile") || die "Cannot open file \"$fastaFile\"\n\n";

my $flag=0;

while (<FASTA_FILE>) {
        # collect sequences in two parallel arrays:
        # one for names and one for aa sequences.
    	    chomp;
    	$thisLine = $_;
    	
	if ($thisLine =~ /^>/ && $flag==0) {
		#print "START LINE DETECTED $thisLine \n";
		#$thisLine .= "\n";	
		push (@seqNames, $thisLine);
		$flag=1;
        	}
	elsif ( $thisLine =~ /^>/ && $flag==1){
		# print "START LINE DETECTED $thisLine \n";
		 push (@seqNames, $thisLine);
		 push (@dnaSeqs,$tempSeq);
		 $tempSeq=();
		}
	else {
        	$tempSeq .= "$thisLine";
    	}
	
}
#print "PUSHING FINAL LINE\n";
push (@dnaSeqs,$tempSeq);
close FASTA_FILE;

print "SEQUENCES IN NAME ARRAY  ".scalar(@seqNames)."\n";
print "SEQUENCES IN SEQS ARRAY  ".scalar(@seqNames)."\n";

#-------------------------------------------------------------------------
# 3. Open fasta and gff outfiles
#-------------------------------------------------------------------------

open (OUT, ">$outFile") || die "Cannot open file \"$outFile\"\n\n";
print OUT "Following FASTA files contain the mimp motif...\n";

open (GFF, ">$outGff") || die "Cannot open file \"$outGff\"\n\n";

#-------------------------------------------------------------------------
# 4. Work through arrays identifying mimps
#-------------------------------------------------------------------------
my $subSeq;
my $mimpCount = 0;
my $i = 0;
$total = @seqNames;



for ($i = 0; $i < $total; $i++) {
#for ($i = 140; $i < 141; $i++) {
    #print "Working on $seqNames[$i] \n";
    $seqName = $seqNames[$i];
    $dnaSeq = $dnaSeqs[$i];
  
	while ($dnaSeq =~/CAGTGGG..GCAA[TA]AA/g ){
		my $mimp_pos = pos($dnaSeq) ;
		print "FOUND MIMP on $seqName at $mimp_pos\n";
         	$mimpCount=$mimpCount+1;;
		print OUT "$seqName --mimp starts at $mimp_pos\.\n";
        	print OUT substr($dnaSeq, ($mimp_pos-16), 80), "\n";
		my $mimpStart = $mimp_pos-15;
		my $mimpEnd = $mimp_pos;
		my $strand = "+";
   	    print_gff($seqName, $mimpStart, $mimpEnd, $strand, $mimpCount);
	}
	while ($dnaSeq =~/TT[TA]TTGC..CCCACTG/g ){
		my $mimp_pos = pos($dnaSeq) ;
		print "FOUND MIMP on $seqName at $mimp_pos\n";
		$mimpCount=$mimpCount+1;;
		print OUT "$seqName --mimp starts at $mimp_pos\.\n";
			print OUT substr($dnaSeq, ($mimp_pos-16), 80), "\n";
		my $mimpStart = $mimp_pos-15;
		my $mimpEnd = $mimp_pos;
		my $strand = "-";
      	print_gff($seqName, $mimpStart, $mimpEnd, $strand, $mimpCount);
    }

    
}

print OUT "There are $mimpCount sequences that contain the consensus mimp motif\n";
print "There are $mimpCount sequences that contain the consensus mimp motif\n";
close (OUT);
close (GFF);

exit;


#-------------------------------------------------------------------------
# Sub 1. Print hit to gff.
#-------------------------------------------------------------------------

sub print_gff {
	my ($seqName, $mimpStart, $mimpEnd, $strand, $mimpCount) = @_;
	$seqName = substr $seqName, 1;			# Removes ">" from name

	my $col1 = "$seqName";					# Sequence id
	my $col2 = "mimp_finder.pl";			# Source
	my $col3 = "MIMP_motif";				# Type
	my $col4 = "$mimpStart";				# Start
	my $col5 = "$mimpEnd";					# End
	my $col6 = ".";							# Score
	my $col7 = "$strand";					# Strand
	my $col8 = ".";							# Phase
	my $col9 = "ID=\"MIMP_$mimpCount\"";	# Attributes

	print GFF join ("\t", $col1, $col2, $col3, $col4, $col5, $col6, $col7, $col8, $col9) . "\n";
}
