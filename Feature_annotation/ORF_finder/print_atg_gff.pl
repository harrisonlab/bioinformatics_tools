#!/usr/bin/perl -w
use strict;
use warnings;

#--------------------------------------------------------------------------
# print_atg_50.pl
#
#
#
#
# Usage: ./print_atg_50.pl <input_fasta_file> direction <aa_outfile> <nucleotide_outfile> <gff_outfile>
#
# Written by Joe Win, Kamoun lab, OSU-OARDC
#
#	Modified by Andrew Armitage, East Malling Research:
#				* extract predicted	nucleotide sequence in addition to aa's
#				* output a gff of ORFs
#				* Do not predict ORFs containing Xs
#
#--------------------------------------------------------------------------



#
# hash for translating DNA
#
my %DNAtoAA = ('GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TGT' => 'C',
	       'TGC' => 'C', 'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
	       'TTT' => 'F', 'TTC' => 'F', 'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G',
	       'GGG' => 'G', 'CAT' => 'H', 'CAC' => 'H', 'ATT' => 'I', 'ATC' => 'I',
	       'ATA' => 'I', 'AAA' => 'K', 'AAG' => 'K', 'TTG' => 'L', 'TTA' => 'L',
	       'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'ATG' => 'M',
	       'AAT' => 'N', 'AAC' => 'N', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
	       'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'CGT' => 'R', 'CGC' => 'R',
	       'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TCT' => 'S',
	       'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	       'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T', 'GTT' => 'V',
	       'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TGG' => 'W', 'TAT' => 'Y',
	       'TAC' => 'Y', 'TAA' => 'Z', 'TAG' => 'Z', 'TGA' => 'Z',
	       'ACN' => 'T', 'CCN' => 'P', 'CGN' => 'R', 'CTN' => 'L',
		   'GCN' => 'A', 'GGN' => 'G', 'GTN' => 'V', 'TCN' => 'S');


my $thisLine = "";
my $oneLongSeq = "";
my $thisSeqIndex = 0;
my $inSequence = 0;
my $seqComplete = 0;
my $totalSeqCount = 0;
my $orf_count = 0;
my $tempSeq = "";
my @sequences;
my $input = shift;
my $direction = shift;
my $aa_out = shift;
my $nuc_out = shift;
my $gff_out = shift;
open(INP, $input) || die "Cannot open file \"$input\"\n\n";
open(AA_OUT, ">$aa_out") || die "Cannot create file \"$aa_out\"\n\n";
open(NUC_OUT, ">$nuc_out") || die "Cannot create file \"$nuc_out\"\n\n";
open(GFF_OUT, ">$gff_out") || die "Cannot create file \"$gff_out\"\n\n";

my $gff1_name = '';
my $gff2_source = 'print_atg.pl';
my $gff3_type = 'CDS';
my $gff4_start = '';
my $gff5_end = '';
my $gff6_score = '.';
my $gff7_strand = '';
my $gff8_phase = '0';
my $gff9_attributes = '\"ID = ;\"';


while (<INP>) {
    chomp;
    $thisLine = $_;
    if ($thisLine =~ /^>/) {
        if ($inSequence) {
            $inSequence = 0;
            $seqComplete = 1;
        } else {
            $inSequence = 1;
            $seqComplete = 0;
        }
        $thisLine .= "\n";
    } else {
        $inSequence = 1;
    }

    if ($inSequence) {
        $tempSeq .= "$thisLine";
    }

    if ($seqComplete) {
        $tempSeq =~ /^(>.*)\n/;
		my $header = $1;
		print_atg_50 ($header, $tempSeq, 1, $direction);
        $tempSeq = "$thisLine";
        $seqComplete = 0;
        $inSequence = 0;
        $totalSeqCount++;
    }
}
$tempSeq =~ /^(>.*)\n/;
my $header = $1;
print_atg_50 ($header,$tempSeq, 1, $direction);
$totalSeqCount++;

exit;


#--------------------------------------------------------------------------
# Subroutines
#--------------------------------------------------------------------------

sub print_atg_50 {
	my ($this_header, $thisSeq, $frame, $direction) = @_;
	$thisSeq =~ s/>.*?\n//;
	$thisSeq =~ tr/[actg]/[ACTG]/;
	my $this_contig = $thisSeq;
	my $contig_length = length ($thisSeq);
	while ($thisSeq =~ /[ATCG]*?(ATG\w+?$)/) {
		my $this_frame = $1;
		$this_contig =~ /$this_frame/;
		my $before_atg = $`;
		my $feature_start = length ($before_atg);
		my $seq_length = length ($this_frame);
		if ($seq_length < 210) {
			last;
		}
		my ($peptide, $nuc) = translate_from_atg($this_frame);
		my $feature_end = ("$feature_start" + (length ($nuc)));
		if (length($peptide) >= 50 and $peptide !~ /X/) {
			$orf_count++;
			print AA_OUT "$this_header","_","$direction","$frame","\n";
			print AA_OUT "$peptide\n";
			print NUC_OUT "$this_header","_","$direction","$frame","\n";
			print NUC_OUT "$nuc\n";
			$gff1_name = substr ($this_header, 1);
			if ($direction eq 'F') {
				$gff7_strand = '+';
				$gff4_start = $feature_start;
				$gff5_end = $feature_end;
				$gff4_start++;
				$gff5_end++;
			} elsif ($direction eq 'R') {
				$gff7_strand = '-';
				$gff4_start = ("$contig_length" - "$feature_end");
				$gff5_end = ("$contig_length" - "$feature_start");
			} else {$gff7_strand = '?';}
			my $feature_name = "ORF_$direction" . "_$orf_count";
			my $feature_id = "$gff1_name" . "_$direction" . "$frame";
			$gff9_attributes = "\"ID = $feature_name\"; \"Name = $feature_id\"";
			print GFF_OUT "$gff1_name\t$gff2_source\t$gff3_type\t$gff4_start\t$gff5_end\t";
			print GFF_OUT "$gff6_score\t$gff7_strand\t$gff8_phase\t$gff9_attributes\n";
			$thisSeq = substr($this_frame, 3);
			$frame++;
		} else {
			$thisSeq = substr($this_frame, 3);
		}
	}
	return;
}

# #--------------------------------------------------------------------------
#
 sub translate_from_atg {
 	my $this_seq = shift;
 	$this_seq =~ tr/[a-z]/[A-Z]/;
 	my $pept;
 	my $nuc;
 	for (my $y = 0; $y < (length($this_seq) - 3); $y += 3) {
 		my $cur_codon = substr($this_seq, $y, 3);
 		if (!defined $DNAtoAA{$cur_codon}) {
 			$pept .= "X";
 			$nuc .= $cur_codon;
			last
 		} else {
 			if ($DNAtoAA{$cur_codon} eq "Z") {last}
 			$pept .= $DNAtoAA{$cur_codon};
			$nuc .= $cur_codon;
 		}
 	}
 	return ("$pept", "$nuc");
}
#
# #--------------------------------------------------------------------------
