#!/home/gomeza/miniconda3/envs/perly_env/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

# taken from https://github.com/ISUgenomics/common_scripts/blob/master/gff2fasta.pl
# after inspiration from https://www.biostars.org/p/46281/ and modified as needed.
# 14/04/16 - MOdified to allow writing of final gene in gff file.

#add a help message here
#my $num_args=$#ARGV + 1;
#if ($num_args != 4) {
#	print "\nUsage: gff2perl Genome.fasta Annotation.gff OutputPrefix \n\n";
#	exit;
#}

$| = 1;    # Flush output
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cds.fasta" );
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].pep.fasta" );
my $outfile_cdna = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cdna.fasta" );
my $outfile_gene = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].gene.fasta" );
my $outfile_upstream3000 = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].upstream3000.fasta" );
my $outfile_exon = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].exon.fasta");

###### Output type description ######
# cds - translated sequence (starting with ATG and ending with a stop codon included)
# cdna - transcribed sequence (devoid of introns, but containing untranslated exons)
# protein - cds translated (includes a * as the stop codon)
# gene - the entire gene sequence (including UTRs and introns)
# upstream3000 - the 3000 upstream region of the gene (likely including the promoter)

### First, index the genome
my $file_fasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");

my $output_cdna;
# my $mergedCDNA_seq;
my $output_cds;
# my $mergedCDS_seq;
my @frame_array;
my $in_frame_start;
my $in_frame_stop;
my $out_cds_in_frame;
my $output_pep;
my $type;

### Second, parse the GFF3
my %CDS;
my %CDNA;
my %EXON;
my $mRNA_name;
my $frame;
my @array;
my $first_gene = 'true';
my $first_mRNA = 'true';
open GFF, "<$ARGV[1]" or die $!;
while ( my $line = <GFF> ) {
    chomp $line;
    # my @array = split( "\t", $line );
    @array = split( "\t", $line );
    $type = $array[2];

    if ($type eq 'gene' || $type eq 'mt_gene' ) {
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        if ($first_gene eq 'true') {
          $frame = $array[6];
          $first_gene = 'false';
        }
        my $gene_name = $attrs[0];
        my $gene_start = $array[3];
        my $gene_end = $array[4];
        my $gene_seq = $db->seq( $array[0], $gene_start, $gene_end );
        my $output_gene = Bio::Seq->new(
            -seq        => $gene_seq,
            -id         => $gene_name,
            -display_id => $gene_name,
            -alphabet   => 'dna',
        );


        # The upstream 3000
        my $upstream_start;
        my $upstream_end;
        if($array[6] eq '+') {
            $upstream_start=$gene_start-3000;
            $upstream_end=$gene_start-1;
        }
        elsif ($array[6] eq '-') {
            $upstream_start=$gene_end+1;
            $upstream_end=$gene_end+3000;
        }
        my $upstream_seq = $db->seq( $array[0], $upstream_start, $upstream_end );
        my $output_upstream3000 = Bio::Seq->new(
            -seq        => $upstream_seq,
            -id         => $gene_name."_upstream3000",
            -display_id => $gene_name."_upstream3000",
            -alphabet   => 'dna',
        );

        # Reverse Complement if the frame is minus
        if($array[6] eq '+') {
        }
        elsif ($array[6] eq '-') {
            $output_gene = $output_gene->revcom();
            $output_upstream3000 = $output_upstream3000->revcom();
        }
        else {
            die "Unknown frame! At line $. of the GFF\n";
        }
	#added an if statement for all outputs requiring there to be sequence information before writing to file otherwise the fasta file contains lots of empty fasta headers
        if (length($gene_seq) != 0) {
            	$outfile_gene->write_seq($output_gene);
      	}
      	if (length($upstream_seq) != 0) {
              $outfile_upstream3000->write_seq($output_upstream3000);
      	}
    }

#CDS
    if ( ( $type eq 'mRNA' || $type eq 'transcript' ) and ( $. > 2 ) ) {
        # CDS: Collect CDSs and extract sequence of the previous mRNA
        my $mergedCDS_seq;
	# WARNING we must sort by $cds_coord[1]
        # if ($first_mRNA eq 'true') {
        #   my @attrs = split( ";", $array[8] );
        #   $attrs[0] =~ s/ID=//;
        #   $mRNA_name = $attrs[0];
        #   $frame=$array[6];
        #   $first_mRNA = 'false';
        # }

        foreach my $key (sort {$a <=> $b} keys %CDS) { # Ascending numeric sort of the starting coordinate
            my $coord = $CDS{$key};
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mergedCDS_seq .= $cds_seq;
        }


        $output_cds = Bio::Seq->new(
            -seq        => $mergedCDS_seq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
        );
        if ($frame eq '-') {
            $output_cds = $output_cds->revcom();
            @frame_array = reverse @frame_array;
        }
      	#translate CDS to peptide for protein sequence

        $in_frame_start = ($frame_array[0] +1);
        $in_frame_stop = $output_cds->length();
        my $cds_subseq = $output_cds->subseq($in_frame_start, $in_frame_stop);

        $out_cds_in_frame = Bio::Seq->new(
            -seq        => $cds_subseq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
          );

        $output_pep = $out_cds_in_frame->translate();
        undef @frame_array;
        # $output_pep = $output_cds->translate();
      	#write to file
      	if (length($mergedCDS_seq) != 0) {
              $outfile_cds->write_seq($output_cds);
      	}
      	if (length($mergedCDS_seq) != 0) {
              $outfile_pep->write_seq($output_pep);
      	}

#exons
#should be able to add exon output here since exons will be useful in gene models for other organisms can be added in the EVM program
        my $mergedEXON_seq;
        foreach my $key (sort {$a <=> $b} keys %EXON) { # Ascending numeric sort of the starting coordinatg
            my $coord = $EXON{$key};
            my @exon_coord = split( " ", $coord );
            my $exon_seq = $db->seq( $exon_coord[0], $exon_coord[1], $exon_coord[2] );
            $mergedEXON_seq .= $exon_seq;
        }

       my $output_exon = Bio::Seq->new(
            -seq        => $mergedEXON_seq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
        );
        if ($frame eq '-') {
            $output_exon = $output_exon->revcom();
        }
	#write to file
        if (length($mergedEXON_seq) != 0) {
        $outfile_exon->write_seq($output_exon);
        }



        # CDNA: Collect UTRs and CDSs and extract sequence of the previous mRNA
        my $mergedCDNA_seq;
        foreach my $key (sort {$a <=> $b} keys %CDNA) { # Ascending numeric sort of the starting coordinate
            my $coord = $CDNA{$key};
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mergedCDNA_seq .= $cds_seq;
        }

        my $output_cdna = Bio::Seq->new(
            -seq        => $mergedCDNA_seq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
        );
        if ($frame eq '-') {
            $output_cdna = $output_cdna->revcom();
        }
        if (length($mergedCDNA_seq) != 0) {
          	$outfile_cdna->write_seq($output_cdna);
      	}


        # Now initialize the next mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
        %CDS = (); %CDNA = (); # Empty the chunk arrays
      	%EXON = (); %EXON = (); #Empty the EXON chunk arrays
    }
    elsif ( $type eq 'mRNA'  || $type eq 'transcript') {    # First mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
    }
    elsif ( $type eq 'CDS' ) {
        my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDS{$array[3]}=$cds_coord;
        $CDNA{$array[3]}=$cds_coord;
        push @frame_array, $array[7];
    }
    elsif ($type eq 'UTR' ) {
        my $utr_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDNA{$array[3]}=$utr_coord;
    }
    elsif ($type eq 'exon' ) {
      	my $exon_coord = $array[0] . " " . $array[3] . " " . $array[4];
      	$EXON{$array[3]}=$exon_coord;
    }
}


#--------
# Printing final mRNA and features
#--------

#CDS
# CDS: Collect CDSs and extract sequence of the previous mRNA
my $mergedCDS_seq;
foreach my $key (sort {$a <=> $b} keys %CDS) { # Ascending numeric sort of the starting coordinate
    my $coord = $CDS{$key};
    my @cds_coord = split( " ", $coord );
    my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
    $mergedCDS_seq .= $cds_seq;
}

$output_cds = Bio::Seq->new(
    -seq        => $mergedCDS_seq,
    -id         => $mRNA_name,
    -display_id => $mRNA_name,
    -alphabet   => 'dna',
);
if ($frame eq '-') {
    $output_cds = $output_cds->revcom();
}
#translate CDS to peptide for protein sequence
$in_frame_start = ($frame_array[0] +1);
$in_frame_stop = $output_cds->length();
my $cds_subseq = $output_cds->subseq($in_frame_start, $in_frame_stop);

$out_cds_in_frame = Bio::Seq->new(
    -seq        => $cds_subseq,
    -id         => $mRNA_name,
    -display_id => $mRNA_name,
    -alphabet   => 'dna',
  );

$output_pep = $out_cds_in_frame->translate();
undef @frame_array;

#write to file
if (length($mergedCDS_seq) != 0) {
      $outfile_cds->write_seq($output_cds);
}
if (length($mergedCDS_seq) != 0) {
      $outfile_pep->write_seq($output_pep);
}

#exons
#should be able to add exon output here since exons will be useful in gene models for other organisms can be added in the EVM program
my $mergedEXON_seq;
foreach my $key (sort {$a <=> $b} keys %EXON) { # Ascending numeric sort of the starting coordinatg
    my $coord = $EXON{$key};
    my @exon_coord = split( " ", $coord );
    my $exon_seq = $db->seq( $exon_coord[0], $exon_coord[1], $exon_coord[2] );
    $mergedEXON_seq .= $exon_seq;
}
my $output_exon = Bio::Seq->new(
    -seq        => $mergedEXON_seq,
    -id         => $mRNA_name,
    -display_id => $mRNA_name,
    -alphabet   => 'dna',
);
if ($frame eq '-') {
    $output_exon = $output_exon->revcom();
}
#write to file
if (length($mergedEXON_seq) != 0) {
$outfile_exon->write_seq($output_exon);
}


# CDNA: Collect UTRs and CDSs and extract sequence of the previous mRNA
my $mergedCDNA_seq;
foreach my $key (sort {$a <=> $b} keys %CDNA) { # Ascending numeric sort of the starting coordinate
    my $coord = $CDNA{$key};
    my @cds_coord = split( " ", $coord );
    my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
    $mergedCDNA_seq .= $cds_seq;
}

$output_cdna = Bio::Seq->new(
    -seq        => $mergedCDNA_seq,
    -id         => $mRNA_name,
    -display_id => $mRNA_name,
    -alphabet   => 'dna',
);
if ($frame eq '-') {
    $output_cdna = $output_cdna->revcom();
}
if (length($mergedCDNA_seq) != 0) {
  	$outfile_cdna->write_seq($output_cdna);
}

close GFF;
