#!/home/gomeza/miniconda3/envs/perly_env/bin/perl

# use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;

# taken from https://github.com/ISUgenomics/common_scripts/blob/master/gff2fasta.pl
# after inspiration from https://www.biostars.org/p/46281/ and modified as needed.
# 14/04/16 - MOdified to allow writing of final gene in gff file.

GetOptions( 'fasta=s' => \my $file_fasta
          , 'gff=s' => \my $gff
          , 'prefix=s' => \my $prefix
          , 'ranges=s{1,}' => \my @ranges
          );

# print "@ranges\n";
print "Selected ranges are:";
for (@ranges) {print "$_\n"}
# die;
# my $range_min = 1;
# my $range_max = 100;
# my $outfile = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].upstream${range_min}-${range_max}.fasta" );
my %outhash;
# $outlines{$range_max} = [split];

###### Output type description ######
# cds - translated sequence (starting with ATG and ending with a stop codon included)
# cdna - transcribed sequence (devoid of introns, but containing untranslated exons)
# protein - cds translated (includes a * as the stop codon)
# gene - the entire gene sequence (including UTRs and introns)
# upstream3000 - the 3000 upstream region of the gene (likely including the promoter)

### First, index the genome
# my $file_fasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");

my @frame_array;
my $in_frame_start;
my $in_frame_stop;
my $type;

### Second, parse the GFF3
my %CDS;
my %CDNA;
my %EXON;
my $mRNA_name;
my $strand;
my @array;
my $first_gene = 'true';
my $first_mRNA = 'true';
open GFF, "$gff" or die $!;
while ( my $line = <GFF> ) {
    chomp $line;
    # my @array = split( "\t", $line );
    @array = split( "\t", $line );
    $type = $array[2];

    if ($type eq 'gene' || $type eq 'mt_gene' ) {
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        if ($first_gene eq 'true') {
          $strand = $array[6];
          $first_gene = 'false';
        }
        my $gene_name = $attrs[0];
        my $contig_name = $array[0];
        my $gene_start = $array[3];
        my $gene_end = $array[4];

        # subroutine to extract upstream region
        # my $range_min = 1;
        # my $range_max = 100;
        for (@ranges) {
          @split_range = split(':', $_);
          $range_min = $split_range[0];
          $range_max = $split_range[1];
          extract_upstream($contig_name, $gene_name, $gene_start, $gene_end, $strand, $range_min, $range_max);
        }
    }
}
close GFF;



sub extract_upstream{
      my $contig_name = shift;
      my $gene_name = shift;
      my $gene_start = shift;
      my $gene_end = shift;
      my $strand = shift;
      my $range_min = shift;
      my $range_max = shift;
      # print "@_\n";
      # print @_[4];
        # The upstream 3000
        my $upstream_start;
        my $upstream_end;
        if($strand eq '+') {
            $upstream_start=$gene_start-$range_max;
            $upstream_end=$gene_start-$range_min;
            # if ($gene_name eq 'g2856') {print "$upstream_start\t$upstream_end\n";}
            # correct for genes on the start/end of contigs
            if ($upstream_start < 1) {$upstream_start = 1}
            if ($upstream_end < 1) {$upstream_end = 1}
        }
        elsif ($strand eq '-') {
            $upstream_start=$gene_end+$range_min;
            $upstream_end=$gene_end+$range_max;
            my $contig_lgth = length($db->seq($array[0]));
            if ($upstream_end > $contig_lgth) {$upstream_start = $contig_lgth}
            if ($upstream_start > $contig_lgth) {$upstream_start = $contig_lgth}
        }
        my $upstream_seq = $db->seq( $array[0], $upstream_start, $upstream_end);
        my $output = Bio::Seq->new(
            -seq        => $upstream_seq,
            -id         => $gene_name,
            -display_id => $gene_name,
            -alphabet   => 'dna',
        );

        # Reverse Complement if the frame is minus
        if($strand eq '+') {
        }
        elsif ($strand eq '-') {
            $output = $output->revcom();
        }
        else {
            die "Unknown frame! At line $. of the GFF\n";
        }
	#added an if statement for all outputs requiring there to be sequence information before writing to file otherwise the fasta file contains lots of empty fasta headers
          if (length($upstream_seq) > 1) {
            #
            # @outlines = $outlines{$range_max};
            # push @outlines, $output;
            # $outhash{$range_max}  = $output;
            # $outfile->write_seq($output);
            # print $output->id."\n".$output->seq()."\n";
            $outhash{$range_max} .= ">".$output->id."\n".$output->seq()."\n";
      	}
  }

for (@ranges) {
  @split_range = split(':', $_);
  $range_min = $split_range[0];
  $range_max = $split_range[1];
  open(my $fh, '>', "$prefix.upstream${range_min}-${range_max}.fasta");
  print $fh $outhash{$range_max};
  close $fh;
}
