#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum);
use Cwd;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Bio::DB::Fasta;

my $gff_f = shift or die;
my $fasta_f = shift or die;

my $gffio = Bio::Tools::GFF->new( -file => $gff_f, -gff_version => 3 );
my $fasta_db = Bio::DB::Fasta->new($fasta_f);

my @outlines;
my $feature;
my $feat_type;
my $location;

my $prev_seq_lgth;
my $prev_seq_id;

my $gene_start;
my $gene_end;
my $gene_strand;
my $gene_id;

my $mRNA_id;

my $exon_start;
my $exon_end;

my $prev_feat = '';

my $i = 1;
my $cds_id;
my $cds_gffstr;
my $cds_strand;
my $prev_cds_start;
my $prev_cds_end;

my $end_feat_B;
my $end_B_gffstr;
my $codon;


while($feature = $gffio->next_feature()) {
        if ($feature->primary_tag() eq 'transcript') {
          $feature->primary_tag('CDS')
        }
        my $seq_id = $feature->seq_id();
        my $source = $feature->source_tag();
        $feat_type = $feature->primary_tag();
        if ($feat_type eq 'gene' ) {
              $i = 1;
              $location = $feature->location();
              $gene_start = $location->start;
              $gene_end = $location->end;
              $gene_strand = $location->strand;
              if ($gene_strand eq '1') {
                $location->end(($gene_end +2));
                if (int($location->end()) >= int($fasta_db->length($seq_id))) {
                  # print ($location->end());
                  # print "\t";
                  # print($fasta_db->length($seq_id)+ "\n");
                  # print "\n";
                  $location->end($fasta_db->length($seq_id))}
              } elsif ($gene_strand eq '-1') {
                $location->start(($gene_start -2));
                if (int($location->start()) <= int(1)) {$location->start('1')}
              }
              # my $seqX=$fasta_db->seq($seq_id, $location->start, $location->end);
              # print("$seqX\n");
              $gene_id = join(', ', $feature->get_tag_values('ID'));
              $gene_id =~ tr/\./_/;
              $feature->remove_tag('ID');
              $feature->add_tag_value('ID', "$gene_id");
              my $gene_gffstr = $gffio->gff_string($feature);

              # Make A start or stop codon for the previous gene
              # print "$cds_strand\n";
              if ($prev_feat eq 'CDS' and $prev_cds_end != $prev_seq_lgth) {
                  $end_feat_B = Bio::SeqFeature::Generic->new( -gff_string => $cds_gffstr);
                  $end_feat_B->start($prev_cds_end -3);
                  $end_feat_B->add_tag_value('Parent', "$mRNA_id");
                  if ($cds_strand eq '-1') {
                    $end_feat_B->primary_tag("start_codon");
                    $codon = $fasta_db->seq($prev_seq_id, ($prev_cds_end), ($prev_cds_end -2));
                    $codon =~ tr/atg/ATG/;
                    # print "codon is: $codon\n";
                    if ($codon =~ /ATG/i) {
                      my $end_B_gffstr = $gffio->gff_string($end_feat_B);
                      print "$end_B_gffstr\n";
                      # print "$codon\n";
                    }
                  } elsif ($cds_strand eq '1') {
                    $end_feat_B->primary_tag("stop_codon");
                    $codon = $fasta_db->seq($prev_seq_id, ($prev_cds_end -2), $prev_cds_end);
                    $codon =~ tr/atgc/ATGC/;
                    # print "codon is: $codon\n";
                    if ($codon =~ /T[AA|AG|GA]/i) {
                      my $end_B_gffstr = $gffio->gff_string($end_feat_B);
                      print "$end_B_gffstr\n";
                      # print "$codon\n";
                    }
                  }
              }

              # Print the features of this gene
              print "$gene_gffstr\n";

              # Make an mRNA feaute for the gene
              my $mRNA_feat = Bio::SeqFeature::Generic->new( -gff_string => $gene_gffstr);
              $mRNA_id = "$gene_id.t1";
              $mRNA_feat->primary_tag("mRNA");
              $mRNA_feat->add_tag_value('ID', "$mRNA_id");
              $mRNA_feat->add_tag_value('Parent', "$gene_id");
              my $mRNA_gffstr = $gffio->gff_string($mRNA_feat);
              print "$mRNA_gffstr\n";
              $prev_feat = 'gene';

        # Modify CDS feautes
        } elsif ($feat_type eq 'CDS') {
            $location = $feature->location();
            $cds_strand = $location->strand;
            if ($gene_strand eq '1') {
              $location->end(($location->end() +2));
              if (int($location->end()) >= int($fasta_db->length($seq_id))) {$location->end($fasta_db->length($seq_id))}
            } elsif ($gene_strand eq '-1') {
              $location->start(($location->start() -2));
              if (int($location->start()) <= int(1)) {$location->start('1')}
            }
            my $this_cds_start = $location->start;
            my $this_cds_end = $location->end;
            my $cds_feat = $feature;
            $cds_feat->remove_tag('Parent');
            $cds_feat->remove_tag('ID');
            $cds_feat->add_tag_value('Parent', "$mRNA_id");
            my $cds_id = "$mRNA_id.CDS$i";
            $i =+1;
            $cds_feat->add_tag_value('ID', "$cds_id");
            $cds_gffstr = $gffio->gff_string($cds_feat);


            # Make a start/stop codon feature if this is the first CDS in a gene.
            if ($prev_feat eq 'gene' and $this_cds_start != 0) {
                my $end_feat_A = Bio::SeqFeature::Generic->new( -gff_string => $cds_gffstr);
                $end_feat_A->end($this_cds_start +3);
                $end_feat_A->add_tag_value('Parent', "$mRNA_id");
                if ( $cds_strand eq '1') {
                  $end_feat_A->primary_tag("start_codon");
                  $codon = $fasta_db->seq($seq_id, $this_cds_start, ($this_cds_start +2));
                  # print "Codon is: $codon\n";
                  $codon =~ tr/atgc/ATGC/;
                  if ($codon =~ /ATG/i) {
                    my $end_A_gffstr = $gffio->gff_string($end_feat_A);
                    print "$end_A_gffstr\n";
                  }
                } elsif ($cds_strand eq '-1') {
                  $end_feat_A->primary_tag("stop_codon");
                  $codon = $fasta_db->seq($seq_id, ($this_cds_start +2), $this_cds_start);
                  # print "$codon\n";
                  $codon =~ tr/atgc/ATGC/;
                  if ($codon =~ /T[AA|AG|GA]/i) {
                    my $end_A_gffstr = $gffio->gff_string($end_feat_A);
                    print "$end_A_gffstr\n";
                    # print "$codon\n";
                  }
                }


            # Make any introns if necessary:
            } elsif ($prev_feat eq 'CDS') {
              my $intron_start = $prev_cds_end +1;
              my $intron_end = $this_cds_start -1;
              my $intron_feat = Bio::SeqFeature::Generic->new( -gff_string => $cds_gffstr);
              $intron_feat->primary_tag("intron");
              $intron_feat->start($intron_start);
              $intron_feat->end($intron_end);
              $intron_feat->add_tag_value('Parent', "$mRNA_id");
              my $intron_gffstr = $gffio->gff_string($intron_feat);
              print "$intron_gffstr\n";
            }

            # Print the CDS feature:
            print "$cds_gffstr\n";

            # Make an exon fubfeature:
            my $exon_feat = Bio::SeqFeature::Generic->new( -gff_string => $cds_gffstr);
            $exon_feat->primary_tag("exon");
            $exon_feat->add_tag_value('Parent', "$mRNA_id");
            my $exon_id = "$mRNA_id.exon$i";
            $exon_feat->add_tag_value('ID', "$exon_id");
            my $exon_gffstr = $gffio->gff_string($exon_feat);
            print "$exon_gffstr\n";

            # Set variables for later use in start/stop codon and intron boundaries
            $prev_feat = 'CDS';
            $prev_cds_start = $this_cds_start;
            $prev_cds_end = $this_cds_end;
        }
        $prev_seq_lgth=$fasta_db->length("$seq_id");
        $prev_seq_id = $seq_id;
}

# Print a start or stop codon for the final gene
$end_feat_B = Bio::SeqFeature::Generic->new( -gff_string => $cds_gffstr);
$end_feat_B->start($prev_cds_end -3);
$end_feat_B->add_tag_value('Parent', "$mRNA_id");
if ($cds_strand eq '-1') {
  $end_feat_B->primary_tag("start_codon");
  $codon = $fasta_db->seq($prev_seq_id, ($prev_cds_end), ($prev_cds_end -2));
  $codon =~ tr/atg/ATG/;
  # print "codon is: $codon\n";
  if ($codon =~ /ATG/i) {
    my $end_B_gffstr = $gffio->gff_string($end_feat_B);
    print "$end_B_gffstr\n";
    # print "$codon\n";
  }
} elsif ($cds_strand eq '1') {
  $end_feat_B->primary_tag("stop_codon");
  $codon = $fasta_db->seq($prev_seq_id, ($prev_cds_end -2), $prev_cds_end);
  $codon =~ tr/atg/ATG/;
  if ($codon =~ /T[AA|AG|GA]/i) {
    my $end_B_gffstr = $gffio->gff_string($end_feat_B);
    print "$end_B_gffstr\n";
    # print "$codon\n";
  }
}
