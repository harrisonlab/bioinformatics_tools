#!/home/gomeza/miniconda3/envs/perly_env/bin/perl
use strict;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Search::Result::BlastResult;

# Perform a tBlASTn search of query genes against a file containing nucleotide data. It will then parse the results
# into a tab separated output

my $usage = "blast2csv.pl <query_file.fa> <blastn/tblastn/tblastx/blastp/blastx> <genomic_contigs.fa> <no.hits_to_report> > <outfile.csv>\n\n";
my $query_file = shift or die $usage;
my $blast_type = shift or die $usage;
my $database = shift or die $usage;
my $no_hits = shift or die $usage;
my $blast_fac;
my $result_obj;
my $report_obj;
my @ao_outlines;
my $seq_obj;
my $method;
my $alphabet;

if ($blast_type eq 'blastn') {$alphabet = 'dna'}
elsif ($blast_type eq 'blastx') {$alphabet = 'dna'}
elsif ($blast_type eq 'blastp') {$alphabet = 'protein'}
elsif ($blast_type eq 'tblastn') {$alphabet = 'protein'}
elsif ($blast_type eq 'tblastx') {$alphabet = 'dna'}
else {die "$usage"}


#-------------------------------------------------------
# 		Step 1.		Set outfile header
#-------------------------------------------------------

my $outline_header = "ID\tSequence\tSequence_lgth\tNo.hits";
for (1 .. $no_hits) {
#	print "$_";
	$outline_header .= "\tHit\tE-value\tHit_lgth\tPer_length\tPer_ID\tHit_strand\tHit_start\tHit_end\tHit_seq";
}
$outline_header .= "\n";
print "$outline_header";

#-------------------------------------------------------
# 		Step 2.		Create BLAST factory
#-------------------------------------------------------

$blast_fac = Bio::Tools::Run::StandAloneBlastPlus->new(
												'-db_name' => 'genome_db',
												'-db_dir' => '.',
												'-create' => 1,
												'-overwrite' => 1,
												'-db_data' => $database,
												'-no_throw_on_crash' => 1
												);

$blast_fac->make_db;


#-------------------------------------------------------
# 		Step 3.		Open query file
#-------------------------------------------------------

my $input_obj = Bio::SeqIO->new('-file' => $query_file, '-format' => 'fasta', '-alphabet' => $alphabet);

#-------------------------------------------------------
# 		Step 4.		Perform blast for each query
#-------------------------------------------------------

while (my $seq = $input_obj->next_seq) {

	$report_obj = $blast_fac->run('-method' => $blast_type, '-query' => $seq, '-method_args' => ['-evalue' => 1e-10]);
	my @ao_hits = $report_obj->hits;

#-------------------------------------------------------
# 		Step 5.		Get query info
#-------------------------------------------------------
 	my $seq_id =  $seq->id;
	my $sequence = $seq->seq;
	my $query_lgth = ($seq->length);
	my $num_hits = $report_obj->num_hits;

#-------------------------------------------------------
# 		Step 6.		Declare and set hit values
#-------------------------------------------------------

	my $hit_id;
	my $hit_lgth;
	my $per_query;
	my $ident;
	my $per_id;
	my $e_value;
 	my $hit_seq;
 	my $outline = "$seq_id\t$sequence\t$query_lgth\t$num_hits";
 	my $strand;
 	my $start;
 	my $end;
 	my @ao_hsps;
 	my $hit_itterations;
 	foreach (@ao_hits) {
 		my $hit = $_;
 		my $hsp = $hit->hsp('best') or last;
 	# 	$hit_id = substr $hit->name(), 4;
		$hit_id = substr $hit->name(),
		$hit_seq = $hsp->seq_str('hit');
 		$hit_lgth = length ($hit_seq);
 		if ($blast_type eq 'tblastx' or $blast_type eq 'blastx') { $hit_lgth = ($hit_lgth * 3); }
 		$per_query = ($hit_lgth / $query_lgth);
 		$per_query = substr $per_query, 0, 4;
		if ($blast_type eq 'tblastx' or $blast_type eq 'blastx') { $per_id = ($hsp->num_identical() / ($query_lgth / 3));
		} else { $per_id = ($hsp->num_identical() / $query_lgth); }
		$per_id = substr $per_id, 0, 4;
		$e_value = $hsp->evalue;
		$strand = $hsp->strand('hit');
		$start = $hsp->start('hit');
		$end = $hsp->end('hit');
		$outline .= "\t$hit_id\t$e_value\t$hit_lgth\t$per_query\t$per_id\t$strand\t$start\t$end\t$hit_seq";

	}

	$outline .= "\n";


#-------------------------------------------------------
# 		Step 7.		Print BLAST info to outfile
#-------------------------------------------------------

 	push @ao_outlines, $outline;

 }

 foreach (@ao_outlines) { print "$_";}

 exit;
