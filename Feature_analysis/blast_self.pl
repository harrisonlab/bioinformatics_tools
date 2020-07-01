#!/home/gomeza/miniconda3/envs/perly_env/bin/perl
use strict;
use Cwd;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;

# Blast a file of candidate effectors against themselves to establish 
# homology between query sequences

#-------------------------------------------------------
# 		Step 1.		Initialise values
#-------------------------------------------------------

my $usage = "blast_self.pl <query_file.fa> <blastn/blastp> > <outfile.csv>";
my $query_file = shift or die $usage;
my $blast_type = shift or die $usage;
my @keys;
my $blast_fac;
my $result_obj;
my $report_obj;
my %hash;
my $hit_id;
my @ao_seqs;
my $seq;
my $read_no = 0;
my $outline_edit;
my @outline_start;
my @ao_arrays;
my %outline_hash;
my $alphabet;


#-------------------------------------------------------
# 		Step 2.		Create BLAST factory
#-------------------------------------------------------
$blast_fac = Bio::Tools::Run::StandAloneBlastPlus->new(
												'-db_name' => 'query_db', 
												'-db_dir' => '.', 
												'-create' => 1, 
												'-overwrite' => 1,  
												'-db_data' => $query_file,
												'-no_throw_on_crash' => 1,
												);

$blast_fac->make_db;												
#-------------------------------------------------------
# 		Step 3.		Collect sequence names from input
#------------------------------------------------------- 
if ($blast_type eq 'blastn') {$alphabet = 'dna'} elsif ($blast_type eq 'blastp') {$alphabet = 'protein'} else {die "$usage"}
my $seq_obj = Bio::SeqIO->new('-file' => $query_file, '-format' => "fasta", '-alphabet' => $alphabet );

while (my $seq = $seq_obj->next_seq) {push @ao_seqs, $seq};

#-------------------------------------------------------
# 		Step 3.		Create initial output line
#------------------------------------------------------- 

foreach (@ao_seqs) {
	$read_no ++;
	my $id = $_->id;
	push @keys, $id; 
	$hash{"$id"} = "$read_no";
	push @outline_start, '-';
}

#-------------------------------------------------------
# 		Step 3.	Create hash to reference for sequence order
#------------------------------------------------------- 

print "header\t";
foreach (@keys) {print "$_\t";}
print "\n";

#-------------------------------------------------------
# 		Step 4.		Re-open input
#-------------------------------------------------------

$seq_obj = Bio::SeqIO->new('-file' => $query_file, '-format' => "fasta", '-alphabet' => $alphabet );

#-------------------------------------------------------
# 		Step 4.		Perform BLAST, collect hits, for each
#					hit collect the name, reference the hash
#					and modify the outline accordingly.
#-------------------------------------------------------

while (my $seq = $seq_obj->next_seq) {
	my @outline_edit;
	my @ao_homologs;
	my $id = $seq->id;
	@outline_edit = @outline_start;
	unshift @outline_edit, $id;
	$report_obj = $blast_fac->run('-method' => $blast_type, '-query' => $seq, '-method_args' => ['-evalue' => 1e-10]);
  	my @ao_hits = $report_obj->hits;
  	foreach (@ao_hits) {				
 		my $hit = $_;
  		if ($hit) {			
 			my $hit_id = substr $hit->name(), 4;
 			my $hit_element = $hash{$hit_id};
			push @ao_homologs, $hit_element;
			splice @outline_edit, "$hit_element", 1, '1';	 			
		}
  	}
	
	foreach (@outline_edit) {
		print "$_\t";
	}		
	print "\n";
}

exit;