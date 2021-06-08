#!/usr/bin/perl -w
use strict;

#S2_OPTiMuS. Perl script file for OPTiMuS.

# OPTimus script as described in Nodvig et al (2015)

#Interface, accepting either zero, two or three arguments. 
#If arguments are given from the command line, the order is: The Fastafile, the setting and if treshold or find is chosen a treshold value or protospacer
my $file;
my $setting;
my $item;
if (scalar(@ARGV)==0){
	print "Please input the name of a fastafile\n";
	$file = <STDIN>;
	chomp $file;
	print "Please input a setting, can be \"threshold\", \"find\" or \"setcover\"\n";
	$setting = <STDIN>;
	chomp $setting;
	$setting = lc($setting);
	if ($setting =~ m/^threshold$/ ) {
		print "Please input a number, for how many organisms a protospacer should be match before being displayed\n";
		$item = <STDIN>;
		chomp $item;
	}
	elsif ($setting =~ m/^find$/){
		print "Please input a protospacer (20bp) that you want to search your sequences for\n";
		$item = <STDIN>;
		chomp $item;
	}
	elsif ($setting =~ m/^setcover$/){
	}
	else {
		die "You didn't input a valid setting, try again\n"
	}
}
elsif (scalar(@ARGV==2) or scalar(@ARGV==3)){
	$file = $ARGV[0];
	$setting = $ARGV[1];
	$setting = lc($setting);	
	if ($setting =~ m/^find$/ or $setting =~ m/^threshold$/) {
		$item = $ARGV[2];
	}
	elsif ($setting =~ m/^setcover$/){
	}
	else {
		die "You didn't input a valid setting, try again\n"
	}
}
else {
	die "Invalid number of arguments\n";
}
if ($setting ne "setcover") {
	unless ($item =~m/^\d+$/ or $item =~m/^[atgcATGC]{20}$/ or $item =~ m/^$/) {
		die "Your third argument needs to be either a number or a protospacer (without PAM)\n";
	}
}
#Reading fastas into two arrays, one with the headers and one with the sequence
open (IN, '<', $file) or die "$!\n";

my $linecount = 0;
my (@header, @sequence);
my $line = <IN>;
while (defined($line)){
	$line =~ s/[\r\n]+$//;
	chomp $line;
	push(@header, $line);
	$linecount++;
	my $seq;
	while (defined($line = <IN>) and $line !~m/^>.+/){ #check if line is a header or sequence
		chomp $line;
		if ($line =~m/^[atgcATGC\s]+$/){ 
			$seq .= $line;
			$linecount++;
		}
		elsif ($line =~m/^$/){ 
		}
		else {
			die "Illegal characters in sequence at line $linecount\n";
		}
	}
	$seq =~s/\s//g;
	$seq = uc($seq);
	push (@sequence, $seq);
}

close IN;

#Ensure numbers of headers and sequences match
unless (scalar(@header)==scalar(@sequence)){
	die "number of headers does not match number of sequences\n";
}

#shortening the headers, removing everything after the first space
foreach my $headline (@header){
	my @tmp = split /\s/, $headline;
	$headline = $tmp[0];
}

#When the threshold setting is chosen
if ($setting eq "threshold" and $item =~ m/^\d+$/){
	#Finding all possible protospacers, in both direction and storing them in a hash of arrays, to find overlapping ones
	my %crispr;
	for (my $i= 0; $i < scalar(@sequence); $i++){
		for (my $j =0; $j < length($sequence[$i]) -23; $j++){
			my $tmp = substr($sequence[$i], $j, 23);
			if ($tmp =~m/^([ATGC]{20})[ATGC]GG$/){
				my $coordinate = $j+1;
				my $headline = $header[$i]."_coor-$coordinate";
				if (exists $crispr{$1}){
					push (@{$crispr{$1}}, $headline);
				}
				else {
					my @tmp;
					push (@tmp, $headline);
					$crispr{$1} = [@tmp];
				}
			}
			if ($tmp =~m/^CC[ATGC]([ATGC]{20})$/){
				my $coordinate = $j+4;
				my $headline = $header[$i]."_coor-$coordinate";
				my $revcomp = reverse($1);
				$revcomp =~tr/ATGC/TACG/;
				if (exists $crispr{$revcomp}){
					push (@{$crispr{$revcomp}}, $headline);
				}
				else {
					my @tmp;
					push (@tmp, $headline);
					$crispr{$revcomp} = [@tmp];
				}
			}
		}
	}

	#Printing the protospacers and their corresponding headers, which are above the threshold
	my $count = 0;
	foreach my $protospacer (sort {scalar(@{$crispr{$b}}) <=> scalar(@{$crispr{$a}})} keys %crispr) {
		unless (scalar(@{$crispr{$protospacer}} < $item)  ){
			print "$protospacer:\n";
			$count++;
			@{$crispr{$protospacer}} = sort @{$crispr{$protospacer}};
			foreach my $organism ( @{$crispr{$protospacer}} ){
				print "$organism\n";
			}
		}
	}
	if ($count == 0){
		print "Sorry no protospacer matching a number of organism above the threshold could be found, try with a lower threshold\n";
	}
}

#When the find setting is chosen
if ($setting eq "find" and $item =~m/^[atgcATGC]{20}$/){
	$item = uc($item);
	#Searches for an existing protospacer and prints the organisms it matches
	#Seaches each fasta for a match and store the hits in a new array
	my @hits;
	for (my $i= 0; $i < scalar(@sequence); $i++){
		for (my $j =0; $j < length($sequence[$i]) -23; $j++){
			my $tmp = substr($sequence[$i], $j, 23);
			if ($tmp =~m/^([ATGC]{20})[ATGC]GG$/){
				if ($1 eq $item) {
					my $match = $header[$i]."_pos: $j";
					push (@hits, $match);
				}
			}
			if ($tmp =~m/^CC[ATGC]([ATGC]{20})$/){
				my $revcomp = reverse($1);
				$revcomp =~tr/ATGC/TACG/;
				if ($revcomp eq $item){
					my $coordinate = $j+4;
					my $match = $header[$i]."_coor-$coordinate";
					push (@hits, $match);
				}
			}
		}
	}
	#prints the hits after sorting the organisms lexically
	if (scalar(@hits)>0) {
		print "Your protospacer matched following organisms:\n";
		@hits = sort @hits;
		foreach my $organism (@hits){
			print "$organism\n";
		}
	}
	else {
		print "Sorry, your protospacer did not match any of your input sequences";
	}
}

#When the setcover setting is chosen
#Finding all possible protospacers, in both direction and storing them in a hash of arrays
if ($setting eq "setcover")  {
	my %crispr;
	for (my $i= 0; $i < scalar(@sequence); $i++){
		for (my $j =0; $j < length($sequence[$i]) -23; $j++){
			my $tmp = substr($sequence[$i], $j, 23);
			if ($tmp =~m/^([ATGC]{20})[ATGC]GG$/){
				if (exists $crispr{$1}){
					push (@{$crispr{$1}}, $header[$i]);
				}
				else {
					my @tmp;
					push (@tmp, $header[$i]);
					$crispr{$1} = [@tmp];
				}
			}
			if ($tmp =~m/^CC[ATGC]([ATGC]{20})$/){
				my $revcomp = reverse($1);
				$revcomp =~tr/ATGC/TACG/;
				if (exists $crispr{$revcomp}){
					push (@{$crispr{$revcomp}}, $header[$i]);
				}
				else {
					my @tmp;
					push (@tmp, $header[$i]);
					$crispr{$revcomp} = [@tmp];
				}
			}
		}
	}
	my @crisprhits;
	my %refcrispr = %crispr;
	my %reforganism;
	
	#Solves the setcover problem with a Greedy-Algorithm
	while (scalar(@header)>scalar(keys %reforganism)){
		my @candidates = ();
		my $max = 0;
		foreach my $protospacer (keys %crispr) { #checking each protospacer and its underlying organism hits against the reference organism list, deleting those that aren't relevant
			
			#For each organism, if it is in the reference list delete it
			for (my $i = 0; $i < scalar(@{$crispr{$protospacer}}); $i++){
				my $organism = $crispr{$protospacer}[$i];
				if (defined $organism) {
					if (exists $reforganism{$organism}){
						delete $crispr{$protospacer}[$i];
					}
				}
			}
		
			#Finding the protospacer which matches the most unmatched organisms
			if (scalar(@{$crispr{$protospacer}})>=$max){ 
				$max = scalar(@{$crispr{$protospacer}});
				@candidates=();
				push (@candidates, $protospacer); # Empties and adds more candidates if a protospacer with more organisms is found
			}
		}
		
		
		
		#Finding the hits matching the most unmatched organism, and if several matches an equal amount of organisms, 
		#then compare already matched organisms for greater coverage. After that it's arbitary
		my $hit = ();
		if (scalar(@candidates) == 1) {
			$hit = $candidates[0];
			push (@crisprhits, $hit);
		}
		elsif (scalar(@candidates) > 1){
			my $altmax = 0; 
			foreach my $candidate (@candidates) {
				if ($candidate ne "") {
					if (exists $refcrispr{$candidate}) {
						if (scalar(@{$refcrispr{$candidate}}) > $altmax) {
							
							$altmax = scalar(@{$refcrispr{$candidate}});
							$hit = $candidate;
						}
					} 
				}
			}
			if (defined $hit) {
				push (@crisprhits, $hit);
			}
		}
		
		#deleting each organism covered by the new hit from the organism reference list
		#The hit being the best protospacer matching the most unmatched organisms
		if (defined $hit) {
			foreach my $organism (@{$crispr{$hit}}) {
				if (defined $organism) {
					if (!exists $reforganism{$organism}){
						$reforganism{$organism} = 1;
					}
				}
			}
		}
	
	}
	
	foreach my $str (@crisprhits) {
		print $str, "\n";
	}
}

__END__