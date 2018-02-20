#!/usr/bin/env perl

# make plottable table for plot_sig_hits.R

use strict;
use warnings;

# input/output
my $pirate = $ARGV[0];
my $treewas = $ARGV[1];
my $output = $ARGV[2];

# column index for PIRATE
my $idx = 19;

# parse PIRATE file
my @headers = ();
my @samples = ();
my $no_samples = "";

my %allele_hash = ();
my %group_hash = ();

my $gene_count = 0;
open INPUT, $pirate or die "Input file did not open.\n";
while(<INPUT>){
	
	my $line=$_;
	chomp $line;
	
	my @line = split(/\t/, $line);
	
	# get genome names
	if(/^allele_name/){
		
		@headers = @line;
		@samples = @headers[19..$#headers];
		$no_samples  = scalar(@samples);
		
	}else{
	
		++$gene_count;
	
		# sanity check 
		die " - ERROR: header not found in file" if @headers == 0; 
			
		# variables		
		my $a_name = $line[0];
		my $g_name = $line[1];
		
		# store presence/absence
		for my $i ($idx..$#line){
			$allele_hash{$a_name}{$headers[$i]} = 1 if $line[$i] ne "";
			$group_hash{$g_name}{$headers[$i]} = 1 if $line[$i] ne "";
		}		
	}
		
}close INPUT;

# parse treewas output
my %plot_hash = ();
my %var_type = ();
open TREEWAS, $treewas or die $!;
while(<TREEWAS>){
	
	my $line = $_;
	chomp $line;
	
	my @vars = split("\t", $line, -1);
	
	$var_type{$vars[0]} = $vars[$#vars];
	
	unless (/^SNP_locus/){
		my @inc = ();
		if ( $allele_hash{ $vars[0] } ){
			@inc = keys( %{ $allele_hash{$vars[0]}} );
		}elsif( $group_hash{ $vars[1] } ){
			@inc = keys (%{ $group_hash{$vars[0]}} );
		}else{
			die " - could not find group for $vars[0]-$vars[1]\n";
		}
		
		# add out plot hash
		for my $i (@inc) { $plot_hash{$vars[0]}{$i}=1 };
	}
}close TREEWAS;

# print to file
my @variants = sort keys (%plot_hash); 

# open output file
open OUT, ">$output" or die $!;

# add col headers
print OUT join("\t", "id", @variants), "\n";

# print variant presence/absence
for my $s (@samples){
	
	# store
	my @out_line = ($s);
	for my $v (@variants){
		my $out_val = "NA";
		$out_val = $var_type{$v} if $plot_hash{$v}{$s};
		push(@out_line, $out_val);
	}
	
	# print
	print OUT join("\t", @out_line), "\n";
	
}

exit
