#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions :config no_ignore_case);
use Pod::Usage;

# convert metadata file to input format for run_treeWAS.R.

=head1  SYNOPSIS

 metadata_to_treewas_input.pl -i [input_file] -c [column_header] -f [field1,field2] -o [path/to/output_file] [opt-args]
	
 -i|--input	input file containing sample id (first row) and 
		metadata column of interest [required]
 -c|--column	column header of metadata field of interest [required]
 -f|--fields	comma delimited list of row values matching an entry 
		in --column [required]
 -o|--output	path to output file [required]
 -d|--delimiter	delimiter used to seperate fields in input file
		[default: \t] 
 -h|--help 	usage information
=cut

# variables 
my $help = 0;

my $file = "";
my $header = "";
my $fields = "";
my $output = "";
my $delim = "\t";

GetOptions(
	'help|?' 	=> \$help,
	'input=s' 	=> \$file,
	'column=s' => \$header,
	'fields=s' => \$fields,
	'output=s' => \$output,
	'delimiter=s' => \$delim,
) or pod2usage(1);
pod2usage(1) if $help == 1;

# check for arguements
die " - ERROR: no input file specified\n" if $file eq "";
die " - ERROR: no column specified\n" if $header eq ""; 
die " - ERROR: no column fields specified\n" if $fields eq ""; 
die " - ERROR: no output file specified\n" if $output eq ""; 

# check input and output files exist/can be created
die " - ERROR: input file ($file) could not be created\n" if !( -e $file);
open OUT, ">$output" or die " - ERROR: output file ($output) could not be created\n";

# variables
my $count = 0;
my $header_idx = ""; 
my @isolate_list = ();
my %meta_store = ();

# split fields
my @f = split(/,/, $fields);

# parse file
open FILE, $file or die $!;
while(<FILE>){

	++$count; 
	
	my $line = $_;
	chomp $line;
	my @line = split( /$delim/, $line );

	# Header - find header field
	if( $count == 1 ){
		
		for my $i ( 0..$#line ){
			$header_idx = $i if ( $line[$i] eq $header );
		}
		
		# sanity check
		die "No header found.\n" if $header_idx eq "";
		
	}
	# Info line
	else{
		
		# Store isolate name
		push(@isolate_list, $line[0]);
		
		#print "$header_idx\t$line[0]\t$line[$header_idx]\n";
		
		# Check if isolate has metadata.
		for my $t (@f){
			if( $line[$header_idx] eq $t ){
				$meta_store{$line[0]} {$t} = 1;
				last;
			}
		}
	
	}
	
}

# print to output
my @headers = join("\t", @f);
print OUT "id\t@headers\n";

for my $i (@isolate_list){
	
	my @out_line = ();
	push( @out_line, $i );		

	for my $t (@f){
		if($meta_store{$i}{$t}){			
			push(@out_line, "1" );						
		}else{
			push(@out_line, "0");
		}
	}
	
	my $out = join("\t", @out_line);
	print OUT "$out\n";
	 
}

print "finished\n";

exit;

