#!/bin/env perl

use FindBin;
use lib $FindBin::Bin;

use Data::Dumper;
use Ort;
use strict;

my $file_name = $ARGV[0];

my $gene_name = $ARGV[1];

my $org_name  = $ARGV[2];

my $usage = "orts_query.pl <orts file> <gene name> [org name]\n\n";

die $usage if scalar( @ARGV ) != 3 && scalar( @ARGV ) != 2;

my $orts = Ort::new( $file_name );
$orts->read();

if( defined($org_name ) ){
	print Dumper( $orts->get_orts( $gene_name, $org_name ) );
}else{
	my @orgs = $orts->get_orgs();
	for my $curr_org ( @orgs ){
		my $org_num = $orts->get_org_num( $curr_org ) + 1;
		my @orts = $orts->get_orts( $gene_name, $curr_org );
		print "$org_num $curr_org:\t";
		for my $curr_ort ( @orts ){
			print "$curr_ort,";
		}
		print "\n";
	}
}
		
	




print "\n";