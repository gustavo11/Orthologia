#!/bin/env perl


use Data::Dumper;
use Ort;
use strict;

my $file_name = $ARGV[0];

my $cluster_num = $ARGV[1];

my $usage = "orts_pretty_print.pl <orts file> <cluster num>\n\n";

die $usage if scalar( @ARGV ) != 2;

my $orts = Ort::new( $file_name );
$orts->read();

	my @orgs = $orts->get_orgs();
	for my $curr_org ( @orgs ){
		my @orts = $orts->get_orts_by_cluster_num_org( $cluster_num, $curr_org );
		print "Org. $curr_org:\t";
		for my $curr_ort ( @orts ){
			print "$curr_ort,";
		}
		print "\n";
	}
		
	




print "\n";