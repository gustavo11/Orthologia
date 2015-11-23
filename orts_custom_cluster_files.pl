#!/bin/env perl

use FindBin;
use lib $FindBin::Bin;

use Data::Dumper;
use Ort;
use strict;

my $file_name = $ARGV[0];
my $orgs_list = $ARGV[1];
my $exclusive = $ARGV[2];

my $usage =
"orts_pretty_print.pl <orts file> <orgs separated by comma> <exclusive=1, non-exclusive=0>\n\n";

die $usage if scalar(@ARGV) != 3;

my @orgs = split ",", $orgs_list;
my $num_orgs = scalar(@orgs);

print STDERR "Reading orts file...\n";
my $orts = Ort::new( $file_name, 'transcript' );
$orts->read();
print STDERR "DONE\n";

foreach my $cluster_num ( $orts->get_cluster_ids() ) {
	print STDERR "Cluster num: $cluster_num...\n";
	#getc();

	my %selected;
	my %all;

	# Iterate through all orgs
	my @all_orgs = $orts->get_orgs();
	for my $curr_org (@all_orgs) {

		#print STDERR "Org: $curr_org ...\n";
		#getc();

		my @orts =
		  $orts->get_orts_by_cluster_num_org( $cluster_num, $curr_org );

		if ( $orts[0] == 0 ) {

			#print STDERR "No ort of $curr_org in cluster $cluster_num";
			#getc();
		}
		else {
			$all{$curr_org} = 1;
		}

	}

	# Iterate through selected orgs
	for my $curr_org (@orgs) {
		my @orts =
		  $orts->get_orts_by_cluster_num_org( $cluster_num, $curr_org );

		if ( $orts[0] == 0 ) {

			#print STDERR "No ort of $curr_org in cluster $cluster_num";
			#getc();
		}
		else {
			$selected{$curr_org} = 1;
		}

	}

	my $num_orts_selected = scalar( keys %selected );
	my $num_orts_all      = scalar( keys %all );

#print STDERR "Num. orgs selected: $num_orgs   presence among selected: $num_orts_selected      presence: $num_orts_all\n";
#getc();

	# If EXCLUSIVE

	# If there are orts for the selected orgs
	# and these are the only one having orts then...

	if ($exclusive) {
		if ( $num_orts_all == $num_orgs && $num_orgs == $num_orts_selected ){
			print $orts->get_cluster_desc_repo_format($cluster_num);
		}else{			
		}
	}else {
		if ( $num_orgs == $num_orts_selected ) {
			print $orts->get_cluster_desc_repo_format($cluster_num);
		}
	}
}

