#!/bin/env perl

use FindBin;
use lib $FindBin::Bin;

use Data::Dumper;
use Ort;
use strict;

my $file_name 	= $ARGV[0];
my $orts_format = $ARGV[1];
my $orgs_list 	= $ARGV[2];
my $exclusive 	= $ARGV[3];

my $usage =
"orts_pretty_print.pl <orts file> <orts format: 'RBH', 'RBH_Calhoun', 'OMCL'> <orgs separated by comma> [ <exclusive=1, non-exclusive=0> or filename listing the NOT set]\n\n";

die $usage if scalar(@ARGV) != 4;

my @orgs = split ",", $orgs_list;
my $num_orgs = scalar(@orgs);


# If third parameter different than 0 and 1, assume that it refers to a filename containing the NOT set
my %not_set;
if( $exclusive ne '1' && $exclusive ne '0' ){
	open IN, "$exclusive" or die "Unable to open $exclusive!\n";
	while(<IN>){
		my $line = $_;
		chomp $line;
		$not_set{$line} = 1; 	
	}
	close(IN);	
}

print STDERR "Reading orts file...\n";
my $orts = Ort::new( $file_name, 'transcript', $orts_format );
$orts->read();
print STDERR "DONE\n";

foreach my $cluster_num ( $orts->get_cluster_ids() ) {
	print STDERR "Cluster num: $cluster_num...\n";
	#getc();

	my %selected;
	my %all;
	my %not_selected;

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
			$not_selected{$curr_org} = 1 if( defined $not_set{$curr_org} );
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
	my $num_orts_not_selected = scalar( keys %not_selected );

	print "Num. orgs selected: $num_orgs   num. selected found: $num_orts_selected      num. non-selected found: $num_orts_not_selected      num. orgs in the cluster: $num_orts_all\n";

	# If EXCLUSIVE

	# If there are orts for the selected orgs
	# and these are the only one having orts then...

	if ($exclusive eq '1') {
		if ( $num_orts_all == $num_orgs && $num_orgs == $num_orts_selected ){
			print $orts->get_cluster_desc_repo_format($cluster_num);
		}else{			
		}
	}else {
		if ( $num_orgs == $num_orts_selected && $num_orts_not_selected == 0 ) {
		#if ( $num_orgs == $num_orts_selected ) {
			print $orts->get_cluster_desc_repo_format($cluster_num);
		}
	}
}

