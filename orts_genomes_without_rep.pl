#!/bin/env perl

use Data::Dumper;
use Ort;
use strict;


my $file_name = $ARGV[0];


my $orts = Ort::new( $file_name );
$orts->read();

print Dumper( $orts->get_orts( $gene_name, $org_name ) );

print "\n";