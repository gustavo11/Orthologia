#!/bin/env perl

use strict;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Getopt::Std;
use GFFFile;
use Ort;

my $MAX_E_VALUE    = 1000;
my $DEFAULT_WIDTH  = 1500;
my $DEFAULT_LENGTH = 1500;
my $DEFAULT_COLOR  = "white";
my $LEFT_BORDER    = 0;
STDOUT->autoflush(1);

my $usage =
    "orts2graph.pl <upstream gene> <downstream gene> <cobra repository config file> <order genomes file>";

die $usage if ( scalar(@ARGV) != 4 );

my $up_gene   = $ARGV[0];
my $down_gene = $ARGV[1];

my $mauve_aln     = $ARGV[2];
my $mauve_bbone   = $ARGV[3];
my $GFF_list_file = $ARGV[4];

my $orts_file = $ARGV[5];

my $order_orts_file = $ARGV[6];

my $start_coord       = $ARGV[7];
my $end_coord         = $ARGV[8];
my $outputFile        = $ARGV[9];
my $protein_desc_file = $ARGV[10];

my $no_projections = $ARGV[11];

#my $addtional_features_file       = $ARGV[9];
my $additional_projections_file   = $ARGV[10];



my $width = $DEFAULT_WIDTH;

# Parsing pos and neg group
my @pos_group = split ",", $pos_group_str;
my @neg_group;
@neg_group = split ",", $neg_group_str if $neg_group_str ne "NULL";

my @order;
my @group;
for ( my $curr_index = 0 ; $curr_index < scalar(@pos_group) ; $curr_index++ ) {
	$pos_group[$curr_index] = $pos_group[$curr_index] - 1;
	push( @order, $pos_group[$curr_index] );
	$group[ $pos_group[$curr_index] ] = "+";
}

for ( my $curr_index = 0 ; $curr_index < scalar(@neg_group) ; $curr_index++ ) {
	$neg_group[$curr_index] = $neg_group[$curr_index] - 1;
	push( @order, $neg_group[$curr_index] );
	$group[ $neg_group[$curr_index] ] = "-";
}

# Instantiating Mauve object
my $mauveFile = MauveFile::new( $mauve_aln, $mauve_bbone );

open GFF_LIST, $GFF_list_file or die "Unable to open GFF file list $GFF_list_file!!";
my @gff_list = <GFF_LIST>;
close(GFF_LIST);

my @maxCoord;
my @minCoord;
my $overallMaxLength = -1;

my $has_orts_file = 1;

$has_orts_file = 0 if ( $order_orts_file eq "NULL" );

my $orts;
my @order_orts_param;
my @order_orts;

# Read Orts file
if ( $has_orts_file == 1 ) {
	$orts = Ort::new($orts_file);
	$orts->read();

	@order_orts_param = split ",", $order_orts_file;

	# Convert "genome order in orts file" from "genome names" to numbers
	if ( $order_orts_param[0] =~ ":" ) {

		foreach my $curr_item (@order_orts_param) {
			my ( $in_mauve, $name_in_orts ) = split ":", $curr_item;

			# If user indicates that genome is not present in orts
			if ( $name_in_orts == -1 ) {
				$order_orts[ $in_mauve - 1 ] = -1;
			}
			else {
				$order_orts[ $in_mauve - 1 ] = $orts->get_org_num($name_in_orts);
			}
		}

	}
	else {
		@order_orts = split ",", $order_orts_file;
		my $max = scalar(@order_orts);
		for ( my $index = 0 ; $index < $max ; $index++ ) {
			if ( $order_orts[$index] != -1 ) {
				$order_orts[$index] = $order_orts[$index] - 1;
			}
		}
	}
}
else {
	for ( my $curr_genome = 0 ; $curr_genome < $mauveFile->get_num_genomes() ; $curr_genome++ ) {
		$order_orts[$curr_genome] = 0;
	}
}
print STDERR "Genome mapping across files\n\n";
for ( my $curr_genome = 0 ; $curr_genome < $mauveFile->get_num_genomes() ; $curr_genome++ ) {
	my $fasta_file = $mauveFile->get_fasta_file_name($curr_genome);

	my $col_orts = "";
	$col_orts = $order_orts[$curr_genome] if $orts_file ne "NULL";

	my $name_in_orts_file = "";
	$name_in_orts_file = $orts->get_org_name($col_orts) if $orts_file ne "NULL";

	print STDERR "Internal num.: $curr_genome\n";
	print STDERR "\tMauve FASTA: $fasta_file\n";
	print STDERR "\tCol. in Orts file: $col_orts\n";
	print STDERR "\tName in Orts file: $name_in_orts_file\n\n";
}

#getc();

# print Dumper( $orts->get_orts( $gene_name, $org_name ) );
# print "\n";

my @genome_segment_length;

# Flags indicating some issues in each genome
my %different_contigs;  # Contig break within the homologous region
my %no_homology; 

# Detect the max and min cord
for ( my $curr_genome = 0 ; $curr_genome < $mauveFile->get_num_genomes() ; $curr_genome++ ) {

	next if $order_orts[$curr_genome] == -1;

	$minCoord[$curr_genome] = $start_coord;
	$maxCoord[$curr_genome] = $end_coord;
	if ( $curr_genome != 0 ) {
		$minCoord[$curr_genome] = $mauveFile->get_homologous_coord( 0, $start_coord, $curr_genome );
		$maxCoord[$curr_genome] = $mauveFile->get_homologous_coord( 0, $end_coord,   $curr_genome );
				
		# Checking if homology region was found
		if ( $minCoord[$curr_genome] == -1 && $maxCoord[$curr_genome] == -1){ 
			$no_homology{$curr_genome} = 1;
			print STDERR "NO homologous region in genome $curr_genome!!!\n"
		} 
	}

	my ( $contig_start, $contig_start_coord ) =
	  $mauveFile->master2contig_coord( $curr_genome, $minCoord[$curr_genome] );
	my ( $contig_end, $contig_end_coord ) =
	  $mauveFile->master2contig_coord( $curr_genome, $maxCoord[$curr_genome] );

	print STDERR
		"Genome $curr_genome  Segment start: $minCoord[ $curr_genome]  Segment end: $maxCoord[ $curr_genome ] \n";
	print STDERR
	  	"Genome contig: $contig_start($contig_start_coord)-$contig_end($contig_end_coord)\n";
	$genome_segment_length[$curr_genome] =
	  abs( $maxCoord[$curr_genome] - $minCoord[$curr_genome] ) + 1;

	if ( $contig_start ne $contig_end ) {
		$different_contigs{$curr_genome} = 1;
	}
	else {
		$overallMaxLength = $genome_segment_length[$curr_genome]
		  if ( $overallMaxLength < $genome_segment_length[$curr_genome] );
	}

}

#print STDERR "Overall length: $overallMaxLength\n";

my $bp_pixel_ratio = $overallMaxLength / $DEFAULT_WIDTH;

my @genome_panel_offset;
for ( my $curr_genome = 0 ; $curr_genome < $mauveFile->get_num_genomes() ; $curr_genome++ ) {
	next if $order_orts[$curr_genome] == -1;

	$genome_panel_offset[$curr_genome] =
	  int( ( $overallMaxLength - $genome_segment_length[$curr_genome] ) / 2 );
}

my $panel = Bio::Graphics::Panel->new(
	-length      => $overallMaxLength,
	-image_class => 'GD::SVG',
	-width       => $width,
	-pad_left    => 40,
	-pad_right   => 200,
	-pad_top     => 40,
	-pad_bottom  => 40,
);

my @genes_per_genome;

open PROT_DESC, ">$protein_desc_file" or die "Unable to open file $protein_desc_file for write !!";

my %features;

foreach my $curr_genome (@order) {
	next if $order_orts[$curr_genome] == -1;

	my $outside_genome_number = $curr_genome + 1;
	my $group_str             = $group[$curr_genome];
	my $genome_name           = $mauveFile->get_fasta_file_name($curr_genome);
	$genome_name =~ s/[\w\W]+\///;

	# If :
	# - there is a contig break within the homologous region
	# - there is no homology  
	# do not try to retrieve genome features from this genome.
	# 
	# An empty ruler will be drawn instead with a caption
	# indicating the corresponding issue
	
	# if contig break...
	if ( $different_contigs{$curr_genome} == 1 ) {
		my $full_length = Bio::SeqFeature::Generic->new(
			-start => 1,
			-end   => 1
		);

		$panel->add_track(
			$full_length,
			-glyph      => 'arrow',
			-arrowstyle => 'filled',
			-tick       => 1,
			-fgcolor    => 'black',
			-font2color => 'black',
			-double     => 1,
			-label => "$genome_name GROUP: $group_str (Homologous region contains a contig break)"
		);

		next;
		
	# if no homology
	}elsif( $no_homology{$curr_genome} == 1 ){
		my $full_length = Bio::SeqFeature::Generic->new(
			-start => 1,
			-end   => 1
		);

		$panel->add_track(
			$full_length,
			-glyph      => 'arrow',
			-arrowstyle => 'filled',
			-tick       => 1,
			-fgcolor    => 'black',
			-font2color => 'black',
			-double     => 1,
			-label => "$genome_name GROUP: $group_str (NO homology)"
		);

		next;
	}
		
	
	#####################################
	# Retrieving homologous coordinates
	
	my $effective_start = $start_coord;
	my $effective_end   = $end_coord;
	if ( $curr_genome != 0 ) {
		$effective_start = $mauveFile->get_homologous_coord( 0, $start_coord, $curr_genome );
		$effective_end   = $mauveFile->get_homologous_coord( 0, $end_coord,   $curr_genome );
	}

	# Adjusting MASTER START and END coord.. START should be lower than END
	if ( $effective_start > $effective_end ) {
		my $temp = $effective_end;
		$effective_end   = $effective_start;
		$effective_start = $temp;
	}

	my ( $contig_start, $contig_start_coord ) =
	  $mauveFile->master2contig_coord( $curr_genome, $effective_start );
	my ( $contig_end, $contig_end_coord ) =
	  $mauveFile->master2contig_coord( $curr_genome, $effective_end );

	my $reverse_complement_str = "";

	# Adjusting START and END coord. START should be lower value than END
	if ( $contig_start_coord > $contig_end_coord ) {
		my $temp = $contig_end_coord;
		$contig_end_coord       = $contig_start_coord;
		$contig_start_coord     = $temp;
		$reverse_complement_str = "(reverse complement)";
	}
	
	
	############################################
	# Adding track representing the genome and 
	# genomic features
	my $full_length = Bio::SeqFeature::Generic->new(
		-start => $genome_panel_offset[$curr_genome],
		-end   => $genome_panel_offset[$curr_genome] + $genome_segment_length[$curr_genome]
	);

	$panel->add_track(
		$full_length,
		-glyph                  => 'arrow',
		-arrowstyle             => 'filled',
		-tick                   => 1,
		-fgcolor                => 'black',
		-font2color             => 'black',
		-double                 => 1,
		-relative_coords        => 1,
		-relative_coords_offset => $contig_start_coord,
		#-label      => "$genome_name($outside_genome_number) GROUP: $group_str ; $effective_start - $effective_end From $contig_start:$contig_start_coord  to  $contig_end:$contig_end_coord"
		-label => "$genome_name GROUP: $group_str CONTIG: $contig_start $reverse_complement_str"
	);

	my $track;

	$track = $panel->add_track(
		-glyph   => 'box',
		-label   => 1,
		-bgcolor => sub {
			my $color = shift->{color};
			return $DEFAULT_COLOR if not defined($color);
			return $color;
		},
		-min_score  => 0,
		-max_score  => 0,
		-font2color => 'blue',
		-sort_order => 'high_score'

		  #		-description => sub {
		  #				    my $feature = shift;
		  #				    my $desc    = $feature->seq_id . " " . $feature->start . "-" . $feature->end;
		  #				    return $desc;
		  #				  }
	);
	
	
	#############################################
	# Retrieving genomic features from GFF file

	print STDERR "Opening GFF file $gff_list[ $curr_genome ] ...";
	my $gffFile = GFFFile::new( $gff_list[$curr_genome] );
	$gffFile->read();

	print STDERR "Retrieving from GFF genes in the range \n";
	print STDERR "\tmaster coord.: $effective_start-$effective_end ...\n";
	print STDERR
"\tcontig coord.: contig $contig_start ($contig_start_coord) - contig $contig_end ($contig_end_coord) ...\n";

	#getc();

	my $gffGenes = $gffFile->get_genes_hash();

	#my @list_genes_order = sort { $a->get_start() <=> $b->get_end() } $gffGenes;

	#for my $currGene (@list_genes_order){


	###############################################
	# Adding genomic features to its corresponding
	# track

	my $num_features_retrieved = 0;
	for my $currGene ( keys %{$gffGenes} ) {
		my $contig = $gffGenes->{$currGene}->get_chrom();
		my $id     = $gffGenes->{$currGene}->get_id();
		my $alias  = $gffGenes->{$currGene}->get_attribute("alias");
		my $name   = $gffGenes->{$currGene}->get_name();
		my $start  = $gffGenes->{$currGene}->get_start();
		my $end    = $gffGenes->{$currGene}->get_end();

		# Use alias if not empty
		# otherwise use id as gene identifier
		$id = $alias if ( $alias ne "" );

#if( $curr_genome != 1 ){
#	print STDERR ">>>>>Genome $curr_genome  Gene: $id Gene contig: $contig    Track start contig: $contig_start($contig_start_coord)   Track end contig: $contig_end($contig_end_coord)\n";
#	getc();
#}

		###################
		### WARNING
		### Quick fix for meeting un/comment next line

		next if ( $contig ne $contig_start && $contig ne $contig_end );
		next if ( $start < $contig_start_coord || $end > $contig_end_coord );

		push( @{ $genes_per_genome[$curr_genome] }, $id );

		###################
		### WARNING
		### Quick fix for meeting

		#################
		my $panel_adjusted_start =
		  $mauveFile->contig2master_coord( $curr_genome, $contig, $start ) - $effective_start;
		my $panel_adjusted_end =
		  $mauveFile->contig2master_coord( $curr_genome, $contig, $end ) - $effective_start;

		################
		#my $panel_adjusted_start =  $start - $effective_start;
		#my $panel_adjusted_end   =  $end - $effective_start;

		print PROT_DESC "$curr_genome\t$contig\t$id\t$name\n";

		if ( $id =~ "contig_break" ) {
			$features{$id}{biofeature} = Bio::SeqFeature::Generic->new(
				-score        => 100,
				-seq_id       => $id,
				-display_name => "Contig break -" . $name,
				-start        => $panel_adjusted_start + $genome_panel_offset[$curr_genome],
				-end          => $panel_adjusted_end + $genome_panel_offset[$curr_genome]
			);
			$features{$id}{biofeature}{color} = "red";
		}
		else {
			$features{$id}{biofeature} = Bio::SeqFeature::Generic->new(
				-score        => 100,
				-seq_id       => $id,
				-display_name => $id . "-" . $name,
				-start        => $panel_adjusted_start + $genome_panel_offset[$curr_genome],
				-end          => $panel_adjusted_end + $genome_panel_offset[$curr_genome]
			);
		}

		$features{$id}{genome} = $curr_genome;

		$track->add_feature( $features{$id}{biofeature} );
		$num_features_retrieved++;

	}

	print STDERR "Num. features retrieved: $num_features_retrieved\n\n";
}
close(PROT_DESC);

<<<<<<< .mine


###############################
#Additional features from file

#Format, separated by tab
#my $format_feat_file  =  "<name> <genome num.> contig:<contig name>:<start>-<end>       #if contig coord";
#$format_feat_file    .= "<name> <genome num.> <start>-<end>	                         #if master coord";

#open ADD_FEAT, $addtional_features_file or die "Unable to open file $addtional_features_file\n";
#my $cont_line;
#while(<ADD_FEAT>){
#	my $line = $_;
#	chomp $line;
#	my ($name,$contig_name,$start,$end);
	
	#if( my ($name,$genome,$contig_name,$start,$end) = ( $line =~ /^(\w\W]+?)\t(\d+)\tcontig:(\w\W]+?):(\d+)-(\d+)/ ) ){
		
	#}elsif( my ($name,$genome,$contig_name,$start,$end) = ( $line =~ /^(\w\W]+?)\tcontig:(\w\W]+?):(\d+)-(\d+)/ ){
	#}

		
	#}else{
	#	die "Line $cont_line of file $addtional_features_file does not meet contain the correct" .
	#	    " format\n$format_feat_file\n";
	#}
#} 



my $additional_projections_file   = $ARGV[10];


my @pos_group =  (3,4,5);
my @neg_group =  (0,1,2);

=======
>>>>>>> .r93
################################
## Coloring based on orhtology
################################

my $debug_coloring = 0;

print STDERR "\n\nColoring genes...\n" if $debug_coloring;
if ( $has_orts_file == 1 && scalar(@pos_group) != 0 && scalar(@neg_group) != 0 ) {

	foreach my $gene_id ( keys %features ) {
		next if ( $gene_id =~ "contig_break" );

		# Color
		my $color_soft_gene_deletion;
		my $color_hard_gene_deletion;
		my @verify_against;

		print STDERR "Coloring $gene_id from genome $features{$gene_id}{genome}...\n"
		  if $debug_coloring;

		if ( scalar( grep( /^$features{$gene_id}{genome}$/, @pos_group ) ) >= 1 ) {
			$color_soft_gene_deletion = "cyan";
			$color_hard_gene_deletion = "blue";
			@verify_against           = @neg_group;
			print STDERR "\tCompared against negative group ...\n" if $debug_coloring;
		}
		elsif ( scalar( grep( /^$features{$gene_id}{genome}$/, @neg_group ) ) >= 1 ) {
			$color_soft_gene_deletion = "yellow";
			$color_hard_gene_deletion = "orange";
			@verify_against           = @pos_group;
			print STDERR "\tCompared against positive group ...\n" if $debug_coloring;
		}
		else {
			next;
		}

		# Check for orthologs along the whole genome
		my $string_orgs = "";
		foreach my $curr_pos (@verify_against) {
			$string_orgs .= $orts->get_org_name( $order_orts[$curr_pos] ) . ":";
		}
		$string_orgs =~ s/:$//;

		print STDERR "\tLooking for orthologs in the complete genomes $string_orgs ...\n"
		  if $debug_coloring;

		if ( not $orts->has_at_least_one_ort( $gene_id, $string_orgs ) ) {
			$features{$gene_id}{biofeature}{color} = $color_hard_gene_deletion;

#print "2. Verify: $string_orgs Gene: $gene_id Genome:" . $features{$gene_id}{genome} . " color: $paint_in\n";
#getc();
			print STDERR "\t\tNO ORTS at all.\n\n" if $debug_coloring;

			#getc();
			next;
		}
		else {
			print STDERR "\t\tORT FOUND among at least one of those genomes\n" if $debug_coloring;
		}

		# Check for orthologs in the vicinity (syntenic orthologs)
		my $found = 0;
		print STDERR "\tLooking for orthologs in region depicted by the figure ...\n"
		  if $debug_coloring;

		# Iterate through genomes
		foreach my $curr_org (@verify_against) {
			my $org_name = $orts->get_org_name( $order_orts[$curr_org] );
			my @ort_genes = $orts->get_orts( $gene_id, $org_name );

			print STDERR "\tChecking genome $curr_org  $org_name ...\n" if $debug_coloring;

			# Checking if orhtolog is present in the figure
			# if that's the case skip coloring
			#if ( $ort_genes[0] != 0 ){
			foreach my $curr_ort (@ort_genes) {
				print STDERR "\t\tChecking if ortholog $curr_ort is in the figure ...\n"
				  if $debug_coloring;

				#found among the features that will be displayed
				if ( defined $features{$curr_ort} ) {
					$found = 1;
					print STDERR "\t\tFOUND among features in the figure\n\n" if $debug_coloring;
					last;
				}
			}

			#}
			last if $found == 1;
		}

		# If no orts found in the figure
		if ( $found == 0 ) {
			$features{$gene_id}{biofeature}{color} = $color_soft_gene_deletion;
			print STDERR "\t\tNOT FOUND among features in the figure\n\n" if $debug_coloring;
		}

		#getc();

	}

}

if ( $no_projections || $has_orts_file == 0 ) {
	open OUTPUT, ">" . $outputFile;
	print OUTPUT $panel->svg;
	close OUTPUT;
	exit;
}

###########################
# Drawing projections

my $gd    = $panel->gd;
my @boxes = $panel->boxes;

my $red   = $panel->translate_color('red');
my $blue  = $panel->translate_color('blue');
my $green = $panel->translate_color('springgreen');
my $white = $panel->translate_color('white');

# Parameters are R,G,B and alpha for transparency.
# R,G,B goes from 0 to 255 (complete saturated)
# The alpha value may range from 0 (opaque) to 127 (transparent).
# The alphaBlending function changes the way this alpha channel affects the resulting image.
my $projection_fill_color = $gd->colorAllocateAlpha( 125, 125, 125, 100 );
my $projection_outline    = $gd->colorAllocateAlpha( 80,  80,  80,  100 );

# Retrieve boxes
my %feature_box;
for my $box (@boxes) {
	my ( $feature, @points ) = @$box;

	my $gene_name = $feature->seq_id();
	next if not defined($gene_name) or $gene_name eq "";

	$points[1] = $points[1] + 10;
	$points[2] = ( $feature->length() / $bp_pixel_ratio ) + $points[0];

	@{ $feature_box{ $feature->seq_id() } } = @points;
}

print STDERR "\n\nDrawing projections ...\n";

for ( my $curr_index = 0 ; $curr_index <= $#order ; $curr_index++ ) {

	my $curr_genome = $order[$curr_index];

	#print STDERR "Src genome: $curr_genome\n";
	#getc();
	next if $order_orts[$curr_genome] == -1;

	foreach my $curr_gene ( @{ $genes_per_genome[$curr_genome] } ) {

		my $gene_name = $curr_gene;
		next if not defined($gene_name) or $gene_name eq "";

		#print STDERR "Gene name: $gene_name\n";
		next if not defined( $feature_box{$gene_name} );

		my @src_rectangle = @{ $feature_box{$gene_name} };

		#print STDERR "Src rectangle: $src_rectangle[ 0 ], $src_rectangle[ 1 ]\n";

		# Draw projection to ortholog in genome just below, if there are no orhtologs in that genome
		# try the next genome until it finds at least one ortholog

		for (
			my $curr_index_dest = $curr_index + 1 ;
			$curr_index_dest <= $#order ;
			$curr_index_dest++
		  )
		{

			my $dest_genome = $order[$curr_index_dest];

			#print STDERR "Dst genome: $dest_genome\n";
			#getc();
			if ( $order_orts[$dest_genome] == -1 ) {

		   #print STDERR "Skipping genome $dest_genome. Not present on cluster of orthologs file\n";
				next;
			}
			next if $different_contigs{$curr_genome} == 1;

			my $org_name = $orts->get_org_name( $order_orts[$dest_genome] );

			my @ort_genes;

			@ort_genes = $orts->get_orts( $gene_name, $org_name );

			#print STDERR "Length ort_genes:" . scalar( @ort_genes ) . "\n";

			#if( @ort_genes = $orts->get_orts( $gene_name, $org_name ) ){

			my $draw_at_least_one = 0;

			foreach my $curr_ort_gene (@ort_genes) {

				#print STDERR "Ort name: $curr_ort_gene\n";
				next if $curr_ort_gene eq "0";

				next if not defined( $feature_box{$curr_ort_gene} );
				my @dst_rectangle = @{ $feature_box{$curr_ort_gene} };

				#print STDERR "Dst rectangle: $dst_rectangle[ 0 ], $dst_rectangle[ 1 ]\n";

				my $poly = new GD::Polygon;
				$poly->addPt( $src_rectangle[0], $src_rectangle[3] );
				$poly->addPt( $dst_rectangle[0], $dst_rectangle[1] );

				$poly->addPt( $dst_rectangle[2], $dst_rectangle[1] );
				$poly->addPt( $src_rectangle[2], $src_rectangle[3] );

				$gd->filledPolygon( $poly, $projection_fill_color );
				$gd->polygon( $poly, $projection_outline );

				$draw_at_least_one = 1;

				#print STDERR $feature->seq_id();
				#print STDERR "  @points[0], @points[1]\n";
				#getc();
			}

			# Go to the next gene if at least one ortholog was found and draw
			last if $draw_at_least_one == 1;
		}
	}
}

open OUTPUT, ">" . $outputFile;
print OUTPUT $panel->svg;
close OUTPUT;

