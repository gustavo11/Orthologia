package Ort;
use strict 'vars';
use strict 'refs';

use Data::Dumper;

sub new {
	my ( $filename, $gene_identifying_feature ) = @_;
	my $self = { filename => $filename, gene_identifying_feature => $gene_identifying_feature };
	
	bless $self, Ort;
	
	return $self;
}

sub get_filename {
	my $self = shift;
	return 	$self->{filename};
}

sub read {
	my $self = shift;
	
	# Checking format of the file
	open ORTS_FILE, $self->{filename} or die "Unable to open file " . $self->{filename} . "\n";
	while( <ORTS_FILE> ){
		my $line = $_;
		chomp $line;
		
		my @cols = split "\t", $line;
		if( $cols[0] eq "#CLUSTER ID" ){
			$self->{format} = "one_cluster_per_line";
			last;
		}else{
			$self->{format} = "one_gene_per_line";
			last;
		}
				
	}
	close( ORTS_FILE );
	
	$self->read_one_cluster_per_line() if $self->{format} eq "one_cluster_per_line";
	$self->read_one_gene_per_line()   if $self->{format} eq "one_gene_per_line";
	
	
}

sub read_one_cluster_per_line {
	my $self = shift;	
		
	my $first_line = 1;

   die "This format is not being used anymore due to non-unique identification of genes " .
       "associated to it (genes from different org with same identifier)\n\n";
			
	open ORTS_FILE, $self->{filename} or die "Unable to open file " . $self->{filename} . "\n";
	while( <ORTS_FILE> ){
		my $line = $_;
		chomp $line;
		
		# Read header, disregarding first 4 columns
		if( $first_line ){
			my @cols = split "\t", $line;
			
			for( my $col_num = 4; $col_num <= $#cols; $col_num++ ){
				my $col_content = $cols[ $col_num ];
				my $org_num = $col_num - 4;
				
				$self->{org_num}{$col_content} = $org_num;	 
				$self->{org_name}[ $org_num ]  = $col_content;
			}
			$first_line = 0;
			next;
		}
		my @cols = split "\t", $line;
		
		my $orts_num =  $cols[ 0 ];
		push( @{$self->{cluster_ids}}, $orts_num );
		
		for( my $col_num = 4; $col_num <= $#cols; $col_num++ ){
			my $org_num  =  $col_num - 4;
			my $org_name =  $self->{org_name}[ $org_num ];
			
			my @genes = split ",", $cols[ $col_num ];
			
			foreach my $curr_gene ( @genes ){
				#print STDERR "Org num.: $org_num\n";
				#print STDERR "Org num.: $org_name\n";
				#getc();
				
				push( @{$self->{orts}[$orts_num]{$org_name}}, $curr_gene );
				$self->{gene_index}{$curr_gene} = $orts_num;
			}
		}
	}
	
	close(ORTS_FILE);
}


#828547707	Ecoli_H112180280	Ecoli_H112180280_POSTPRODIGAL_2	7000006964752101	7000006964752100	None	mannose-1-phosphate guanylyltransferase 1
#828547707	Ecoli_TY_2482_BGI	Ecoli_TY_2482_BGI_POSTPRODIGAL_2	7000006964764261	7000006964764260	None	mannose-1-phosphate guanylyltransferase 1
#828547707	EscCol_55989_GBD3	EscCol_55989_GBD3_POSTPRODIGAL_1	7000006964790339	7000006964790338	None	mannose-1-phosphate guanylyltransferase 1
#828547707	Esch_coli_04-8351_V1	Esch_coli_04-8351_V1_POSTPRODIGAL_1	7000006961080847	7000006961080846	EUDG_02241	mannose-1-phosphate guanylyltransferase 1

sub read_one_gene_per_line {
	my $self = shift;	
		
	open ORTS_FILE, $self->{filename} or die "Unable to open file " . $self->{filename} . "\n";
	while( <ORTS_FILE> ){
		my $line = $_;
		chomp $line;
		my ($orts_num, $org_name, $set, $transcript, $gene, $locus_name, $func_annot )  = split "\t", $line;
				
		$self->{cluster_ids_hash}{$orts_num} = 1;
				
		my $gene_identifier;
		if( $self->{gene_identifying_feature} eq "transcript" ){
			$gene_identifier = $transcript;
		}elsif( $self->{gene_identifying_feature} eq "gene" ){
			$gene_identifier = $gene;			
		}elsif( $self->{gene_identifying_feature} eq "locus_name" ){
			$gene_identifier = $locus_name;			
		}else{
			die "Unrecognizable gene identifying feature \'" . $self->{gene_identifying_feature} . "\'\n";
		}	
		
		#print STDERR "Orts num: $orts_num  Org name: >>>>>$org_name<<<<<    Identifier: $gene_identifier\n";
		#getc(); 
				
		push( @{ $self->{orts}[ $orts_num ]{ $org_name } }, $gene_identifier );
		$self->{gene_index}{$gene_identifier} = $orts_num;
		$self->{gene_annot}{$gene_identifier} = $func_annot;
		
		$self->{org_num}{$org_name} = 1;
	}	
	close(ORTS_FILE);
	
	# Associating a number to each organisms/genome
	my $cont_org = 0;
	foreach my $curr_org ( keys %{$self->{org_num}} ){
		$self->{org_num}{$curr_org} = $cont_org;
		$self->{org_name}[ $cont_org ]  = $curr_org;
		$cont_org++;
	}
	
	@{$self->{cluster_ids}} = sort {$a <=> $b} keys( %{$self->{cluster_ids_hash}} ); 
}

sub get_gene_annot{
	my $self = shift;
	my ($gene_name) = @_;
	return $self->{gene_annot}{$gene_name};
}

sub combined_annot_cluster{
	my $self = shift;	
	my ($orts_num) = @_;
	
	my %all_annot;		
	foreach my $org_name ( keys %{$self->{orts}[$orts_num]} ){
		foreach my $ortholog ( @{$self->{orts}[$orts_num]{$org_name}} ){
			my $curr_annot = $self->{gene_annot}{$ortholog};
			$curr_annot =~ s/name="//g;
			$curr_annot =~ s/"//g;
			$all_annot{$curr_annot} = 1;			
		}
	}	
	return join( ' && ', keys %all_annot ); 		
}	


sub combined_annot_from_genes{
	my $self = shift;
	my ($refArrGeneName) = @_;
	
	my %all_annot;
	foreach my $curr_gene ( @{$refArrGeneName} ){		
		my $curr_annot = $self->{gene_annot}{$curr_gene};
		$curr_annot =~ s/name="//g;
		$curr_annot =~ s/"//g;
		
		$all_annot{$curr_annot} = 1;
	}
	
	return join( ' && ', keys %all_annot ); 		
}	

sub get_orgs{
	my $self = shift;
	return @{$self->{org_name}};
}	

sub get_cluster_ids{
	my $self = shift;
	return @{$self->{cluster_ids}};
}	


sub get_org_name{
	my $self = shift;
	my ($org_num) = @_;
	
	return $self->{org_name}[ $org_num ];
}

sub get_org_num{
	my $self = shift;
	my ($org_name) = @_;
	
	return -1 if not defined $self->{org_num}{$org_name};	
	
	return $self->{org_num}{$org_name};
}

sub get_num_orgs{
	my $self = shift;
	return scalar( @{$self->{org_name}} );
}	

sub get_orts{
	my $self = shift;
	my ($gene_name, $org_name) = @_;
	
	my $orts_num = $self->{gene_index}{$gene_name}; 
	
	#print STDERR "Org name. $org_name Orts num:" . $orts_num . "\n";
	if( not defined( $orts_num ) ){
		print STDERR "Not able to find gene >>>$gene_name<<< in the ortholog clusters file.\n";
		return "";
	}elsif( not defined ( $self->{orts}[$orts_num]{$org_name} ) ){
		print STDERR "Not able to find orthologs of >>>$org_name<<< on cluster num. $orts_num\n";
		return "";		
	}elsif( not defined( $self->get_org_num( $org_name ) ) ){
		print STDERR "Genome >>>$org_name<<< not found in the ortholog clusters file.\n";
		return "";		
	}		
		 
	return 	@{$self->{orts}[$orts_num]{$org_name}};
}

sub get_core{
	my $self = shift;
	my ($gene_name, $org_name) = @_;
	
	my $orts_num = -1; 
	$orts_num = $self->{gene_index}{$gene_name}; 
	
	#print STDERR "Org name. $org_name Orts num:" . $orts_num . "\n";
	if( $orts_num == -1 || not defined ( $self->{orts}[$orts_num]{$org_name} ) ){
		print STDERR "Not able to find ortholog of gene $gene_name on $org_name\n";
		return 0;
	}
	 
	return 	@{$self->{orts}[$orts_num]{$org_name}};
}

sub print_cluster{
	my $self = shift;
	my ($gene_name) = @_;
	
	my $orts_num = $self->{gene_index}{$gene_name}; 
	
	#print STDERR "Org name. $org_name Orts num:" . $orts_num . "\n";
	if( not defined( $orts_num ) ){
		print STDERR "Not able to find gene >>>$gene_name<<< in the ortholog clusters file.\n";
	}
	
	foreach my $org_name ( keys %{$self->{orts}[$orts_num]} ){
		print "Org: $org_name\n";
		foreach my $ortholog ( @{$self->{orts}[$orts_num]{$org_name}} ){
			print "\t$ortholog\n";
		}
	}	
}
	
	

sub get_orts_by_cluster_num_org{
	my $self = shift;
	my ($orts_num, $org_name) = @_;
		
	#print STDERR "Org name. $org_name Orts num:" . $orts_num . "\n";
	if( not defined ( $self->{orts}[$orts_num]{$org_name} ) ){
		print STDERR "Not able to find ortholog on $org_name for cluster $orts_num\n";
		return 0;
	}
	 
	return 	@{$self->{orts}[$orts_num]{$org_name}};
}	



sub has_at_least_one_ort{
	my $self = shift;
	my ($gene_name, $orgs_list ) = @_;
	my @orgs =  split ":", $orgs_list;
		
	foreach my $curr_org ( @orgs ){
		return 1 if $self->get_orts( $gene_name, $curr_org ) != 0;
	}
	
	return 0;
}



return 1;