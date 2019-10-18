package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqRegionJson;

use strict;
use warnings;

use JSON;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    metadata_name => 'seq_region',
  };
}

sub prepare_data {
  my ($self) = @_;

  my $dba = $self->core_dba();

  # Get seq_regions
  my $sa = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "slice" );
  my $syna = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "seqregionsynonym" );

  # Get all top level seq regions
  my $slices = $sa->fetch_all('toplevel');

  my @seq_regions;
  foreach my $slice (@$slices) {
    my $syns = $syna->get_synonyms( $slice->get_seq_region_id() );
    my $seq_region = {
      name => $slice->seq_region_name(),
      coord_system_level => $slice->coord_system_name(),
      length => $slice->length(),
    };

    # Additional metadata
    # Synonyms? Array
    $seq_region->{synonyms} = [ map { $_->name } @$syns ] if @$syns;

    # Is circular? Boolean
    $seq_region->{circular} = JSON::true if $slice->is_circular;

    # alternate codon table? integer
    my ($codon_table) = @{$slice->get_all_Attributes('codon_table')};
    $seq_region->{codon_table} = int($codon_table->value()) if $codon_table;

    # SO_term? string
    my ($so_term) = @{$slice->get_all_Attributes('SO_term')};
    $seq_region->{SO_term} = $so_term->value() if $so_term;

    push @seq_regions, $seq_region;
  }

  $dba->dbc()->disconnect_if_idle();

  # Write data to json
  return \@seq_regions;
}

1;
