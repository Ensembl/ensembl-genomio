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
  my $sa = $dba->get_adaptor('Slice');
  my $csa = $dba->get_adaptor('CoordSystem');
  my $syna = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "seqregionsynonym" );

  # Get all coord system seq regions
  my @coord_ids = $self->get_coords($dba);

  my @seq_regions;
  for my $coord_id (@coord_ids) {
    my $coord_level_cs = $csa->fetch_by_dbID($coord_id);
    my $slices = $sa->fetch_all($coord_level_cs->name, undef, 1);

    print(scalar(@$slices) . " " . $coord_level_cs->name . "s (" . ($coord_level_cs->version || "no version") . ")\n");

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

      # non reference? Boolean
      my ($non_ref) = @{$slice->get_all_Attributes('non_ref')};
      $seq_region->{non_ref} = JSON::true if $non_ref;

      # SO_term? string
      my ($so_term) = @{$slice->get_all_Attributes('SO_term')};
      $seq_region->{SO_term} = $so_term->value() if $so_term;

      push @seq_regions, $seq_region;
    }
  }

  $dba->dbc()->disconnect_if_idle();

  print(scalar(@seq_regions) . " seq_regions\n");

  # Write data to json
  return \@seq_regions;
}

sub get_coords {
  my ($self, $dba) = @_;

  my $coords_sql = '
    SELECT coord_system_id
    FROM coord_system
    WHERE attrib LIKE "%default_version%";
  ';
  
  my $dbh = $dba->dbc->db_handle();

  my $array = $dbh->selectcol_arrayref($coords_sql);

  return @$array;
}

1;
