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
    location_so_term => [
            "apicoplast_chromosome",
            "chloroplast_chromosome",
            "chromoplast_chromosome",
            "cyanelle_chromosome",
            "leucoplast_chromosome",
            "macronuclear_chromosome",
            "micronuclear_chromosome",
            "mitochondrial_chromosome",
            "nuclear_chromosome"
          ]
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
  my $kba = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "KaryotypeBand" );

  # Get location SO terms
  my %location_so_term = map { $_ => 1 } @{$self->param("location_so_term")};

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
      if (@$syns) {
        $seq_region->{synonyms} = [ sort { $a->{name} cmp $b->{name} } map { { name => $_->name, source => $_->dbname } } @$syns ];
      }

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
      if ($so_term) {
        # Special case: location SO terms
        if ($location_so_term{$so_term->value}) {
          $seq_region->{location} = $so_term->value();
        } else {
          die("No SO_term allowed: '" . $so_term->value . "'");
          #$seq_region->{SO_term} = $so_term->value();
        }
      }

      # karyotype bands
      my $karyo_bands = $self->get_karyotype_bands($slice, $kba);
      $seq_region->{karyotype_bands} = $karyo_bands if (scalar(@$karyo_bands) > 0);

      push @seq_regions, $seq_region;
    }
  }

  $dba->dbc()->disconnect_if_idle();

  print(scalar(@seq_regions) . " seq_regions\n");

  # Sort
  @seq_regions = sort { $a->{name} cmp $b->{name} } @seq_regions;

  # Write data to json
  return \@seq_regions;
}

sub get_coords {
  my ($self, $dba) = @_;

  my $species = $self->param('species');

  my $coords_sql = "
    SELECT coord_system_id
    FROM coord_system
      LEFT JOIN meta USING(species_id)
    WHERE attrib LIKE '%default_version%'
      AND meta_key = 'species.production_name'
      AND meta_value = '$species';
    ";

  my $dbh = $dba->dbc->db_handle();

  my $array = $dbh->selectcol_arrayref($coords_sql);

  return @$array;
}

sub get_karyotype_bands {
  my ($self, $slice, $kba) = @_;

  #return [] if (!$slice->adaptor()->is_toplevel( $slice->get_seq_region_id() ));
  return [] if (!$slice->is_toplevel); # "At risk", can became legacy

  my $res = [];
  foreach my $band ( @{ $kba->fetch_all_by_Slice($slice) } ) {
    my $o = {
      start => $band->seq_region_start,
      end => $band->seq_region_end,
      name => $band->name,
    };
    my $stain = $band->stain;
    if (defined $stain) {
      $o->{stain} = $stain;
      $o->{structure} = "telomere" if ($stain eq 'TEL');
      $o->{structure} = "centromere" if ($stain eq 'ACEN');
    }
    push @$res, $o;
  }

  return $res;
}

1;
