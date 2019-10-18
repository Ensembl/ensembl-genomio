package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeJson;

use strict;
use warnings;

use JSON;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    metadata_name => 'genome',
  };
}

sub prepare_data {
  my ($self) = @_;

  my $sub_dir = $self->create_dir('json');
  $self->info(
          "Processing " . $self->production_name() . " genome into $sub_dir" );

  # Get meta table
  my $ma = Bio::EnsEMBL::Registry->get_adaptor($self->production_name, 'Core', 'MetaContainer');
  my $dba = $self->core_dba();

  my $meta = {
    species => {},
    assembly => {},
    genebuild => {}
  };

  my %meta_list = (
    species => [ qw(taxonomy_id production_name scientific_name strain) ],
    assembly => [ qw(accession date name version) ],
    genebuild => [ qw(version) ],
    provider => [ qw(name url) ],
  );
  for my $domain (keys %meta_list) {
    for my $key (@{$meta_list{$domain}}) {
      my @values = @{ $ma->list_value_by_key($domain . '.' . $key) };
      my ($value) = @values;
      next if not defined $value;
      $value = int($value) if $value =~ /^\d+$/;
      $meta->{$domain}->{$key} = $value;
    }
  }

  $dba->dbc()->disconnect_if_idle();

  return $meta;
}

1;
