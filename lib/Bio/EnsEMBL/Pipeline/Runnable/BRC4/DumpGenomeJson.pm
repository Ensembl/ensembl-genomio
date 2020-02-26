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
  my $species = $self->param_required("species");
  $self->info("Processing $species genome into $sub_dir");

  # Get meta table
  my $dba = $self->core_dba();
  my $ma = $dba->get_adaptor('MetaContainer');

  my $meta = {
    species => {},
    assembly => {},
    genebuild => {}
  };

  my %meta_list = (
    species => [ qw(taxonomy_id production_name scientific_name strain display_name division alias) ],
    assembly => [ qw(accession date name version) ],
    genebuild => [ qw(version method start_date) ],
    provider => [ qw(name url) ],
  );
  my %integer = map {$_ => 1} qw(species.taxonomy_id assembly.version);

  for my $domain (keys %meta_list) {
    for my $key (@{$meta_list{$domain}}) {
      my @values = @{ $ma->list_value_by_key($domain . '.' . $key) };
      next if scalar(@values) < 1;
      # Special case: assembly version should be a number,
      # extracted from the end of the assembly
      if ("$domain.$key" eq 'assembly.version') {
        die if @values > 1;
        my $av = $values[0];
        $av =~ s/^\D+(\d+)$/$1/;
        @values = ($av);
      }
      @values = map { int($_) } @values if $integer{$domain . '.' . $key};
      $meta->{$domain}->{$key} = (scalar(@values) > 1)? \@values : $values[0];
    }
  }

  # Check the assembly version
  $meta = $self->check_assembly_version($meta);

  $dba->dbc()->disconnect_if_idle();

  return $meta;
}

sub check_assembly_version {
  my ($self, $meta) = @_;

  my $assembly = $meta->{assembly};

  # Is there a version?
  my $version = $assembly->{version} // $assembly->{name};
  if ($version) {
    # Version is an integer
    if (not $version =~ /\D/) {
      $meta->{assembly}->{version} = int($meta->{assembly}->{version});
      return $meta;
    }

    # Version is not an integer, but ends in one
    if ($version =~ /[A-z]+([0-9]+)$/) {
      $meta->{assembly}->{version} = int($1);
      return $meta;
    
    # There is a version but I can't get a number out of it
    } else {
      die("Can't extract version number from assembly version: $version");
    }
  } else {
    use Data::Dumper;
    die("No assembly version found in meta: " . Dumper($meta));
  }
  return;
}

1;
