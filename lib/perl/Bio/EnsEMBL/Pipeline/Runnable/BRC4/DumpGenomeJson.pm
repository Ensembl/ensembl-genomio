=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


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
    assembly => {}
  };

  my %meta_list = (
    species => [ qw(taxonomy_id production_name scientific_name strain display_name division alias annotation_source) ],
    assembly => [ qw(accession date name version provider_name provider_url) ],
    genebuild => [ qw(version method start_date method_display) ],
    annotation => [ qw(provider_name provider_url) ],
    BRC4 => [ qw(organism_abbrev component) ],
    added_seq => [ qw(region_name) ],
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
    if ($version =~ /[A-z_\.]+([0-9]+)$/) {
      $meta->{assembly}->{version} = int($1);
      return $meta;
    
    # There is a version but I can't get a number out of it
    } else {
      # Last resort: try to get the version from the assembly accession
      my $acc = $assembly->{accession};
      if ($acc and $acc =~ /\.(\d+$)/) {
        my $version = $1;
        $meta->{assembly}->{version} = $version;
        return $meta;
      } else {
        die("Can't extract version number from assembly version: $version (accession = $acc)");
      }
    }
  } else {
    use Data::Dumper;
    die("No assembly version found in meta: " . Dumper($meta));
  }
  return;
}

1;
