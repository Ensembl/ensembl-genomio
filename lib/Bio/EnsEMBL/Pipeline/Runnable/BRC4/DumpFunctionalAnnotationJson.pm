package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFunctionalAnnotationJson;

use strict;
use warnings;

use JSON;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    metadata_name => 'functional_annotation',
  };
}

sub prepare_data {
  my ($self) = @_;

  my $dba = $self->core_dba();

  # Get genes
  my @features;
  my $ga = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "gene" );
  my $ta = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "transcript" );
  my $pa = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "translation" );

  my @items;
  push @items, @{$ga->fetch_all()};
  push @items, @{$ta->fetch_all()};
  push @items, @{$pa->fetch_all()};

  foreach my $item (@items) {

    my $type = '';
    if (ref($item) =~ /Gene/) {
      $type = 'gene';
    } elsif (ref($item) =~ /Transcript/) {
      $type = 'transcript';
    } elsif (ref($item) =~ /Translation/) {
      $type = 'translation';
    }

    # Basic metadata
    my %feat = (
      id => $item->stable_id,
      object_type => $type,
    );
    
    # Gene specific metadata
    if ($type eq 'gene') {
      my $syns = get_synonyms($item);
      $feat{synonyms} = $syns if $syns and @$syns;
      $feat{description} = $item->description if $item->description;
      $feat{is_pseudogene} = JSON::true if $item->biotype eq 'pseudogene';
    }

    # Xrefs (if any)
    my ($xrefs, $onto) = get_xrefs($item);
    $feat{xrefs} = $xrefs if $xrefs and @$xrefs;
    $feat{ontology_terms} = $onto if $onto and @$onto;

    push @features, \%feat;
  }

  $dba->dbc()->disconnect_if_idle();

  return \@features;
}

sub get_synonyms {
  my ($gene) = @_;

  my $disp = $gene->display_xref();
  return if not $disp;

  my $name = $disp->display_id;
  my @syns;
  push @syns, { synonym => $name, default => JSON::true } if $name;

  for my $syn (@{ $disp->get_all_synonyms() }) {
    push @syns, $syn;
  }

  return \@syns;
}

sub get_xrefs {
  my ($gene) = @_;

  my $entries = $gene->get_all_DBEntries();

  my @xrefs;
  my @onto;
  ENTRY: for my $entry (@$entries) {
    my $dbname = $entry->dbname;
    my $id = $entry->display_id;

    if ($dbname =~ /[GS]O/) {
      push @onto, $id;
    } else {
      my $xref = { dbname => $dbname, id => $id };
      push @xrefs, $xref;
    }
  }
  return \@xrefs, \@onto;
}

1;
