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
  my $ga = $dba->get_adaptor('Gene');
  my $ta = $dba->get_adaptor('Transcript');
  my $pa = $dba->get_adaptor('Translation');

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
      $feat{version} = $item->version if $item->version;
      $feat{is_pseudogene} = JSON::true if $item->biotype eq 'pseudogene';
    }

    # Transcript specific metadata
    if ($type eq 'transcript') {
      $feat{version} = $item->version if $item->version;
    }

    # Xrefs (if any)
    my $xrefs = get_xrefs($item);
    $feat{xrefs} = $xrefs if $xrefs and @$xrefs;

    push @features, \%feat;
  }

  $dba->dbc()->disconnect_if_idle();

  # Sort for easier file comparison
  @features = sort { $a->{id} cmp $b->{id} } @features;

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
  my ($feature) = @_;

  my $entries = $feature->get_all_DBEntries();

  my @xrefs;
  my %found_entries;
  ENTRY: for my $entry (@$entries) {
    push @xrefs, create_xref($entry);
    $found_entries{$entry->dbID} = 1;
  }

  # Check that the display_xref is among the xref,
  # add it to the xref otherwise
  if ($feature->can('display_xref')) {
    my $display_entry = $feature->display_xref;
    if ($display_entry and not $found_entries{$display_entry->dbID}) {
      push @xrefs, create_xref($display_entry);
    }
  }

  return \@xrefs;
}

sub create_xref {
  my ($entry) = @_;

  my $dbname = $entry->dbname;
  my $id = $entry->display_id;

  my $xref = { dbname => $dbname, id => $id };
  $xref->{description} = $entry->description if ($entry->description);
  $xref->{info_type} = $entry->info_type if ($entry->info_type);
  $xref->{info_text} = $entry->info_text if ($entry->info_text);
  return $xref;
}

1;
