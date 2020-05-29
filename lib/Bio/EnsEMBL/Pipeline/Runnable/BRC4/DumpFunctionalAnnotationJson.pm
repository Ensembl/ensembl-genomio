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
  
  $self->load_external_db_map();

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
    $feat{version} = $item->version if $item->version;

    # Gene specific metadata
    if ($type eq 'gene') {
      my $syns = $self->get_synonyms($item);
      $feat{synonyms} = $syns if $syns and @$syns;
      $feat{is_pseudogene} = JSON::true if $item->biotype eq 'pseudogene';
      $feat{description} = $item->description if $item->description;
    }

    # Transcript specific metadata
    if ($type eq 'transcript') {
      $feat{description} = $item->description if $item->description;
    }


    # Xrefs (if any)
    my $xrefs = $self->get_xrefs($item);
    $feat{xrefs} = $xrefs if $xrefs and @$xrefs;

    push @features, \%feat;
  }

  $dba->dbc()->disconnect_if_idle();

  # Sort for easier file comparison
  @features = sort { $a->{id} cmp $b->{id} } @features;

  return \@features;
}

sub get_synonyms {
  my ($self, $gene) = @_;

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
  my ($self, $feature) = @_;

  my $entries = $feature->get_all_DBEntries();

  my @xrefs;
  my %found_entries;
  ENTRY: for my $entry (@$entries) {
    push @xrefs, $self->create_xref($entry);
    $found_entries{$entry->dbID} = 1;
  }

  # Check that the display_xref is among the xref,
  # add it to the xref otherwise
  if ($feature->can('display_xref')) {
    my $display_entry = $feature->display_xref;
    if ($display_entry and not $found_entries{$display_entry->dbID}) {
      push @xrefs, $self->create_xref($display_entry);
    }
  }

  return \@xrefs;
}

sub create_xref {
  my ($self, $entry) = @_;

  my $dbname = $entry->dbname;
  my $id = $entry->display_id;
  
  # Replace dbname from external_db map
  my $db_map = $self->param('db_map');
  if ($db_map and $db_map->{$dbname}) {
    $dbname = $db_map->{$dbname};
  }

  my $xref = { dbname => $dbname, id => $id };
  $xref->{description} = $entry->description if ($entry->description);
  $xref->{info_type} = $entry->info_type if ($entry->info_type and $entry->info_type ne 'NONE');
  $xref->{info_text} = $entry->info_text if ($entry->info_text);
  return $xref;
}

sub load_external_db_map {
  my ($self) = @_;
  
  my %map;
  my $map_path = $self->param("external_db_map");
  if ($map_path) {
    open my $mapfh, "<", $map_path or die "$!";
    while (my $line = readline $mapfh) {
      chomp $line;
      next if $line =~ /^\*$/ or $line =~ /^#/;
      # We use the mapping in reverse order because we dump
      my ($to, $from) = split("\t", $line);
      $map{$from} = $to;
    }
    close $mapfh;
  }
  $self->param('db_map', \%map);
}

1;
