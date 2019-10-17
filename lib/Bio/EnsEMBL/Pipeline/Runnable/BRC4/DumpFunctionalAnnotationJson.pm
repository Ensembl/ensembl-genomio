package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFunctionalAnnotationJson;

use strict;
use warnings;

use JSON;
use File::Path qw(make_path);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub fetch_input {
  my ($self) = @_;

  $self->param( 'dba', $self->core_dba() );

  return;
}

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;

  print $self->param('species');
  my $output = $self->write_json();

  $self->dataflow_output_id({ "functional_annotation_json" => $output }, 2);
  return;
}

sub write_json {
  my ($self) = @_;

  my $sub_dir = $self->create_dir('json');
  $self->info(
          "Processing " . $self->production_name() . " into $sub_dir" );
  my $dba = $self->core_dba();

  # Get genes
  my @features;
  my $ga = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "gene" );

  my @items;
  push @items, @{$ga->fetch_all()};

  foreach my $item (@items) {

    my $type = '';
    if (ref($item) =~ /Gene/) {
      $type = 'gene';
    } elsif (ref($item) =~ /Transcript/) {
      $type = 'transcript';
    } elsif (ref($item) =~ /Translation/) {
      $type = 'translation';
    }
    my $syns = get_synonyms($item) // [];
    my ($xrefs, $onto) = get_xrefs($item);
    $xrefs //= [];
    $onto  //= [];

    my $feat = {
      id => $item->stable_id,
      object_type => $type,
      description => $item->description // '',
      synonyms => $syns,
      xrefs => $xrefs,
      ontology_terms => $onto,
    };

    push @features, $feat;
  }

  $dba->dbc()->disconnect_if_idle();

  # Write data to json
  $self->write_json_file($sub_dir, \@features);
} ## end sub write_json

sub write_json_file {
  my ($self, $sub_dir, $data) = @_;
  my $json_file_path =
    $sub_dir . '/' . $self->production_name() . '_functional_annotation.json';
  $self->info("Writing to $json_file_path");
  open my $json_file, '>', $json_file_path or
    die "Could not open $json_file_path for writing";
  print $json_file encode_json($data);
  close $json_file;
  $self->info("Write complete");
  return $json_file_path;
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
    my $xref = "$dbname:$id";

    if ($dbname =~ /[GS]O/) {
      push @onto, $xref;
    } else {
      push @xrefs, $xref;
    }
  }
  return \@xrefs, \@onto;
}

1;
