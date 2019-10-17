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

  foreach my $gene (@{$ga->fetch_all()}) {
    my $feat = {
      id => $gene->stable_id,
      object_type => 'gene',
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

1;
