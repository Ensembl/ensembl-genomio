package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase;

use strict;
use warnings;

use JSON;
use File::Path qw(make_path);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub param_defaults {
  return {
    metadata_name => 'metadata',
  };
}

sub run {
  my ($self) = @_;

  my $metadata_name =  $self->param('metadata_name');
  my $data = $self->prepare_data();

  # Write data to json
  my $output_file = $self->write_json_file($data);

  $self->dataflow_output_id(
    {
      "metadata_json" => $output_file,
      "metadata_type" => $metadata_name
    }, 2);
  return;
}

sub prepare_data {
  my ($self) = @_;
  return {}
}

sub write_json_file {
  my ($self, $data) = @_;

  my $metadata_name =  $self->param('metadata_name');
  my $sub_dir = $self->create_dir('json');
  $self->info("Processing " . $self->production_name() . " into $sub_dir" );
  my $json_file_path =
    $sub_dir . '/' . $self->production_name() . '_' . $metadata_name . '.json';
  $self->info("Writing to $json_file_path");
  open my $json_file, '>', $json_file_path or
    die "Could not open $json_file_path for writing";
  print $json_file encode_json($data);
  close $json_file;
  $self->info("Write complete");
  return $json_file_path;
}

1;
