package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase;

use strict;
use warnings;
use autodie;

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
  my $species = $self->param_required("species");
  $self->info("Processing $species into $sub_dir" );
  my $json_file_path =
    $sub_dir . '/' . $species . '_' . $metadata_name . '.json';
  $self->info("Writing to $json_file_path");
  
  # Print pretty JSON
  my $json = JSON->new;
  open my $jsonfh, '>', $json_file_path;
  print $jsonfh $json->pretty->encode($data);
  close $jsonfh;
  
  $self->info("Write complete");
  return $json_file_path;
}

1;
