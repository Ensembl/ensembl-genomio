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
  return if not $output_file;

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

  # Check that data is not empty
  if (ref($data) eq 'HASH' and scalar(%$data) == 0) {
    return;
  }
  elsif (ref($data) eq 'ARRAY' and scalar(@$data) == 0) {
    return;
  }
  elsif (ref($data) eq '' and not defined $data) {
    return;
  }

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
  print $jsonfh $json->pretty->canonical(1)->encode($data);
  close $jsonfh;
  
  $self->info("Write complete");
  return $json_file_path;
}

1;
