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

=head1 NAME

Base config pipeline for BRC/Metazoa.

=cut

package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_base_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf');


sub resource_classes {
  my ($self) = @_;

  my $data_queue = $self->o('datamove_queue_name');
  my $queue = $self->o('queue_name');
  my $short = "1:00:00";
  my $long = "24:00:00";
  my $week = "7-00:00:00";

  my %resources = (
    %{$self->SUPER::resource_classes},
    'small'    => $self->make_resource({"queue" => $queue, "memory" => 100, "time" => $short}),
  );

  # Additional names in the form xGB e.g. "2GB"
  # From 2GB to 64GB
  my @mems = (2, 4, 8, 12, 16, 32, 64, 128, 256);
  my $time = $long;

  for my $mem (@mems) {
    my $name = "${mem}GB";
    $resources{$name} = $self->make_resource({"queue" => $queue, "memory" => $mem * 1000, "time" => $time});
  }

  # same for weekly classes
  $time = $week;
  for my $mem (@mems) {
    my $name = "${mem}GB_week";
    $resources{$name} = $self->make_resource({"queue" => $queue, "memory" => $mem * 1000, "time" => $time});
  }

  return \%resources;
}

1;
