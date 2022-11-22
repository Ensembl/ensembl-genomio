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

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub make_resource {
  my ($self, $conf) = @_;

  my $res = {
    'LSF' => _lsf_resource($conf),
    'SLURM' => _slurm_resource($conf),
  };
  return $res;
}

sub _lsf_resource {
  my ($conf) = @_;

  my $mem = $conf->{memory};
  my $cpus = $conf->{cpus};
  my $queue = $conf->{queue};
  my $time = $conf->{time};

  if ($time and $time =~ /(\d+):(\d\d):(\d\d)/) {
    $time = "$1:$2";
  }
  
  my @res_params = ("-M $mem", "-R 'rusage[mem=$mem]'");
  push @res_params, "-q $queue" if $queue;
  push @res_params, "-We $time" if $time;
  push @res_params, "-n $cpus" if $cpus;
  push @res_params, $conf->{lsf_params} if $conf->{lsf_params};
  my $res_string = join(" ", @res_params);
  return $res_string;
}

sub _slurm_resource {
  my ($conf) = @_;

  my $mem = $conf->{memory};
  my $cpus = $conf->{cpus};
  my $queue = $conf->{queue};
  my $time = $conf->{time};

  # Prepare memory string
  die("Memory needed for resource") if not $mem;
  my $rmem;
  if ($mem > 1000) {
    $mem = int($mem/1000);
    $rmem = $mem . 'g';
  } else {
    $rmem = $mem . "m";
  }
  
  my @res_params = ("--mem=$rmem");
  push @res_params, "--time=$time" if $time;
  push @res_params, "-c $cpus" if $cpus;
  push @res_params, "--partition=$queue" if $queue;
  push @res_params, $conf->{slurm_params} if $conf->{slurm_params};
  my $res_string = join(" ", @res_params);
  return $res_string;
}

sub resource_classes {
  my ($self) = @_;

  my $data_queue = $self->o('datamove_queue_name');
  my $queue = $self->o('queue_name');
  my $short = "1:00:00";
  my $long = "24:00:00";

  my %resources = (
    'default'  => $self->make_resource({"queue" => $queue, "memory" => 4_000, "time" => $short}),
    'normal'   => $self->make_resource({"queue" => $queue, "memory" => 4_000, "time" => $long}),
    'small'    => $self->make_resource({"queue" => $queue, "memory" => 100, "time" => $short}),
    'datamove' => $self->make_resource({"queue" => $data_queue, "memory" => 100, "time" => $short}),
  );

  # Additional names in the form xGB e.g. "2GB"
  # From 2GB to 64GB
  my @mems = (2, 4, 8, 12, 16, 32, 64);
  my $time = $long;

  for my $mem (@mems) {
    my $name = "${mem}GB";
    $resources{$name} = $self->make_resource({"queue" => $queue, "memory" => $mem * 1000, "time" => $time});
  }

  return \%resources;
}

1;
