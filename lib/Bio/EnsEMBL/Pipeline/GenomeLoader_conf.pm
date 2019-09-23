package Bio::EnsEMBL::Pipeline::GenomeLoader_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);
use FindBin;

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{$self->SUPER::pipeline_wide_parameters},
        'debug' => $self->o('debug'),
        'tmp_dir' => $self->o('tmp_dir'),
    };
}

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    pipeline_name => 'genome_loader',
    email => $ENV{USER} . '@ebi.ac.uk',

    debug => 0,
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return
  [
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -flow_into  => {
        '1->A' => 'SpeciesFactory',
        'A->1' => 'Cleanup',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name        => 'SpeciesFactory',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        => {
        cmd => "ls",
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },
  ];
}

sub resource_classes {

  my ($self) = @_;
  
  return {
    'default'           => { 'LOCAL' => '' },
    'normal'            => {'LSF' => '-q production-rh74 -M 1000 -R "rusage[mem=1000]"' },
    'bigmem'           => {'LSF' => '-q production-rh74 -M 4000 -R "rusage[mem=4000]"' },
    'biggermem'           => {'LSF' => '-q production-rh74 -M 32000 -R "rusage[mem=32000]"' },
  }
}

1;
