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
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'SpeciesFactory',
        'A->1' => 'Cleanup',
      },
    },
    {
      -logic_name => 'Cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
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
      -flow_into  => {
        '1->A' => ['CreateDB', 'GetData'],
        'A->1' => 'LoadData',
      },
    },

    {
      -logic_name => 'CreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'GetData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'CheckData',
      },
    },

    {
      -logic_name => 'CheckData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'LoadData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'PrepareAssemblyData',
        'A->1' => 'LoadMetadata',
      },
    },

    {
      -logic_name => 'PrepareAssemblyData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'LoadAssemblyData',
      },
    },

    {
      -logic_name => 'LoadAssemblyData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'SetupAssemblyMetadata',
      },
    },

    {
      -logic_name => 'SetupAssemblyMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'LoadMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'FillMetadata',
        'A->1' => 'ConstructRepeatLib',
      },
    },

    {
      -logic_name => 'FillMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'FillTaxonomy',
      },
    },

    {
      -logic_name => 'FillTaxonomy',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'ConstructRepeatLib',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
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
