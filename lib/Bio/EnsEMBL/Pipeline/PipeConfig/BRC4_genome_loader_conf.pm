package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf;

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

    ##############################
    # Config to be set by the user
    # VB release id
    release => undef,

    # Ensembl version (deduced from the environment?)
    version => undef,

    # Meta configuration directory
    meta_dir => undef,

    # Working directory
    pipeline_dir => undef,

    # Server configuration (where the databases will be created)
    db_host => undef,
    db_port => undef,
    db_user => undef,
    db_pass => undef,
    ##############################

    # Basic pipeline configuration
    pipeline_name => 'genome_loader',
    email => $ENV{USER} . '@ebi.ac.uk',

    debug => 0,
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return
  [
    # Starting point
    # Create the tmp pipeline directory
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -parameters        => {
        pipeline_dir => $self->o('pipeline_dir'),
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'SpeciesList',
        'A->1' => 'Cleanup',
      },
    },

    {
      # Delete the temp working directory
      -logic_name => 'Cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      # Create a thread for each species and define its main parameters
      # Reads all the json files in a meta directory and creates a thread for each valid one
      # Also creates a separate tmp_dir = work_dir/#species#
      # All parameters output become available to all downstream analyses
      #
      # Output:
      # species =#genus#_#species_name#_#GCA#
      # db =  #species#_core_#release#_#version#_#assembly#
      # meta (extracted from the species metadata file)
      # tmp_dir = #pipeline_dir#/#species#
      -logic_name        => 'SpeciesList',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters        => {
        meta_dir => $self->o('meta_dir'),
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => ['CreateDB', 'GetData'],
        'A->1' => { 'LoadData' => INPUT_PLUS() },
      },
    },

    {
      # Init the Ensembl core
      -logic_name => 'CreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      # Retrieve the sequences files
      # TODO: expand how we should get the data:
      #    cp from a local dir, from ftp?
      #    Do we get the files path from the meta conf file?
      #    How much can this be automated?
      #
      # Output:
      # - fasta file (only one?)
      # - AGP?;qa
      #
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
      # Check and sanitize the sequences files
      # - Unzip
      # - Remove ambiguous IUPAC nucleotides
      # - Check vs the list of seq_regions (should we get a separate list?)
      -logic_name => 'CheckData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [],
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      # Head analysis for the loading of data
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
      # Head analysis for the loading of the metadata
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
