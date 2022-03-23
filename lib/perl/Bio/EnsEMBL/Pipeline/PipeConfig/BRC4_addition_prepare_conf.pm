package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_addition_prepare_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::ApiVersion qw(software_version);

use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use FindBin;
use Class::Inspector;

my $package_path = Class::Inspector->loaded_filename(__PACKAGE__);
my $package_dir = dirname($package_path);
my $root_dir = "$package_dir/../../../../../..";

my $schema_dir = "$root_dir/schema";

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },
    ## default LSF queue
    queue_name => 'standard',

    ############################################
    # Config to be set by the user
    gb_accession => $self->o('gb_accession'),
    
    # Gene ids prefix, to make sure we don't just use gene names
    ids_prefix => "",

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_addition_prepare' . $self->o('pipeline_tag'),

    # Working directory
    pipeline_dir => 'addition_prepare',
    output_dir => $self->o('output_dir'),

    debug => 0,
    ensembl_mode => 0,

    # If the genes appear to be split (multiple parts), merge them as one region
    merge_split_genes => 0,
    
    # Enforce a strong gene ID pattern (replace by GeneID if available)
    validate_gene_id => 0,

    ############################################
    # Config unlikely to be changed by the user

    ## Metadata parameters
    'schemas' => {
      'seq_region' => catfile($schema_dir, "seq_region_schema.json"),
      'seq_attrib' => catfile($schema_dir, "seq_attrib_schema.json"),
      'functional_annotation' => catfile($schema_dir, "functional_annotation_schema.json"),
      'genome' => catfile($schema_dir, "genome_schema.json"),
      'manifest' => catfile($schema_dir, "manifest_schema.json"),
    },

    ## gff3 parameters
    'gt_exe'          => 'gt',
    'gff3_tidy'       => $self->o('gt_exe').' gff3 -tidy -sort -retainids -force',
    'gff3_validate'   => $self->o('gt_exe').' gff3validator',
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    debug          => $self->o('debug'),
    schemas      => $self->o('schemas'),
    pipeline_dir   => $self->o('pipeline_dir'),
    work_dir => catdir($self->o('pipeline_dir'), "process_files", "#gb_accession#"),
  };
}
sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('pipeline_dir'),
    ];
}

sub pipeline_analyses {
  my ($self) = @_;

  return
  [
    # Starting point
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'Download_genbank',
        'A->1' => 'Manifest_maker',
      },
    },

    {
      -logic_name => 'Download_genbank',
      -module     => 'ensembl.brc4.runnable.download_genbank',
      -language    => 'python3',
      -parameters     => {
        gb_accession => $self->o('gb_accession'),
        download_dir => catdir($self->o('pipeline_dir'), "download"),
      },
      -max_retry_count => 0,
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into => {
        2 => "Extract_from_gb"
      }
    },

    {
      -logic_name     => 'Extract_from_gb',
      -module         => 'ensembl.brc4.runnable.extract_from_gb',
      -language => 'python3',
      -parameters     => {
        gb_file => '#gb_file#',
        ids_prefix => $self->o('ids_prefix'),
        pepkey => 'fasta_pep',
        dnakey => 'fasta_dna',
      },
      -analysis_capacity => 1,
      -failed_job_tolerance => 0,
      -rc_name        => 'default',
      -flow_into  => {
        2 => ['Process_gff3',
          { 'Check_json_schema' => { 'metadata_type' => 'genome', 'metadata_json' => '#genome_data#' } },
          { 'Check_json_schema' => { 'metadata_type' => 'seq_region', 'metadata_json' => '#seq_region#' } },
          '?accu_name=manifest_files&accu_address={pepkey}&accu_input_variable=fasta_pep',
          '?accu_name=manifest_files&accu_address={dnakey}&accu_input_variable=fasta_dna',
        ],
      }
    },

    {
      -logic_name    => "Process_gff3",
      -module     => 'ensembl.brc4.runnable.process_gff3',
      -language    => 'python3',
      -parameters  => {
        in_gff3 => "#gff3#",
        merge_split_genes => $self->o('merge_split_genes'),
        validate_gene_id => $self->o('validate_gene_id'),
      },
      -max_retry_count => 0,
      -failed_job_tolerance => 100,
      -analysis_capacity   => 5,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
          2 => 'GFF3_validation',
          3 => "Check_json_schema",
        },
    },
   
   { -logic_name     => 'GFF3_validation',
     -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters     => {
       temp_gff3 => "#gff3#" . ".tmp",
       cmd => "mv #gff3# #temp_gff3#" .
       " && " . $self->o('gff3_tidy') . " -o #gff3# #temp_gff3#" .
       " && " . $self->o('gff3_validate') . ' #gff3#',
       file_name => "gff3",
     },
      -max_retry_count => 0,
     -analysis_capacity => 10,
     -batch_size        => 10,
     -rc_name           => 'default',
     -flow_into  => '?accu_name=manifest_files&accu_address={file_name}&accu_input_variable=gff3',
   },

    { -logic_name     => 'Check_json_schema',
      -module         => 'ensembl.brc4.runnable.schema_validator',
      -language => 'python3',
      -parameters     => {
        json_file => '#metadata_json#',
        json_schema => '#schemas#',
        hash_key => "#metadata_type#",
      },
      -analysis_capacity => 2,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => { 1 => '?accu_name=manifest_files&accu_address={hash_key}&accu_input_variable=metadata_json' },
    },

    # Collate files to their final dir
    { -logic_name  => 'Manifest_maker',
      -module      => 'ensembl.brc4.runnable.manifest',
      -language    => 'python3',
      -max_retry_count => 0,
      -analysis_capacity   => 5,
      -rc_name         => 'default',
      -parameters     => {
        output_dir => $self->o('output_dir'),
        genome_name => $self->o('gb_accession'),
        species => 'species',
      },
      -flow_into       => { '2' => 'Manifest_check' },
    },

    {
      -logic_name     => 'Manifest_check',
      -module         => 'ensembl.brc4.runnable.schema_validator',
      -language => 'python3',
      -parameters     => {
        metadata_type => 'manifest',
        json_file => '#manifest#',
        json_schema => '#schemas#',
      },
      -analysis_capacity => 2,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into       => ['Integrity_check', "Manifest_stats"],
    },

    {
      -logic_name  => 'Integrity_check',
      -module      => 'ensembl.brc4.runnable.integrity',
      -language    => 'python3',
      -parameters     => {
        ensembl_mode => $self->o('ensembl_mode'),
      },
      -failed_job_tolerance => 100,
      -analysis_capacity   => 5,
      -rc_name         => 'default',
      -max_retry_count => 0,
#      -flow_into       => 'Manifest_stats',
    },
    {
      -logic_name  => 'Manifest_stats',
      -module      => 'ensembl.brc4.runnable.manifest_stats',
      -language    => 'python3',
      -analysis_capacity   => 1,
      -rc_name         => 'default',
      -max_retry_count => 0,
    },
  ];
}

sub resource_classes {
    my $self = shift;
    return {
      'default' => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 4000   -R "rusage[mem=4000]"'},
      '8GB'     => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 8000   -R "rusage[mem=8000]"'},
      '15GB'    => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 15000  -R "rusage[mem=15000]"'},
      '32GB '   => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 32000  -R "rusage[mem=32000]"'},
      '64GB'    => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 64000  -R "rusage[mem=64000]"'},
      '128GB'   => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 128000 -R "rusage[mem=128000]"'},
      '256GB  ' => {'LSF' => '-q ' . $self->o("queue_name") . ' -M 256000 -R "rusage[mem=256000]"'},
    }
}

1;
