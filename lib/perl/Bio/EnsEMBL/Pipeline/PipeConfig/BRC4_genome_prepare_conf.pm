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

Genome prepare pipeline for BRC/Metazoa.

See full documentation in docs/BRC4_genome_preapre_conf.md

=cut

package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_prepare_conf;

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
my $metazoa_script_dir = "$root_dir/scripts/gff_metaparser";

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },
    ## default LSF queue
    queue_name => 'standard',

    ############################################
    # Config to be set by the user
    # MZ/BRC4 release id
    db_prefix => "",

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_genome_prepare' . $self->o('pipeline_tag'),

    # Working directory
    pipeline_dir => 'genome_prepare',
    data_dir => $self->o('data_dir'),
    output_dir => $self->o('output_dir'),

    debug => 0,
    ensembl_mode => 0,
    parser_conf => catfile($metazoa_script_dir, "conf/gff_metaparser.conf"),
    parser_patch => catfile($metazoa_script_dir, "conf/gff_metaparser/brc4.patch"),

    # If the genes appear to be split (multiple parts), merge them as one region
    merge_split_genes => 0,
    
    # Do not include those seq_regions (apply to all genomes, this should be seldom used)
    exclude_seq_regions => [],
    
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
    'schemas'      => $self->o('schemas'),
    pipeline_dir   => $self->o('pipeline_dir'),
    exclude_seq_regions   => $self->o('exclude_seq_regions'),

    download_dir   => catdir($self->o('pipeline_dir'), "download", '#accession#'),
    work_dir       => catdir($self->o('pipeline_dir'), "process_files", "#accession#"),
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
      -flow_into  => 'Genome_factory',
    },

    {
      # Create a thread for each species = manifest file
      # Output:
      -logic_name        => 'Genome_factory',
      -module         => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters        => {
        data_dir => $self->o('data_dir'),
        inputcmd => "find #data_dir# -type f -name '*.json'",
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -max_retry_count => 0,
      -flow_into => {
        '2' => { 'Process_genome_metadata' => { genome_json => '#_0#' } },
      },
    },

    {
      -logic_name => 'Process_genome_metadata',
      -module     => 'ensembl.brc4.runnable.process_genome_data',
      -language    => 'python3',
      -parameters     => {
        json_path => '#genome_json#',
      },
      -max_retry_count => 0,
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into => {
        2 => "Check_genome_schema"
      }
    },

    { -logic_name     => 'Check_genome_schema',
      -module         => 'ensembl.brc4.runnable.schema_validator',
      -language => 'python3',
      -parameters     => {
        json_file => '#genome_json#',
        json_schema => '#schemas#',
        metadata_type => "genome",
        hash_key => "#metadata_type#",
      },
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => {
        1 => 'Download_assembly_data'
      },
    },

    {
      -logic_name     => 'Download_assembly_data',
      -module         => 'ensembl.brc4.runnable.download_assembly_data',
      -language => 'python3',
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -rc_name        => 'default',
      -flow_into  => {
        '2->A' => 'Process_data',
        'A->2' => 'Manifest_maker',
      },
    },

    {
      -logic_name => 'Process_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into  => [
        WHEN("#gff3_raw#", ['Ungzip_gff3', 'Process_fasta_pep']),
        'Process_seq_region',
        'Process_fasta_dna',
      ],
    },

    {
      -logic_name     => 'Ungzip_gff3',
      -module         => 'ensembl.brc4.runnable.ungzip',
      -language => 'python3',
      -parameters  => {
        input => "#gff3_raw#",
        output => "#work_dir#/flat.gff3",
        out_name => "gff3_flat",
      },
      -analysis_capacity => 5,
      -failed_job_tolerance => 100,
      -rc_name        => 'default',
      -flow_into  => { 2 => 'Process_gff3' },
    },

    {
      -logic_name    => "Process_gff3",
      -module     => 'ensembl.brc4.runnable.process_gff3',
      -language    => 'python3',
      -parameters  => {
        in_gff3 => "#gff3_flat#",
        merge_split_genes => $self->o('merge_split_genes'),
        validate_gene_id => $self->o('validate_gene_id'),
      },
      -max_retry_count => 0,
      -failed_job_tolerance => 100,
      -analysis_capacity   => 5,
      -rc_name    => '8GB',
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

    {
      -logic_name => 'Process_seq_region',
      -module     => 'ensembl.brc4.runnable.process_seq_region',
      -language    => 'python3',
      -analysis_capacity   => 5,
      -rc_name    => 'default',
      -parameters     => {
        file_name => "seq_region",
      },
      -max_retry_count => 0,
      -failed_job_tolerance => 100,
      -flow_into  => {
        2 => 'Check_json_schema'
      },
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

    {
      -logic_name => 'Process_fasta_dna',
      -module     => 'ensembl.brc4.runnable.process_fasta',
      -language    => 'python3',
      -analysis_capacity   => 2,
      -rc_name    => 'default',
      -parameters     => {
        file_name => "fasta_dna",
      },
      -flow_into  => { 2 => '?accu_name=manifest_files&accu_address={file_name}&accu_input_variable=fasta_dna' },
    },

    {
      -logic_name => 'Process_fasta_pep',
      -module     => 'ensembl.brc4.runnable.process_fasta',
      -language    => 'python3',
      -analysis_capacity   => 2,
      -rc_name    => 'default',
      -parameters     => {
        file_name => "fasta_pep",
        in_genbank => '#gbff#',
        peptide => 1
      },
      -flow_into  => { 2 => '?accu_name=manifest_files&accu_address={file_name}&accu_input_variable=fasta_pep' },
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
      -rc_name         => '8GB',
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
