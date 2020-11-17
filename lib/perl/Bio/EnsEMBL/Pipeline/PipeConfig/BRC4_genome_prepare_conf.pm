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
        '2' => { 'Check_genome_schema' => { genome_json => '#_0#' } },
      },
    },

    {
      -logic_name     => 'Check_genome_schema',
      -module         => 'ensembl.brc4.runnable.schema_validator',
      -language => 'python3',
      -parameters     => {
        json_file => '#genome_json#',
        json_schema => '#schemas#',
        metadata_type => 'genome'
      },
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => {1 => { 'Read_genome_data' => INPUT_PLUS() } },
    },

    {
      -logic_name     => 'Read_genome_data',
      -module         => 'ensembl.brc4.runnable.read_json',
      -language => 'python3',
      -parameters     => {
        json_path => '#genome_json#',
        name => "genome_data"
      },
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => { 2 => 'Get_accession' },
    },

    {
      -logic_name     => 'Get_accession',
      -module         => 'ensembl.brc4.runnable.say_accession',
      -language => 'python3',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => { 2 => { 'Download_assembly_data' => INPUT_PLUS() } },
    },

    {
      -logic_name     => 'Download_assembly_data',
      -module         => 'ensembl.brc4.runnable.download_assembly_data',
      -language => 'python3',
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -rc_name        => 'default',
      -flow_into  => {
        '2->A' => { 'Process_data' => INPUT_PLUS() },
        'A->2' => { 'Manifest_maker' => INPUT_PLUS() },
      },
    },

    {
      -logic_name => 'Process_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into  => [
        WHEN("#gff3_raw#", ['Ungzip_gff3', 'Process_fasta_pep']),
        'Process_genome_metadata',
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

#    {
#      -logic_name    => "Process_gff3",
#      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#      -parameters  => {
#        gff3       => catfile("#work_dir#", "gene_models.gff3"),
#        functional_annotation => catfile("#work_dir#", "functional_annotation.json"),
#        seq_region_raw => catfile("#work_dir#", "seq_region_raw.json"),
#        gff3_name => "gff3",
#        func_name => "functional_anotation",
#        parser_conf => $self->o("parser_conf"),
#        parser_patch => $self->o("parser_patch") ? " --conf_patch " . $self->o("parser_patch") : "",
#
#        cmd => "python3 $metazoa_script_dir/gff3_meta_parse.py" .
#          " --dump_used_options" .
#          " --conf #parser_conf#" .
#          "#parser_patch#" .
#          " --gff_out #gff3#" . 
#          " --fann_out #functional_annotation#" .
#          " --seq_region_out #seq_region_raw#" .
#          " #gff3_flat#",
#      },
#      -failed_job_tolerance => 100,
#      -analysis_capacity   => 5,
#      -rc_name    => '8GB',
#      -meadow_type       => 'LSF',
#      -flow_into  => [
#          { 'GFF3_validation' => { gff3 => "#gff3#" } },
#          { "Check_json_schema" => { metadata_type => 'functional_annotation', metadata_json => '#functional_annotation#' } }
#        ],
#    },
    {
      -logic_name    => "Process_gff3",
      -module     => 'ensembl.brc4.runnable.process_gff3',
      -language    => 'python3',
      -parameters  => {
        in_gff3 => "#gff3_flat#",
        merge_split_genes => $self->o('merge_split_genes'),
      },
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
      -logic_name => 'Process_genome_metadata',
      -module     => 'ensembl.brc4.runnable.process_genome_data',
      -language    => 'python3',
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into => {
        2 => "Check_json_schema"
      }
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
      'default'  	=> {'LSF' => '-q production-rh74 -M 4000   -R "rusage[mem=4000]"'},
      '8GB'       => {'LSF' => '-q production-rh74 -M 8000   -R "rusage[mem=8000]"'},
      '15GB'      => {'LSF' => '-q production-rh74 -M 15000  -R "rusage[mem=15000]"'},
      '32GB'  	 	=> {'LSF' => '-q production-rh74 -M 32000  -R "rusage[mem=32000]"'},
      '64GB'  	 	=> {'LSF' => '-q production-rh74 -M 64000  -R "rusage[mem=64000]"'},
      '128GB'  	 	=> {'LSF' => '-q production-rh74 -M 128000 -R "rusage[mem=128000]"'},
      '256GB'  	 	=> {'LSF' => '-q production-rh74 -M 256000 -R "rusage[mem=256000]"'},
	}
}

1;
