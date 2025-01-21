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


package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_base_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::ApiVersion qw(software_version);

use File::Basename;
use File::Spec::Functions qw(catdir catfile);
use FindBin;
use Class::Inspector;

my $package_path = Class::Inspector->loaded_filename(__PACKAGE__);
my $package_dir = dirname($package_path);
my $root_dir = "$package_dir/../../../../../..";

my $scripts_dir = "$root_dir/scripts";
my $config_dir = "$root_dir/config";

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    ## default LSF/Slurm queue name 
    #   is defined at Bio::EnsEMBL::EGPipeline::PipeConfig::EGGeneric_conf
    #   via Bio::EnsEMBL::EGPipeline::PrivateConfDetails::Impl
    #   (https://github.com/Ensembl/ensembl-production-imported/blob/main/lib/perl/Bio/EnsEMBL/EGPipeline/PrivateConfDetails/Impl.pm.example)
    #
    # queue_name => "standard",

    ############################################
    # Config to be set by the user
    # MZ/BRC4 release id
    release => $self->o('release'),
    db_prefix => "",

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_genome_loader' . $self->o('pipeline_tag'),
    email => $ENV{USER} . '@ebi.ac.uk',

    # Registry must contain the production db and taxonomy db
    # as well as the newly created dbs (e.g. via a prefix)
    registry      => $self->o('registry'),
    taxonomy_pass => '',

    # Working directory
    pipeline_dir => 'genome_loader_' . $self->o('release') . '_' . $self->o('ensembl_version'),

    # Meta configuration directory
    config_dir => $self->o('config_dir'),

    # Skip manifest checking (only if you know the checks are passed)
    check_manifest => 1,

    debug => 0,
    
    # External_db name map file
    external_db_map_name => 'default.txt',
    external_db_map => catfile($config_dir, 'external_db_map', $self->o('external_db_map_name')),

    # Do not load xrefs that we generate ourselves
    skip_ensembl_xrefs => 1,


    ############################################
    # Config unlikely to be changed by the user

    ensembl_version => software_version(),

    # Coordinate system order
    cs_order => 'ensembl_internal,chunk,contig,supercontig,non_ref_scaffold,scaffold,primary_assembly,superscaffold,linkage_group,chromosome',
    prune_agp => 0,
    unversion_scaffolds => 0,
    sr_syn_src_name  => 'ensembl_internal_synonym', # 50803
    division => 'EnsemblMetazoa',
    cs_tag_for_ordered => undef, # add cs_tag attr for sr from assembly.display_chromosome_order
    no_contig_ena_attrib => 0,

    # GenomeTools aliases (from eg-pipelines: Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf)
    gt_exe        => 'gt',
    gff3_tidy     => $self->o('gt_exe').' gff3 -tidy -sort -retainids',
    gff3_validate => $self->o('gt_exe').' gff3validator',

    # GFF3 cleaning params
    gff3_ignore_file => catfile(dirname(__FILE__), qw/gff3.ignore/),
    gff3_autoapply_manual_seq_edits => 1,

    # LoadGFF3 params
    load_pseudogene_with_CDS => 0,
    gff3_load_gene_source       => 'EnsemblMetazoa',
    gff3_load_logic_name        => 'brc4_import',
    #gff3_load_logic_name        => 'refseq_import_visible',
    gff3_load_analysis_module   => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::LoadGFF3',
    gff3_load_production_lookup => 1,
    # feature types
    gff3_use_polypeptides  => 0, # ignore 'polypeptides' lines !
    gff3_types_complete  => 1,

    # genes and transcripts versions
    default_feature_version => 1,
    no_feature_version_defaults => 0,

    # Enable brc4 features
    brc4_mode => 1,

    # ignore final stop codons check in Integrity
    ignore_final_stops => 0,

    # Change seq_region.name to the seq_region_synonym value from the mentioned source
    seq_name_code => "EBI_seq_region_name",

    # if loaded from RefSeq(GCF) change seq_region names to GenBank(GCA)
    swap_gcf_gca => 0,

    # default xref display_db
    xref_display_db_default => 'BRC4_Community_Annotation',
    xref_load_logic_name => 'brc4_import',

    # add_sequence mode (instead of creating db from scratch)
    add_sequence => 0,

    # run ProdDBsync parts before adding ad-hoc sequences (add_sequence  mode on)
    prod_db_sync_before_adding => 1,

    # default resource class name for Manifest_integrity
    manifest_integrity_rc_name => '8GB',

    # default resource class for 'LoadSequenceData' and 'AddSequence' steps
    load_sequence_data_rc_name => '8GB',

    # default resource class for LoadFunctionalAnnotation step
    load_func_ann_rc_name => '8GB',

    # default resource class for LoadGFF3 step
    load_gff3_rc_name => '16GB',

    # default resource class for CanonicalTranscriptsAttribs step
    canonical_transript_attribs_rc_name => 'default',

    # size of chunks to split contigs into (0 -- no splitting)
    sequence_data_chunk => 0,
    # coord system name for chunks
    chunk_cs_name => 'ensembl_internal',
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    debug          => $self->o('debug'),
    check_manifest => $self->o('check_manifest'),
    pipeline_dir   => $self->o('pipeline_dir'),
    ensembl_root_dir => $self->o('ensembl_root_dir'),

    taxonomy_url => $self->o('taxonomy_url'),
    taxonomy_pass => $self->o('taxonomy_pass'),
    dbsrv_url    => $self->o('dbsrv_url'),

    cs_order     => $self->o('cs_order'),
    prune_agp    => $self->o('prune_agp'),
    unversion_scaffolds => $self->o('unversion_scaffolds'),
    sr_syn_src_name  => $self->o('sr_syn_src_name'),
    division    => $self->o('division'),
    cs_tag_for_ordered => $self->o('cs_tag_for_ordered'),
    no_contig_ena_attrib => $self->o('no_contig_ena_attrib'),

    gt_exe        => $self->o('gt_exe'),
    gff3_tidy     => $self->o('gff3_tidy'),
    gff3_validate => $self->o('gff3_validate'),

    gff3_ignore_file            => $self->o('gff3_ignore_file'),
    gff3_autoapply_manual_seq_edits => $self->o('gff3_autoapply_manual_seq_edits'),

    gff3_load_gene_source       => $self->o('gff3_load_gene_source'),
    gff3_load_logic_name        => $self->o('gff3_load_logic_name'),
    gff3_load_analysis_module   => $self->o('gff3_load_analysis_module'),
    gff3_load_production_lookup => $self->o('gff3_load_production_lookup'),

    gff3_use_polypeptides => $self->o('gff3_use_polypeptides'),
    gff3_types_complete => $self->o('gff3_types_complete'),

    default_feature_version     => $self->o('default_feature_version'),
    no_feature_version_defaults => $self->o('no_feature_version_defaults'),

    load_pseudogene_with_CDS => $self->o('load_pseudogene_with_CDS'),
    brc4_mode => $self->o('brc4_mode'),
    ignore_final_stops => $self->o('ignore_final_stops'),

    swap_gcf_gca => $self->o('swap_gcf_gca'),

    xref_display_db_default => $self->o('xref_display_db_default'),
    xref_load_logic_name => $self->o('xref_load_logic_name'),

    external_db_map_name => $self->o('external_db_map_name'),
    external_db_map => $self->o('external_db_map'),

    add_sequence => $self->o('add_sequence'),
    prod_db_sync_before_adding => $self->o('prod_db_sync_before_adding'),

    manifest_integrity_rc_name => $self->o('manifest_integrity_rc_name'),
    load_sequence_data_rc_name => $self->o('load_sequence_data_rc_name'),
    load_func_ann_rc_name => $self->o('load_func_ann_rc_name'),
    load_gff3_rc_name => $self->o('load_gff3_rc_name'),
    canonical_transript_attribs_rc_name => $self->o('canonical_transript_attribs_rc_name'),

    sequence_data_chunk => $self->o('sequence_data_chunk'),
    chunk_cs_name        => $self->o('chunk_cs_name'),
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

# Ensures output parameters get propagated implicitly
sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

# Override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return 
      ' -reg_conf ' . $self->o('registry'),
  ;
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
      -flow_into  => {
        '1' => 'FillDBParams',
      },
    },
    {
      # fill db params from urls
      # output ..._host ..._port ...
      -logic_name => 'FillDBParams',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DbUrlToParams',
      -parameters => {
        db_urls => {
          dbsrv    => '#dbsrv_url#',
          taxonomy => '#taxonomy_url#',
        },
      },
      -rc_name    => 'default',
      -flow_into => {
        '1' => 'Manifest_factory',
      },
    },

    {
      # Create a thread for each species = manifest file
      # Output:
      -logic_name        => 'Manifest_factory',
      -module         => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters        => {
        data_dir => $self->o('data_dir'),
        inputcmd => "find -L #data_dir# -type f -name manifest.json",
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -max_retry_count => 0,
      -flow_into => {
        2 => WHEN('#check_manifest#', {
          'Manifest_check' => { manifest => '#_0#' }
        }, ELSE({
          'Prepare_genome' => {manifest => '#_0#' }
        })),
      },
    },

    # Checking files
    {
      -logic_name => 'Manifest_check',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -flow_into  => {
        '1->A' => 'json_schema_factory',
        'A->1' => 'Manifest_integrity',
      },
    },

    {
      -logic_name => 'json_schema_factory',
      -module     => 'ensembl.brc4.runnable.json_schema_factory',
      -language => 'python3',
      -rc_name    => 'default',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -flow_into  => {
        2 => 'check_json_schema',
      },
    },

    {
      -logic_name     => 'check_json_schema',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        # schemas_json_validate will automatically fetch the JSON schema file corresponding to metadata_type
        cmd => 'mkdir -p #log_path#; '
             . 'echo "checking #json# against #metadata_type#" > #log_path#/#metadata_type#.log; '
             . 'schemas_json_validate --json_file #json# --json_schema #metadata_type# '
             . '   >> #log_path#/#metadata_type#.log 2>&1 ',
        log_path => $self->o('pipeline_dir') . '/check_schemas',
        json => '#metadata_json#',
      },
      -analysis_capacity => 2,
      -failed_job_tolerance => 10, # in %
      -batch_size     => 50,
      -rc_name        => 'default',
    },

    {
      # Check the integrity of the manifest before loading anything
      -logic_name => 'Manifest_integrity',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        cmd => 'mkdir -p #log_path#; '
             . 'manifest_check_integrity --manifest_file #manifest# #ignore_final_stops_param#'
             . '   > #log_path#/check.log 2>&1 ',
        log_path => $self->o('pipeline_dir') . '/check_integrity',
        ignore_final_stops_param => '#expr( #ignore_final_stops# ? "--ignore_final_stops" : "")expr#',
      },
      -analysis_capacity   => 10,
      -rc_name         => $self->o('manifest_integrity_rc_name'),
      -max_retry_count => 0,
      -failed_job_tolerance => 10, # in %
      -flow_into => 'Prepare_genome',
    },

    {
      # Prepare all the metadata:
      # - species
      # - db_name
      # - manifest_metadata
      -logic_name => 'Prepare_genome',
      -module     => 'ensembl.brc4.runnable.prepare_genome',
      -language => 'python3',
      -parameters => {
        release => $self->o('release'),
        ensembl_version => $self->o('ensembl_version'),
        db_prefix => $self->o('db_prefix'),
      },
      -analysis_capacity   => 1,
      -max_retry_count => 0,
      -failed_job_tolerance => 100,
      -rc_name    => 'default',
      -flow_into  => {
        '2->A' => 'CreateDB',
        'A->2' => 'Finalize_database',
      },
    },

    {
      -logic_name => 'CreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -parameters => {
        has_gff3 => '#expr( #manifest_data#->{"gff3"} )expr#',
      },
      -analysis_capacity => 1,
      -batch_size     => 50,
      -flow_into  => {
        '1' => WHEN('#add_sequence#', 'ProdDbSyncAndAddSequence', ELSE('CleanUpAndCreateDB')),
      },
    },

    {
      -logic_name => 'ProdDbSyncAndAddSequence',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -flow_into  => {
        '1' => WHEN('#prod_db_sync_before_adding#', 'UpdateAnalysisDescription', ELSE('AddSequence')),
      },
    },

    {
      -logic_name        => 'UpdateAnalysisDescription',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
      -parameters => {
        species => '#expr( #genome_data#->{"species"}->{"production_name"} )expr#',
        group => 'core',
      },
      -analysis_capacity => 2,
      -rc_name    => 'default',
      -max_retry_count   => 3,
      -flow_into         => ['UpdateControlledTables'],
    },

    {
      -logic_name        => 'UpdateControlledTables',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateControlledTables',
      -parameters => {
        species => '#expr( #genome_data#->{"species"}->{"production_name"} )expr#',
        group => 'core',
      },
      -analysis_capacity => 2,
      -rc_name    => 'default',
      -max_retry_count   => 3,
      -flow_into         => ['AddSequence'],
    },

    {
      -logic_name => 'AddSequence',
      -module     => 'ensembl.brc4.runnable.load_sequence_data',
      -language => 'python3',
      -parameters        => {
        work_dir => $self->o('pipeline_dir') . '/#db_name#/add_sequence',
        load_additional_sequences => $self->o('add_sequence'),
        # N.B. chunking will work correctly only if it was used for initial loading
        sequence_data_chunk => $self->o('sequence_data_chunk'),
        chunk_cs_name        => $self->o('chunk_cs_name'),
      },
      -analysis_capacity   => 10,
      -rc_name         => $self->o('load_sequence_data_rc_name'),
      -max_retry_count => 0,
      -flow_into  => WHEN('#has_gff3#', 'Load_gene_models')
    },

    {
      # Init the Ensembl core
      -logic_name => 'CleanUpAndCreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url'),
        sql     => [
          'DROP DATABASE IF EXISTS #db_name#;' ,
          'CREATE DATABASE #db_name#;' ,
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into => 'LoadDBSchema',
    },

    {
      -logic_name => 'LoadDBSchema',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        input_file => $self->o('ensembl_root_dir') . '/ensembl/sql/table.sql',
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into  => 'FillMetadata',
    },

    {
      -logic_name => 'FillMetadata',
      -module     => 'ensembl.brc4.runnable.fill_metadata',
      -language => 'python3',
      -parameters        => {
        work_dir => $self->o('pipeline_dir') . '/#db_name#/fill_metadata',
        division => $self->o('division'),
        ignore => [ qw/ assembly.version / ],
        copy => { 'assembly.name' => 'assembly.default' },
      },
      -max_retry_count => 0,
      -analysis_capacity   => 10,
      -rc_name    => 'default',
      -flow_into  => 'PopulateControlledTables',
    },
    
    {
      -logic_name        => 'PopulateControlledTables',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateControlledTables',
      -parameters => {
        species => '#expr( #genome_data#->{"species"}->{"production_name"} )expr#',
        group => 'core',
      },
      -analysis_capacity => 2,
      -rc_name    => 'default',
      -max_retry_count   => 3,
      -flow_into         => ['LoadSequenceData']
    },

    {
      -logic_name => 'LoadSequenceData',
      -module     => 'ensembl.brc4.runnable.load_sequence_data',
      -language => 'python3',
      -parameters        => {
        work_dir => $self->o('pipeline_dir') . '/#db_name#/load_sequence',
        cs_order => $self->o('cs_order'),
        prune_agp => $self->o('prune_agp'),
        unversion_scaffolds => $self->o('unversion_scaffolds'),
        sr_syn_src  => $self->o('sr_syn_src_name'),
        external_db_map => $self->o('external_db_map'),
        cs_tag_for_ordered => $self->o('cs_tag_for_ordered'),
        no_contig_ena_attrib => $self->o('no_contig_ena_attrib'),
        swap_gcf_gca => $self->o('swap_gcf_gca'),
        sequence_data_chunk => $self->o('sequence_data_chunk'),
        chunk_cs_name        => $self->o('chunk_cs_name'),
      },
      -analysis_capacity   => 10,
      -rc_name         => $self->o('load_sequence_data_rc_name'),
      -max_retry_count => 0,
      -flow_into  => 'FillTaxonomy',
    },

    {
      -logic_name => 'FillTaxonomy',
      -module      => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
      -parameters  => {
        'taxon_id' => '#expr( #genome_data#->{"species"}->{"taxomy_id"} )expr#',
        'target_db' => {
          -host => "#dbsrv_host#",
          -port => "#dbsrv_port#",
          -pass => "#dbsrv_pass#",
          -user => "#dbsrv_user#",
          -dbname => "#db_name#",
        },
        'taxonomy_db' => {
          -host => "#taxonomy_host#",
          -port => "#taxonomy_port#",
          -user => "#taxonomy_user#",
          -pass => "#taxonomy_pass#",
          -dbname => "#taxonomy_dbname#",
        },
      },
      -rc_name    => 'default',
      -analysis_capacity   => 10,
      -flow_into  => WHEN('#has_gff3#' => 'Load_gene_models'),
    },

    {
      -logic_name => 'Load_gene_models',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
        'has_func_annotation' => '#expr( #manifest_data#->{"functional_annotation"} )expr#',
      },
      -rc_name    => 'default',
      -flow_into  => {
        '1->A' => 'GFF3CleanIgnored',
        'A->1' => WHEN('#has_func_annotation#', 'LoadFunctionalAnnotation', ELSE 'Finalize_gene_models'),
      },
    },

    {
      -logic_name => 'GFF3CleanIgnored',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . ' less #gff3_ignore_file# | '
            . '   perl -pe \'s/#.*$/\n/; s/^\s+//; s/ +$//; s/ +/ /g;\' | '
            . '   grep -vP \'^\s*$\' | '
            . '   perl -pe \'s/^/\t/; s/$/\t/\' > #gff3_ignore_pat#; '
            . ' less #gff3_orig_file# | grep -viFf #gff3_ignore_pat# > #gff3_clean_file# ',
        'log_path' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/clean',
        'gff3_orig_file' => '#expr( #manifest_data#->{"gff3"} )expr#',
        'gff3_clean_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/clean/clean.gff3',
        'gff3_ignore_file' => $self->o('gff3_ignore_file'),
        'gff3_ignore_pat' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/clean/gff3.ignore.pat',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 10,
      -flow_into => 'GFF3Tidy',
    },

    {
      -logic_name => 'GFF3Tidy',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . ' #gff3_tidy# #gff3_clean_file# > #gff3_tidy_file# 2> #log_path#/stderr ',
        'log_path' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/tidy',
        'gff3_tidy' => $self->o('gff3_tidy'),
        'gff3_clean_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/clean/clean.gff3',
        'gff3_tidy_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/tidy/tidy.gff3',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 10,
      -flow_into => [ 'GFF3Validate' ],
    },

    {
      -logic_name => 'GFF3Validate',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . ' #gff3_validate# #gff3_tidy_file# '
            . '   > #log_path#/stdout 2> #log_path#/stderr ',
        'log_path' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/validate',
        'gff3_validate' => $self->o('gff3_validate'),
        'gff3_tidy_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/tidy/tidy.gff3',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 10,
      -flow_into => [ 'DNAFastaGetTopLevel' ],
    },

    {
      -logic_name => 'DNAFastaGetTopLevel',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . ' cat  #gff3_tidy_file# | '
            . '   awk -F "\t" \'$0 !~ /^#/ && $3 != "seq_region" {print $1}\' | sort | uniq | '
            . '   perl #dumper# '
            . '     --host #dbsrv_host# --port #dbsrv_port# --user #dbsrv_user# --pass #dbsrv_pass# --dbname #db_name# '
            . '     > #dna_fasta_file# '
            . '     2> #log_path#/stderr ',
        'log_path' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/dna_fasta',
        'gff3_tidy_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/tidy/tidy.gff3',
        'dna_fasta_file' => $self->o('pipeline_dir') . '/#db_name#/load_gff3/dna_fasta/toplevel.fasta',
        'dumper' => "$scripts_dir/get_dna_fasta_for.pl",
      },
      -rc_name    => '4GB',
      -analysis_capacity   => 10,
      -flow_into => 'LoadGFF3AnalysisSetup',
    },

    {
      -logic_name => 'LoadGFF3AnalysisSetup',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AnalysisSetup',
      -parameters        => {
        logic_name         => $self->o('gff3_load_logic_name'),
        keep_logic_name    => "#add_sequence#",
        module             => $self->o('gff3_load_analysis_module'),
        db_url             => '#dbsrv_url#' . '#db_name#',
        production_lookup  => $self->o('gff3_load_production_lookup'),
        delete_existing    => "#expr(#add_sequence# ? 0 : 1)expr#",
      },
      -rc_name    => 'default',
      -analysis_capacity => 10,
      -max_retry_count   => 0,
      -flow_into => [ 'LoadGFF3' ],
    },

    {
      -logic_name => 'LoadGFF3',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::LoadGFF3',
      -parameters => {
        #species         => '#expr( #genome_data#->{"species"}->{"production_name"} )expr#',
        gff3_file       => $self->o('pipeline_dir') . '/#db_name#/load_gff3/tidy/tidy.gff3',
        fasta_file      => $self->o('pipeline_dir') . '/#db_name#/load_gff3/dna_fasta/toplevel.fasta',
        gene_source     => $self->o('gff3_load_gene_source'),
        logic_name      => $self->o('gff3_load_logic_name'),
        types_complete  => $self->o('gff3_types_complete'),
        polypeptides    => $self->o('gff3_use_polypeptides'), # it's better to ignore 'polypeptides' lines
        load_pseudogene_with_CDS => $self->o('load_pseudogene_with_CDS'),
        # dbparams
        db_url          => '#dbsrv_url#' . '#db_name#',
        # condition
        has_fasta_peptide => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
        # log
        log             => $self->o('pipeline_dir') . '/#db_name#/load_gff3/gff3loader.log',
      },
      -failed_job_tolerance => 10,
      -max_retry_count   => 0,
      -rc_name    => $self->o('load_gff3_rc_name'),
      -analysis_capacity   => 10,
      -flow_into => WHEN('#has_fasta_peptide#' => [ 'FixModels' ]),
    },

    {
      -logic_name => 'FixModels',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::FixModels',
      -parameters => {
        logic_name      => $self->o('gff3_load_logic_name'),
        db_url          => '#dbsrv_url#' . '#db_name#',
        protein_fasta_file      => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
        log             => $self->o('pipeline_dir') . '/#db_name#/load_gff3/fix_models.log',
      },
      -max_retry_count   => 0,
      -rc_name    => '16GB',
      -analysis_capacity   => 10,
      -flow_into => [ 'ApplySeqEdits' ],
    },

    {
      -logic_name => 'ApplySeqEdits',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ApplySeqEdits',
      -parameters => {
        logic_name      => $self->o('gff3_load_logic_name'),
        db_url          => '#dbsrv_url#' . '#db_name#',
        protein_fasta_file      => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
        log             => $self->o('pipeline_dir') . '/#db_name#/load_gff3/apply_seq_edits.log',
      },
      -max_retry_count   => 0,
      -rc_name    => '16GB',
      -analysis_capacity   => 10,
      -flow_into => [ 'ReportSeqEdits' ],
    },

    {
      -logic_name => 'ReportSeqEdits',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ReportSeqEdits',
      -parameters => {
        logic_name      => $self->o('gff3_load_logic_name'),
        db_url          => '#dbsrv_url#' . '#db_name#',
        protein_fasta_file      => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
          biotype_report_filename      => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/biotypes.txt',
          seq_edit_tt_report_filename  => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/seq_edit_tt.txt',
          seq_edit_tn_report_filename  => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/seq_edit_tn.txt',
          protein_seq_report_filename  => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/proteins.txt',
          protein_seq_fixes_filename   => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/proteins_fixes.txt',
      },
      -max_retry_count   => 0,
      -rc_name    => '16GB',
      -analysis_capacity   => 10,
      -flow_into => WHEN('#gff3_autoapply_manual_seq_edits#' => [ 'ApplyPatches' ]),
    },

    {
      -logic_name => 'ApplyPatches',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        input_file => $self->o('pipeline_dir') . '/#db_name#/load_gff3/reports/proteins_fixes.txt',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 10,
    },

    {
      -logic_name  => 'LoadFunctionalAnnotation',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . 'perl #fann_loader# '
            . '  --host #dbsrv_host# --port #dbsrv_port# --user #dbsrv_user# --pass #dbsrv_pass# --dbname #db_name# '
            . '  -json #fann_json_file# '
            . '  #default_feat_v# '
            . '  #default_db_display# '
            . '  -analysis_name #xref_load_logic_name# '
            . '  -external_db_map ' . $self->o('external_db_map')
            . '  -skip_ensembl_xrefs ' . $self->o('skip_ensembl_xrefs')
            . '  > #log_path#/stdout '
            . '  2> #log_path#/stderr ',
        'log_path'       => $self->o('pipeline_dir') . '/#db_name#/load_functional_annotation',
        'fann_loader'    => "$scripts_dir/load_fann.pl",
        'fann_json_file' => '#expr( #manifest_data#->{"functional_annotation"} )expr#',
        'default_feat_v' => '#expr( #no_feature_version_defaults# ? "": "-feature_version_default ".#default_feature_version# )expr#',
        'default_db_display' => '#expr( #xref_display_db_default# ? "-display_db_default ".#xref_display_db_default# : "" )expr#',
      },
      -max_retry_count   => 0,
      -rc_name    => $self->o('load_func_ann_rc_name'),
      -analysis_capacity   => 5,
      -flow_into => 'Finalize_gene_models',
    },

    # Finish up annotations
    {
      -logic_name => 'Finalize_gene_models',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -flow_into => WHEN('#no_feature_version_defaults#' => 'MetaCoord',
                    ELSE 'FixFeatureVersion'
      ),
    },

    {
      # Add versions for features if needed
      -logic_name => 'FixFeatureVersion',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#dbsrv_url#' . '#db_name#',
        sql     => [
          'UPDATE gene SET version = #default_feature_version# WHERE version IS NULL;',
          'UPDATE transcript SET version = #default_feature_version# WHERE version IS NULL;',
          'UPDATE translation SET version = #default_feature_version# WHERE version IS NULL;',
          'UPDATE exon SET version = #default_feature_version# WHERE version IS NULL;',
        ],
      },
      -rc_name    => 'default',
      -flow_into => 'MetaCoord',
    },

    {
      -logic_name    => "MetaCoord",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'base_dir'       => $self->o('ensembl_root_dir'),
        'cmd' => 'cd ' . $self->o('pipeline_dir') . ' ;'
            . ' perl #base_dir#/ensembl/misc-scripts/meta_coord/update_meta_coord.pl'
            . ' --dbhost #dbsrv_host#'
            . ' --dbport #dbsrv_port#'
            . ' --dbuser #dbsrv_user#'
            . ' --dbpass #dbsrv_pass#'
            . ' --dbpattern #db_name#'
            . ' ; cd -',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 1,
      -flow_into => 'Frameshifts',
    },

    {
      -logic_name    => "Frameshifts",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'base_dir'       => $self->o('ensembl_root_dir'),
        'cmd' => 'perl #base_dir#/ensembl/misc-scripts/frameshift_transcript_attribs.pl '
            . ' --dbhost #dbsrv_host# --dbport #dbsrv_port# --dbuser #dbsrv_user# --dbpass #dbsrv_pass# --dbpattern #db_name# ',
      },
      -rc_name    => 'default',
      -analysis_capacity   => 2,
      -flow_into => 'CanonicalTranscripts',
    },

    {
      -logic_name    => "CanonicalTranscripts",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #dump_path#; ' .
          ' perl #base_dir#/ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl '
            . ' --dbhost #dbsrv_host# --dbport #dbsrv_port# --dbuser #dbsrv_user# --dbpass #dbsrv_pass# --dbname #db_name# '
            . ' --write --coord_system_name toplevel '
            . ' --log #dump_path#/set_canonical_tr.log '
            . ' > #dump_path#/stdout 2> #dump_path#/stderr ',
        'base_dir'       => $self->o('ensembl_root_dir'),
        'dump_path' => $self->o('pipeline_dir') . '/#db_name#/canonical_transcripts',
      },
      -rc_name    => $self->o('canonical_transript_attribs_rc_name'),
      -analysis_capacity   => 2,
      -flow_into => 'CanonicalTranscriptsAttribs',
    },

    {
      # workaround to add `is_canonical` attributes for canonical transcrits
      -logic_name => 'CanonicalTranscriptsAttribs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        sql     => [
          'INSERT IGNORE INTO transcript_attrib (transcript_id, attrib_type_id, value) ' .
          '  SELECT g.canonical_transcript_id, at.attrib_type_id, 1 ' .
          '    FROM gene g, attrib_type at ' .
          '    WHERE at.code = "is_canonical"; ',
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
    },

    # Finalize database
    {
      -logic_name => 'Finalize_database',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -flow_into  => "PopulateAnalysis",
    },

    {
      -logic_name        => 'PopulateAnalysis',
      -module            => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
      -parameters => {
        species => '#expr( #genome_data#->{"species"}->{"production_name"} )expr#',
        group => 'core',
      },
      -analysis_capacity => 2,
      -rc_name    => 'default',
      -max_retry_count   => 0,
      -flow_into  => {
        '1->A' => 'UpdateSeqRegionNames',
        'A->1' => 'UpdateMetaData',
      },
    },

    {
      -logic_name => 'UpdateSeqRegionNames',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -max_retry_count   => 0,
      -flow_into => WHEN('#swap_gcf_gca#' => 'Swap_GCF_GCA',
                    ELSE 'Change_seq_region_name'
                    ),
    },
   
    {
      # SWAP GCF TO GCA seq_region names
      -logic_name => 'Swap_GCF_GCA',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        sql     => [
          # depends on the order of INSDC syns if there are few of them,
          #   which should't so for the new load
          'UPDATE seq_region sr ' .
          '  LEFT JOIN seq_region_synonym srs USING(seq_region_id) ' .
          '  LEFT JOIN external_db edb USING(external_db_id) ' .
          '  SET sr.name = srs.synonym ' .
          '  WHERE srs.external_db_id IS NOT NULL ' .
          '  AND edb.db_name = "INSDC"; ',
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -flow_into  => WHEN('#brc4_mode#', 'AddBRC4SeqRegionAttr'),
    },

    {
      # Use a different seq_region name
      -logic_name => 'Change_seq_region_name',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        seq_attrib_name => $self->o('seq_name_code'),
        sql     => [
          "UPDATE seq_region s ".
          "LEFT JOIN seq_region_attrib sa USING(seq_region_id) ".
          "LEFT JOIN attrib_type a USING(attrib_type_id) ".
          "SET s.name=value ".
          "WHERE sa.attrib_type_id IS NOT NULL ".
          "AND a.code = '#seq_attrib_name#';",
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
    },

    {
      # Add "(EBI|BRC4)_seq_region_name" seq_region_attrib(s) ("swap_gcf_gca" case)
      -logic_name => 'AddBRC4SeqRegionAttr',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        ebi_seq_region_name_attrib => 'EBI_seq_region_name',
        brc4_seq_region_name_attrib => 'BRC4_seq_region_name',
        sql     => [
          'INSERT IGNORE INTO seq_region_attrib (seq_region_id, attrib_type_id, value) ' .
          '  SELECT sr.seq_region_id, at.attrib_type_id, sr.name ' .
          '    FROM seq_region sr, attrib_type at ' .
          '    WHERE at.code in ("#ebi_seq_region_name_attrib#", "#brc4_seq_region_name_attrib#"); ',
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
    },

    {
      # Updating "added_seq.region_name" meta data
      -logic_name => 'UpdateMetaData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        ebi_seq_region_name_attrib => 'EBI_seq_region_name',
        brc4_seq_region_name_attrib => 'BRC4_seq_region_name',
        sql     => [
          # delete exisiting "added_seq.region_name"
          'DELETE FROM meta ' .
          '  WHERE meta_key = "added_seq.region_name" AND species_id = 1; ',

          # insert ignore "added_seq.region_name"
          'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) ' .
          '  SELECT 1, "added_seq.region_name", sr.name ' .
          '    FROM seq_region sr, seq_region_attrib sra, attrib_type at ' .
          '    WHERE sr.seq_region_id = sra.seq_region_id ' .
          '      AND sra.attrib_type_id = at.attrib_type_id ' .
          '      AND at.code = "added_seq_accession"; ',
        ],
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
    },

  ];
}

1;
