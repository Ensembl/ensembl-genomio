package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_loader_conf;

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

my $scripts_dir = "$root_dir/scripts";
my $schema_dir = "$root_dir/schema";
my $data_dir = "$root_dir/data";

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    ############################################
    # Config to be set by the user
    # MZ/BRC4 release id
    release => $self->o('release'),
    db_prefix => "",

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_genome_loader_' .
      $self->o('release') . '_' . $self->o('ensembl_version') . $self->o('pipeline_tag'),
    email => $ENV{USER} . '@ebi.ac.uk',

    # Working directory
    pipeline_dir => 'genome_loader_' . $self->o('release') . '_' . $self->o('ensembl_version'),

    # Meta configuration directory
    data_dir => $self->o('data_dir'),

    # Skip manifest checking (only if you know the checks are passed)
    check_manifest => 1,

    debug => 0,

    ## Metadata parameters
    'schemas' => {
      'seq_region' => catfile($schema_dir, "seq_region_schema.json"),
      'seq_attrib' => catfile($schema_dir, "seq_attrib_schema.json"),
      'functional_annotation' => catfile($schema_dir, "functional_annotation_schema.json"),
      'genome' => catfile($schema_dir, "genome_schema.json"),
      'manifest' => catfile($schema_dir, "manifest_schema.json"),
    },
    
    # External_db name map file
    external_db_map => catfile($data_dir, 'external_db_map_default.txt'),

    ############################################
    # Config unlikely to be changed by the user

    ensembl_version => software_version(),

    # Coordinate system order
    cs_order => 'chunk,contig,supercontig,non_ref_scaffold,scaffold,primary_assembly,superscaffold,linkage_group,chromosome',
    prune_agp => 0,
    unversion_scaffolds => 0,
    sr_syn_src_name  => 'ensembl_internal_synonym', # 50803
    division => 'EnsemblMetazoa',
    cs_tag_for_ordered => undef, # add cs_tag attr for sr from assembly.display_chromosome_order

    # GenomeTools aliases (from eg-pipelines: Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf)
    gt_exe        => 'gt',
    gff3_tidy     => $self->o('gt_exe').' gff3 -tidy -sort -retainids',
    gff3_validate => $self->o('gt_exe').' gff3validator',

    # GFF3 cleaning params
    gff3_ignore_file => catfile(dirname(__FILE__), qw/gff3.ignore/),
    gff3_autoapply_manual_seq_edits => 1,

    # LoadGFF3 params
    load_pseudogene_with_CDS => 1,
    gff3_load_gene_source       => 'EnsemblMetazoa',
    gff3_load_logic_name        => 'gff3_genes',
    #gff3_load_logic_name        => 'refseq_import_visible',
    gff3_load_analysis_module   => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::LoadGFF3',
    gff3_load_production_lookup => 1,
    # feature types
    gff3_use_polypeptides  => 0, # ignore 'polypeptides' lines !
    gff3_gene_types   => [ qw/
        gene pseudogene miRNA_gene ncRNA_gene
        rRNA_gene snoRNA_gene snRNA_gene tRNA_gene
      /],
    gff3_mrna_types   => [ qw/
        mRNA transcript pseudogenic_transcript
        pseudogenic_rRNA pseudogenic_tRNA
        ncRNA lincRNA lncRNA miRNA pre_miRNA
        RNase_MRP_RNA RNAse_P_RNA rRNA snoRNA
        snRNA sRNA SRP_RNA tRNA scRNA
        lnc_RNA
      /],
    gff3_exon_types   => [qw/ exon pseudogenic_exon /],
    gff3_cds_types    => [qw/ CDS /],
    gff3_utr_types    => [qw/ five_prime_UTR three_prime_UTR /],
    gff3_ignore_types => [ qw/
        misc_RNA RNA
        match match_part
        sequence_feature
        cDNA_match nucleotide_match protein_match
        polypeptide protein
        chromosome supercontig contig
        region biological_region
        regulatory_region repeat_region
        golden_path_region intron orthologous_to
      /],
    gff3_types_complete  => 1,

    # genes and transcripts versions
    default_feature_version => 1,
    no_feature_version_defaults => 0,

    # disable brc4 features
    no_brc4_stuff => 0,

    # ignore final stop codons check in Integrity
    ignore_final_stops => 0,

    # Rename seq_region name in the seq_region table with this attribute
    seq_name_code => "EBI_seq_region_name"
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    debug          => $self->o('debug'),
    check_manifest => $self->o('check_manifest'),
    'schemas'      => $self->o('schemas'),
    pipeline_dir   => $self->o('pipeline_dir'),
    ensembl_root_dir => $self->o('ensembl_root_dir'),

    proddb_url   => $self->o('proddb_url'),
    taxonomy_url => $self->o('taxonomy_url'),
    dbsrv_url    => $self->o('dbsrv_url'),

    cs_order     => $self->o('cs_order'),
    prune_agp    => $self->o('prune_agp'),
    unversion_scaffolds => $self->o('unversion_scaffolds'),
    sr_syn_src_name  => $self->o('sr_syn_src_name'),
    division    => $self->o('division'),
    cs_tag_for_ordered => $self->o('cs_tag_for_ordered'),

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
    gff3_gene_types => $self->o('gff3_gene_types'),
    gff3_mrna_types => $self->o('gff3_mrna_types'),
    gff3_ignore_types => $self->o('gff3_ignore_types'),
    gff3_exon_types => $self->o('gff3_exon_types'),
    gff3_cds_types => $self->o('gff3_cds_types'),
    gff3_utr_types => $self->o('gff3_utr_types'),
    gff3_types_complete => $self->o('gff3_types_complete'),

    default_feature_version     => $self->o('default_feature_version'),
    no_feature_version_defaults => $self->o('no_feature_version_defaults'),

    load_pseudogene_with_CDS => $self->o('load_pseudogene_with_CDS'),
    no_brc4_stuff => $self->o('no_brc4_stuff'),
    ignore_final_stops => $self->o('ignore_final_stops'),
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
          proddb   => '#proddb_url#',
          taxonomy => '#taxonomy_url#',
        },
      },
      -meadow_type       => 'LSF',
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
        inputcmd => "find #data_dir# -type f -name manifest.json",
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
      -analysis_capacity => 1,
      -batch_size     => 50,
      -flow_into  => {
        2 => 'check_json_schema',
      },
    },

    { -logic_name     => 'check_json_schema',
      -module         => 'ensembl.brc4.runnable.schema_validator',
      -language => 'python3',
      -parameters     => {
        json_file => '#metadata_json#',
        json_schema => '#schemas#',
      },
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
    },

    {
      # Check the integrity of the manifest before loading anything
      -logic_name => 'Manifest_integrity',
      -module     => 'ensembl.brc4.runnable.integrity',
      -language   => 'python3',
      -parameters => {
        ignore_final_stops => $self->o('ignore_final_stops'),
      },
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -failed_job_tolerance => 100,
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
      -meadow_type       => 'LSF',
      -flow_into  => {
        '2' => 'CleanUpAndCreateDB',
      },
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
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into => { 1 => { 'LoadDBSchema' => INPUT_PLUS() } },
    },

    {
      -logic_name => 'LoadDBSchema',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        input_file => $self->o('ensembl_root_dir') . '/ensembl/sql/table.sql',
      },
      -analysis_capacity   => 1,
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into  => 'PopulateProductionTables',
    },

    {
      -logic_name    => "PopulateProductionTables",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #dump_path#; ' .
          ' perl #base_dir#/ensembl-production/scripts/production_database/populate_production_db_tables.pl '
            . ' --host #dbsrv_host# --port #dbsrv_port# --user #dbsrv_user# --pass #dbsrv_pass# --database #db_name# '
            . ' --mhost #proddb_host# --mport #proddb_port# --muser #proddb_user# --mdatabase #proddb_dbname# '
            . ' --dumppath #dump_path# --dropbaks '
            . ' > #dump_path#/stdout 2> #dump_path#/stderr ',
        'base_dir'       => $self->o('ensembl_root_dir'),
        'dump_path' => $self->o('pipeline_dir') . '/#db_name#/create_core/fill_production_db_tables',
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => 'LoadSequenceData',
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
        no_brc4_stuff => $self->o('no_brc4_stuff'),
      },
      -analysis_capacity   => 4,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
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
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
          -dbname => "#taxonomy_dbname#",
        },
        has_gff3 => '#expr( #manifest_data#->{"gff3"} )expr#',
      },
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into  => WHEN('#has_gff3#' => 'Load_gene_models'),
    },

    {
      -logic_name => 'Load_gene_models',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
        'has_func_annotation' => '#expr( #manifest_data#->{"functional_annotation"} )expr#',
      },
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
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
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
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
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
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
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into => 'LoadGFF3AnalysisSetup',
    },

    {
      -logic_name => 'LoadGFF3AnalysisSetup',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AnalysisSetup',
      -parameters        => {
        logic_name         => $self->o('gff3_load_logic_name'),
        module             => $self->o('gff3_load_analysis_module'),
        db_url             => '#dbsrv_url#' . '#db_name#',
        production_lookup  => $self->o('gff3_load_production_lookup'),
        prodb_url          => '#proddb_url#',
        delete_existing    => 1,
      },
      -rc_name    => 'default',
      -meadow_type => 'LSF',
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
        # feature types
        gene_types      => $self->o('gff3_gene_types'),
        mrna_types      => $self->o('gff3_mrna_types'),
        exon_types      => $self->o('gff3_exon_types'),
        cds_types       => $self->o('gff3_cds_types'),
        utr_types       => $self->o('gff3_utr_types'),
        ignore_types    => $self->o('gff3_ignore_types'),
        types_complete  => $self->o('gff3_types_complete'),
        polypeptides    => $self->o('gff3_use_polypeptides'), # it's better to ignore ignore 'polypeptides' lines
        load_pseudogene_with_CDS => $self->o('load_pseudogene_with_CDS'),
        # dbparams
        db_url          => '#dbsrv_url#' . '#db_name#',
        # condition
        has_fasta_peptide => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
      },
      -max_retry_count   => 0,
      -rc_name    => '15GB',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into => WHEN('#has_fasta_peptide#' => [ 'FixModels' ]),
    },

    {
      -logic_name => 'FixModels',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::FixModels',
      -parameters => {
        logic_name      => $self->o('gff3_load_logic_name'),
        db_url          => '#dbsrv_url#' . '#db_name#',
        protein_fasta_file      => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
      },
      -max_retry_count   => 0,
      -rc_name    => '15GB',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into => [ 'ApplySeqEdits' ],
    },

    {
      -logic_name => 'ApplySeqEdits',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ApplySeqEdits',
      -parameters => {
        logic_name      => $self->o('gff3_load_logic_name'),
        db_url          => '#dbsrv_url#' . '#db_name#',
        protein_fasta_file      => '#expr( #manifest_data#->{"fasta_pep"} )expr#',
      },
      -max_retry_count   => 0,
      -rc_name    => '15GB',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
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
      -rc_name    => '15GB',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
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
      -meadow_type       => 'LSF',
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
            . '  -external_db_map ' . $self->o('external_db_map')
            . '  > #log_path#/stdout '
            . '  2> #log_path#/stderr ',
        'log_path'       => $self->o('pipeline_dir') . '/#db_name#/load_functional_annotation',
        'fann_loader'    => "$scripts_dir/load_fann.pl",
        'fann_json_file' => '#expr( #manifest_data#->{"functional_annotation"} )expr#',
        'default_feat_v' => '#expr( #no_feature_version_defaults# ? "": "-feature_version_default ".#default_feature_version# )expr#',
      },
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into => 'Finalize_gene_models',
    },

    # Finish up annotations
    {
      -logic_name => 'Finalize_gene_models',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
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
      -meadow_type       => 'LSF',
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
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 2,
      -flow_into => 'PopulateAnalysis',
    },

    # Finalize database
    {
      -logic_name    => "PopulateAnalysis",
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #dump_path#;'
            . ' perl #base_dir#/ensembl-production/scripts/production_database/populate_analysis_description.pl '
            . ' --host #dbsrv_host# --port #dbsrv_port# --user #dbsrv_user# --pass #dbsrv_pass# --database #db_name# '
            . ' --mhost #proddb_host# --mport #proddb_port# --muser #proddb_user# --mdatabase #proddb_dbname# '
            . ' --dumppath #dump_path#'
            . ' --dropbak',
        'base_dir'       => $self->o('ensembl_root_dir'),
        'dump_path' => $self->o('pipeline_dir') . '/#db_name#/fill_production_analysis',
      },
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into => 'Change_seq_region_name',
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
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
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
