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
my $scripts_dir = "$package_dir/../../../../../scripts";

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    ##############################
    # Config to be set by the user
    # VB release id
    release => $self->o('release'),
    db_prefix => "",

    # Ensembl version (deduced from the environment?)
    ensembl_version => software_version(),

    # Meta configuration directory
    data_dir => $self->o('data_dir'),

    # Working directory
    tmp_dir => 'tmp',

    check_manifest => 1,
    cs_order => 'chunk,contig,supercontig,non_ref_scaffold,scaffold,superscaffold,linkage_group,chromosome',
    prune_agp => 0,
    unversion_scaffolds => 0,
    sr_syn_src_name  => 'ensembl_internal_synonym', # 50803
    division => 'EnsemblMetazoa',

    # GenomeTools aliases (from eg-pipelines: Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf)
    gt_exe        => 'gt',
    gff3_tidy     => $self->o('gt_exe').' gff3 -tidy -sort -retainids',
    gff3_validate => $self->o('gt_exe').' gff3validator',

    # GFF3 cleaning params
    gff3_ignore_file => catfile(dirname(__FILE__), qw/gff3.ignore/),
    gff3_clean_additionally => 0,
    gff3_autoapply_manual_seq_edits => 1,

    # LoadGFF3 params
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
        snRNA sRNA SRP_RNA tRNA
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

    ##############################

    # Basic pipeline configuration
    pipeline_tag => '',
    pipeline_name => 'brc4_genome_loader_' .
      $self->o('release') . '_' . $self->o('ensembl_version') . $self->o('pipeline_tag'),
    email => $ENV{USER} . '@ebi.ac.uk',
    pipeline_dir => 'genome_loader_' .  $self->o('release') . '_' . $self->o('ensembl_version'),

    debug => 0,
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    debug          => $self->o('debug'),
    tmp_dir        => $self->o('tmp_dir'),
    check_manifest => $self->o('check_manifest'),
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

    gt_exe        => $self->o('gt_exe'),
    gff3_tidy     => $self->o('gff3_tidy'),
    gff3_validate => $self->o('gff3_validate'),

    gff3_ignore_file            => $self->o('gff3_ignore_file'),
    gff3_clean_additionally     => $self->o('gff3_clean_additionally'),
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
  };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('tmp_dir'),
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
    # Create the tmp pipeline directory
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'FillDBParams',
        'A->1' => 'Cleanup',
      },
    },

    {
      # Delete the temp working directory
      -logic_name => 'Cleanup',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity   => 1,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
      -flow_into => 'Manifest_factory',
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
          'Manifest_integrity' => { manifest => '#_0#' }
        }, ELSE({
          'Prepare_genome' => {manifest => '#_0#' }
        })),
      },
    },

    {
      # Check the integrity of the manifest before loading anything
      -logic_name => 'Manifest_integrity',
      -module     => 'Integrity',
      -language => 'python3',
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -flow_into => 'Prepare_genome',
    },

    {
      # Prepare all the metadata:
      # - species
      # - db_name
      # - manifest_metadata
      -logic_name => 'Prepare_genome',
      -module     => 'PrepareGenome',
      -language => 'python3',
      -parameters => {
        release => $self->o('release'),
        ensembl_version => $self->o('ensembl_version'),
        db_prefix => $self->o('db_prefix'),
      },
      -analysis_capacity   => 1,
      -max_retry_count => 0,
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '2->A' => 'CleanUpAndCreateDB',
        'A->2' => 'LoadData',
      },
    },

    {
      # Init the Ensembl core
      -logic_name => 'CleanUpAndCreateDB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url'),
        sql     => [
          'SET GLOBAL max_allowed_packet=2147483648;',
          'DROP DATABASE IF EXISTS #db_name#;' ,
          'CREATE DATABASE #db_name#;' ,
        ],
      },
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into => [ 'LoadDBSchema' ],
    },

    {
      -logic_name => 'LoadDBSchema',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('dbsrv_url') . '#db_name#',
        input_file => $self->o('ensembl_root_dir') . '/ensembl/sql/table.sql',
      },
      -meadow_type       => 'LSF',
      -rc_name    => 'default',
      -flow_into  => [ 'PopulateProductionTables' ],
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
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },

    {
      # Head analysis for the loading of data
      -logic_name => 'LoadData',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'LoadSequenceData',
        'A->1' => 'LoadMetadata',
      },
    },

    {
      -logic_name => 'LoadSequenceData',
      -module     => 'LoadSequenceData',
      -language => 'python3',
      -parameters        => {
        work_dir => $self->o('pipeline_dir') . '/#db_name#/load_sequence',
        cs_order => $self->o('cs_order'),
        prune_agp => $self->o('prune_agp'),
        unversion_scaffolds => $self->o('unversion_scaffolds'),
        sr_syn_src  => $self->o('sr_syn_src_name'),
      },
      -analysis_capacity   => 2,
      -rc_name         => '8GB',
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
    },

    {
      # Head analysis for the loading of the metadata
      -logic_name => 'LoadMetadata',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'FillMetadata',
        'A->1' => 'ProcessRepeats',
      },
    },

    {
      -logic_name => 'FillMetadata',
      -module     => 'FillMetadata',
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
      -flow_into  => {
        '1' => 'FillTaxonomy',
      },
    },

    {
      -logic_name => 'FillTaxonomy',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        'cmd' => 'mkdir -p #log_path#; '
            . ' perl #base_dir#/ensembl-pipeline/scripts/load_taxonomy.pl '
            . '   --dbhost #dbsrv_host# --dbport #dbsrv_port# '
            . '   --dbuser #dbsrv_user# --dbpass #dbsrv_pass# '
            . '   --dbname #db_name# '
            . '   --taxondbhost #taxonomy_host# --taxondbport #taxonomy_port# ' # no taxondbuser
            . '   --taxondbname #taxonomy_dbname# '
            . '   --taxon_id #taxonomy_id# '
            . '   > #log_path#/stdout 2> #log_path#/stderr ',
        'taxonomy_id' => '#expr( #genome_data#->{"species"}->{"taxonomy_id"} )expr#',
        'base_dir' => $self->o('ensembl_root_dir'),
        'log_path' => $self->o('pipeline_dir') . '/#db_name#/create_core/fill_taxonomy',
      },
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
    },

    {
      -logic_name => 'ProcessRepeats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'ConstructRepeatLib',
        'A->1' => 'LoadGFF3Models',
      },
    },

    {
      -logic_name => 'ConstructRepeatLib',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into => [ 'AnnotateRepeats' ],
    },

    {
      -logic_name => 'AnnotateRepeats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'LoadGFF3Models',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'GFF3CleanAndLoad',
        'A->1' => 'LoadFunctionalAnnotation',
      },
    },

    {
      -logic_name => 'GFF3CleanAndLoad',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1->A' => 'GFF3CleanIgnored',
        'A->1' => 'LoadGFF3AnalysisSetup',
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
      -flow_into => WHEN('#gff3_clean_additionally#' =>
                      [ 'GFF3CleanAdditionally' ],
                    ELSE
                      [ 'GFF3Tidy' ]
                    ),
    },

    {
      -logic_name => 'GFF3CleanAdditionally',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
      -flow_into => [ 'GFF3Tidy' ],
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
        # dbparams
        db_url          => '#dbsrv_url#' . '#db_name#',
      },
      -max_retry_count   => 0,
      -rc_name    => '15GB',
      -meadow_type       => 'LSF',
      -analysis_capacity   => 5,
      -flow_into => [ 'FixModels' ],
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
      -logic_name => 'LoadFunctionalAnnotation',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LSF',
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
