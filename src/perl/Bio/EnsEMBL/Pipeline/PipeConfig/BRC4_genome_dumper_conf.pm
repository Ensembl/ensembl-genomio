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

 Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_dumper_conf;

=head1 DESCRIPTION

=head1 AUTHORS

 ckong@ebi.ac.uk
 mbarba@ebi.ac.uk

=cut
package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_dumper_conf;

use strict;
use warnings;
use File::Spec::Functions qw/catfile catdir/;
use Data::Dumper;
use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Class::Inspector;
use File::Basename;

my $package_path = Class::Inspector->loaded_filename(__PACKAGE__);
my $package_dir = dirname($package_path);
my $root_dir = "$package_dir/../../../../../..";

my $config_dir = "$root_dir/config";
my $runnables_dir = "$root_dir/src/python/ensembl/brc4/runnable";

use base ('Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_base_conf');

sub default_options {
	my ($self) = @_;

	return {
       ## inherit other stuff from the base class
       %{ $self->SUPER::default_options() },

       ## default LSF queue nane 
       queue_name => "standard",

	   ## General parameters
       'registry'      => $self->o('registry'),   
       'pipeline_name' => "brc4_genome_dumper",
	     'email'         => $self->o('ENV', 'USER').'@ebi.ac.uk',
       'output_dir'    => './output',
       'tmp_dir'    => './tmp',
       
       # No impact to the final files (so arbritrary value is fine),
       # but required by the gff3 and fasta dumper from Ensembl production
       'release' => 50,
      
       # Disable all to be able to select each part separately
       'dump_all' => 1,

       'do_fasta_dna' => $self->o('dump_all'),
       'do_fasta_pep' => $self->o('dump_all'),
       'do_gff' => $self->o('dump_all'),
       'do_agp' => $self->o('dump_all'),
       'do_events' => $self->o('dump_all'),
       # Json meta
        'do_func' => $self->o('dump_all'),
        'do_genome' => $self->o('dump_all'),
        'do_seq_attr' => $self->o('dump_all'),
        'do_seq_reg' => $self->o('dump_all'),

      # Other dumps
       'sql_dir' => "",

	   ## 'job_factory' parameters
	   'species'     => [], 
	   'antispecies' => [],
       'division' 	 => [], 
	   'run_all'     => 0,	

	   ## gff3 & gtf parameter
       'abinitio'        => 0,
       'gene' => 1,

       ## gff3 parameters
       'gt_exe'          => 'gt',
       'gff3_tidy'       => $self->o('gt_exe').' gff3 -tidy -sort -retainids -force',
       'gff3_validate'   => $self->o('gt_exe').' gff3validator',

       'feature_type'    => ['Gene', 'Transcript', 'SimpleFeature'], #'RepeatFeature'
       'per_chromosome'  => 0,
       'include_scaffold'=> 1,
	   'logic_name'      => [],
	   'db_type'	     => 'core',
       'out_file_stem'   => undef,
       'xrefs'           => 0,

       # GFF3 specifications
       # Ensembl compatibility
       ensembl_mode => 0,

	   ## fasta parameters
      dump_level => 'seqlevel', # Alternative: toplevel

       # types to emit
       'dna_sequence_type_list'  => ['dna'],
       'pep_sequence_type_list'  => ['cdna', 'ncrna'],

       # Do/Don't process these logic names
       'process_logic_names' => [],
       'skip_logic_names'    => [],

      # Map back the external db names
      external_db_map => catfile($config_dir, 'external_db_map', 'default.txt'),
	};
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
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

# these parameter values are visible to all analyses, 
# can be overridden by parameters{} and input_id{}
sub pipeline_wide_parameters {  
    my ($self) = @_;
    return {
            %{$self->SUPER::pipeline_wide_parameters},  # here we inherit anything from the base class
		    'pipeline_name' => $self->o('pipeline_name'), #This must be defined for the beekeeper to work properly
            'base_path'     => $self->o('tmp_dir'),
            'output_dir'     => $self->o('output_dir'),
            'do_fasta_dna'      => $self->o('do_fasta_dna'),
            'do_fasta_pep'      => $self->o('do_fasta_pep'),
            'do_gff'      => $self->o('do_gff'),
            'do_agp'      => $self->o('do_agp'),
            'do_gff'      => $self->o('do_gff'),
            'do_events'      => $self->o('do_events'),
            'do_genome'      => $self->o('do_genome'),
            'do_func'      => $self->o('do_func'),
            'do_seq_reg'      => $self->o('do_seq_reg'),
            'do_seq_attr'      => $self->o('do_seq_attr'),
            'sql_dir'      => $self->o('sql_dir'),
            'schemas'      => $self->o('schemas'),
            #'remove_features_prefix'      => $self->o('remove_features_prefix'),
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    return [
    {
      -logic_name        => 'Start',
      -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -input_ids      => [ {} ],
      -parameters        => {
        cmd => 'mkdir -p '. $self->o('tmp_dir') . ' ; mkdir -p '. $self->o('output_dir'),
      },
      -rc_name           => 'default',
      -max_retry_count   => 0,
      -flow_into       => 'Species_factory',
    },

	 {
       -logic_name     => 'Species_factory',
       -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
       -parameters     => {
                             species     => $self->o('species'),
                             antispecies => $self->o('antispecies'),
                             division    => $self->o('division'),
                             run_all     => $self->o('run_all'),
                          },
	    -hive_capacity   => -1,
      -rc_name 	       => 'default',
      -max_retry_count => 0,
      -flow_into       => {
                           '2->A' => 'Files_makers',
                           'A->2' => 'Manifest_maker',
                           '2' => WHEN('#sql_dir#', 'Dump_genome_sql'),
                         },
    },
    { -logic_name  => 'Manifest_maker',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::Manifest',
      -max_retry_count => 0,
      -analysis_capacity   => 1,
      -rc_name         => 'default',
      -flow_into       => { '2' => 'Manifest_check' },
    },

     { -logic_name     => 'Manifest_check',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters     => {
         json_file => '#manifest#',
         metadata_type => 'manifest',
         json_schema => '#expr(${#schemas#}{#metadata_type#})expr#',
         cmd => 'jsonschema -i #json_file# #json_schema#',
       },
       -max_retry_count => 0,
       -analysis_capacity => 1,
       -batch_size     => 50,
	     -rc_name        => 'default',
      -flow_into       => { '1' => 'Integrity_check' },
     },

    { -logic_name  => 'Integrity_check',
      -module      => 'ensembl.brc4.runnable.integrity',
      -language    => 'python3',
      -parameters     => {
        ensembl_mode => $self->o('ensembl_mode'),
      },
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
    },
 	
     { -logic_name     => 'Files_makers',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -hive_capacity  => -1,
       -rc_name 	   => 'default',
       -analysis_capacity => 1,
       -flow_into      => {'1' => [
           'Fasta_makers',
           'Metadata_makers',
           WHEN('#do_gff#', 'GFF3_maker'),
           WHEN('#do_agp#', 'AGP_maker'),
           WHEN('#do_events#', 'Extract_connection'),
         ] }
     },

### Ids events
    {
      -logic_name => 'Extract_connection',
      -module     => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::ExtractConnection',
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
      -flow_into  => { 2 => 'Events_dumper' },
    },

    {
      -logic_name      => 'Events_dumper',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters      => {
        events_dir => catdir('#base_path#', 'events', '#species#'),
        output_file => catfile('#events_dir#', 'ids_events.tab'),
        hash_key => 'events',
        cmd => "mkdir -p #events_dir#;"
          . " python $runnables_dir/dump_stable_ids.py"
          . " --host #host# "
          . " --port #port# "
          . " --user #user# "
          . " --password '#password#' "
          . " --dbname #dbname# "
          . " --output_file #output_file#",
      },
      -max_retry_count => 0,
      -analysis_capacity => 2,
      -rc_name         => 'default',
      -flow_into       => {
        1 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=output_file'
      },
    },


### GFF3
     { -logic_name     => 'GFF3_maker',
       -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
       -parameters     => {
          feature_type       => $self->o('feature_type'),
          per_chromosome     => $self->o('per_chromosome'),
          include_scaffold   => $self->o('include_scaffold'),
          logic_name         => $self->o('logic_name'),
          db_type            => $self->o('db_type'),
	      abinitio           => $self->o('abinitio'),
        gene               => $self->o('gene'),
	      out_file_stem      => $self->o('out_file_stem'),
	      xrefs              => $self->o('xrefs'),        
        release  => $self->o('release'),
        },
      -max_retry_count => 0,
       -hive_capacity  => 20,
       -priority        => 5,
       -rc_name 	   => 'default',
       -flow_into      => {
         '-1' => 'GFF3_maker_highmem',
         '1'  => 'GFF3_BRC4_filtering'
       },
     },

	 { -logic_name     => 'GFF3_maker_highmem',
       -module         => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
       -parameters     => {
          feature_type       => $self->o('feature_type'),
          per_chromosome     => $self->o('per_chromosome'),
          include_scaffold   => $self->o('include_scaffold'),
          logic_name         => $self->o('logic_name'),
          db_type            => $self->o('db_type'),
    	  abinitio           => $self->o('abinitio'),
        gene               => $self->o('gene'),
	      out_file_stem      => $self->o('out_file_stem'),
	      xrefs              => $self->o('xrefs'),        
        release  => $self->o('release'),
        },
      -max_retry_count => 0,
	    -hive_capacity  => 20,
      -priority        => 5,
      -rc_name        => '32GB',
      -flow_into      => { '1'  => 'GFF3_BRC4_filtering' },
	 },	

### GFF3:post-processing
   # This only allows the one type of file necessary for BRC4
   { -logic_name  => 'GFF3_BRC4_filtering',
     -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::FilterGFF3',
      -max_retry_count => 0,
     -batch_size     => 10,
     -rc_name        => 'default',
     -flow_into      => { 2 =>'GFF3_BRC4_specifications' },
   },

   # BRC4 specifications alterations
   { -logic_name  => 'GFF3_BRC4_specifications',
     -module      => 'ensembl.brc4.runnable.gff3_specifier',
     -language    => 'python3',
     -parameters     => {
       gff_file => '#filtered_gff_file#',
       ensembl_mode => $self->o("ensembl_mode"),
     },
      -max_retry_count => 0,
     -batch_size     => 10,
     -rc_name        => 'default',
     -flow_into      => { 2 =>'GFF3_validation' },
   },
   
   # Sort the features
   # The validation step also adds the complete list of seq_regions to the header
   { -logic_name     => 'GFF3_validation',
     -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters     => {
       cmd => $self->o('gff3_tidy').' -gzip -o #final_gff_file# #specifications_gff_file# && ' .
       $self->o('gff3_validate').' #final_gff_file#',
       hash_key => "gff3",
     },
      -max_retry_count => 0,
     -analysis_capacity => 10,
     -batch_size        => 10,
     -rc_name           => 'default',
     -flow_into  => { 1 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=final_gff_file' },
   },


### FASTA (cdna, cds, dna, pep, ncrna)
    { -logic_name  => 'Fasta_makers',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 1,
      -flow_into  => [
        WHEN('#do_fasta_dna#', 'Fasta_DNA'),
        WHEN('#do_fasta_pep#', 'Fasta_peptide'),
      ],
      -rc_name         => 'default',
    },
    # 
    { -logic_name  => 'Fasta_DNA',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaDNA',
      -parameters => {
        hash_key => "fasta_dna",
        dump_level => $self->o('dump_level'),
      },
      -flow_into  => { 2 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=fasta_file' },
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
    },

    { -logic_name  => 'Fasta_peptide',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaPeptide',
      -parameters => { hash_key => "fasta_pep", },
      -flow_into  => { 2 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=fasta_file' },
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
    },

### AGP
    { -logic_name  => 'AGP_maker',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpAGP',
      -parameters => { hash_key => "agp" },
      -flow_into  => { 2 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=agp_files' },
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
    },

### METADATA
    { -logic_name  => 'Metadata_makers',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 1,
      -flow_into  => [
        WHEN('#do_seq_reg#', 'Seq_region'),
        WHEN('#do_seq_attr#', 'Seq_attrib'),
        WHEN('#do_func#', 'Functional_annotation'),
        WHEN('#do_genome#', 'Genome'),
      ],
      -rc_name         => 'default',
    },

    { -logic_name  => 'Seq_region',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqRegionJson',
      -parameters     => {
        external_db_map => $self->o('external_db_map'),
        dump_level => $self->o('dump_level'),
      },
      -flow_into  => { 2 => ['Check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'Seq_attrib',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqAttribJson',
      -flow_into  => { 2 => ['Check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'Functional_annotation',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFunctionalAnnotationJson',
      -parameters     => {
        external_db_map => $self->o('external_db_map'),
      },
      -flow_into  => { 2 => ['Check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'Genome',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeJson',
      -flow_into  => { 2 => ['Check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    {
      -logic_name     => 'Check_json_schema',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters  => {
        log_path => $self->o('tmp_dir') . '/check_schemas',
        json => '#metadata_json#',
        cmd => 'mkdir -p #log_path#; '
             . 'echo "checking #json# against #metadata_type#" > #log_path#/#metadata_type#.log; '
             . 'schemas_json_validate --json_file #json# --json_schema #metadata_type# '
             . '   >> #log_path#/#metadata_type#.log 2>&1 ',
        hash_key => "#metadata_type#",
      },
      -analysis_capacity => 2,
      -failed_job_tolerance => 10,
      -batch_size     => 50,
      -rc_name        => 'default',
      -flow_into  => { 1 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=metadata_json' },
    },

    { -logic_name  => 'Dump_genome_sql',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeSQL',
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },
  ];
}

1;

