=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
use File::Spec::Functions qw/catfile/;
use Data::Dumper;
use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Class::Inspector;
use File::Basename;

my $package_path = Class::Inspector->loaded_filename(__PACKAGE__);
my $package_dir = dirname($package_path);
my $schema_dir = "$package_dir/../../../../../schema";

use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
	my ($self) = @_;

	return {
       ## inherit other stuff from the base class
       %{ $self->SUPER::default_options() },

	   ## General parameters
       'registry'      => $self->o('registry'),   
       'release'       => $self->o('release'),
       'pipeline_name' => "brc4_genome_dumper",
	     'email'         => $self->o('ENV', 'USER').'@ebi.ac.uk',
       'output_dir'    => './output',
       'tmp_dir'    => './tmp',

       'do_fasta' => 1,
       'do_gff' => 1,
       'do_meta' => 1,
       'do_agp' => 1,

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


      ## Metadata parameters
      'schemas' => {
        'seq_region' => catfile($schema_dir, "seq_region_schema.json"),
        'seq_attrib' => catfile($schema_dir, "seq_attrib_schema.json"),
        'functional_annotation' => catfile($schema_dir, "functional_annotation_schema.json"),
        'genome' => catfile($schema_dir, "genome_schema.json"),
        'manifest' => catfile($schema_dir, "manifest_schema.json"),
      },
      # Map back the external db names
      external_db_map => undef,
	};
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
      # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      'mkdir -p '.$self->o('tmp_dir'),
      'mkdir -p '.$self->o('output_dir'),
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
            'release'       => $self->o('release'),
            'do_fasta'      => $self->o('do_fasta'),
            'do_gff'      => $self->o('do_gff'),
            'do_agp'      => $self->o('do_agp'),
            'do_gff'      => $self->o('do_gff'),
            'do_meta'      => $self->o('do_meta'),
            'schemas'      => $self->o('schemas'),
            #'remove_features_prefix'      => $self->o('remove_features_prefix'),
    };
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

sub pipeline_analyses {
    my ($self) = @_;
    
    return [
    {  -logic_name => 'Start',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -rc_name       => 'default',
      -input_ids      => [ {} ],
      -flow_into       => {
                           '1->A' => 'Species_factory',
                           'A->1' => 'Cleanup',
                         },
    },

    {  -logic_name => 'Cleanup',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -hive_capacity => 10,
       -rc_name       => 'default',
       -analysis_capacity => 1,
    },

	 { -logic_name     => 'Species_factory',
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
                           '2->A' => 'backbone_job_pipeline',
                           'A->2' => 'manifest',
                         },
    },
    { -logic_name  => 'manifest',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::Manifest',
      -max_retry_count => 0,
      -analysis_capacity   => 1,
      -rc_name         => 'default',
      -flow_into       => { '2' => 'check_manifest' },
    },

     { -logic_name     => 'check_manifest',
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
      -flow_into       => { '1' => 'integrity' },
     },

    { -logic_name  => 'integrity',
      -module      => 'Integrity',
      -language    => 'python3',
      -parameters     => {
        ensembl_mode => $self->o('ensembl_mode'),
      },
      -analysis_capacity   => 5,
      -rc_name         => '8GB',
      -max_retry_count => 0,
    },
 	
     { -logic_name     => 'backbone_job_pipeline',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -hive_capacity  => -1,
       -rc_name 	   => 'default',
       -analysis_capacity => 1,
       -flow_into      => {'1' => [
           WHEN('#do_fasta#', 'fasta'),
           WHEN('#do_gff#', 'gff3'),
           WHEN('#do_agp#', 'agp'),
           WHEN('#do_meta#', 'metadata'),
         ] }
     },

### GFF3
     { -logic_name     => 'gff3',
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
        },
      -max_retry_count => 0,
       -hive_capacity  => 20,
       -priority        => 5,
       -rc_name 	   => 'default',
       -flow_into      => { '-1' => 'gff3_highmem', '1'  => 'gff3_BRC4_filter' },
     },

	 { -logic_name     => 'gff3_highmem',
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
        },
      -max_retry_count => 0,
	    -hive_capacity  => 20,
      -priority        => 5,
      -rc_name        => '32GB',
      -flow_into      => { '1'  => 'gff3_BRC4_filter' },
	 },	

### GFF3:post-processing
   # This only allows the one type of file necessary for BRC4
   { -logic_name  => 'gff3_BRC4_filter',
     -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::FilterGFF3',
      -max_retry_count => 0,
     -batch_size     => 10,
     -rc_name        => 'default',
     -flow_into      => { 2 =>'gff3_BRC4_specifier' },
   },

   # BRC4 specifications alterations
   { -logic_name  => 'gff3_BRC4_specifier',
     -module      => 'GFF3Specifier',
     -language    => 'python3',
     -parameters     => {
       gff_file => '#filtered_gff_file#',
       ensembl_mode => $self->o("ensembl_mode"),
     },
      -max_retry_count => 0,
     -batch_size     => 10,
     -rc_name        => 'default',
     -flow_into      => { 2 =>'gff3_validation' },
   },
   
   # Sort the features
   # The validation step also adds the complete list of seq_regions to the header
   { -logic_name     => 'gff3_validation',
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
    { -logic_name  => 'fasta',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 1,
      -flow_into  => ['fasta_dna', 'fasta_pep'],
      -rc_name         => 'default',
    },
    # 
    { -logic_name  => 'fasta_dna',
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

    { -logic_name  => 'fasta_pep',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaPeptide',
      -parameters => { hash_key => "fasta_pep", },
      -flow_into  => { 2 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=fasta_file' },
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
    },

### AGP
    { -logic_name  => 'agp',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpAGP',
      -parameters => { hash_key => "agp" },
      -flow_into  => { 2 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=agp_files' },
      -max_retry_count => 0,
      -hive_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
    },

### METADATA
    { -logic_name  => 'metadata',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -analysis_capacity => 1,
      -flow_into  => [
        'seq_region',
        'seq_attrib',
        'functional_annotation',
        'genome',
      ],
      -rc_name         => 'default',
    },

    { -logic_name  => 'seq_region',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqRegionJson',
      -parameters     => {
        external_db_map => $self->o('external_db_map'),
      },
      -flow_into  => { 2 => ['check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'seq_attrib',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqAttribJson',
      -flow_into  => { 2 => ['check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'functional_annotation',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFunctionalAnnotationJson',
      -parameters     => {
        external_db_map => $self->o('external_db_map'),
      },
      -flow_into  => { 2 => ['check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name  => 'genome',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeJson',
      -flow_into  => { 2 => ['check_json_schema'] },
      -max_retry_count => 0,
      -hive_capacity  => 20,
      -rc_name         => 'default',
    },

    { -logic_name     => 'check_json_schema',
      -module         => 'SchemaValidator',
      -language => 'python3',
      -parameters     => {
        json_file => '#metadata_json#',
        json_schema => '#schemas#',
        hash_key => "#metadata_type#",
      },
      -flow_into  => { 1 => '?accu_name=manifest&accu_address={hash_key}&accu_input_variable=metadata_json' },
      -analysis_capacity => 1,
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -batch_size     => 50,
      -rc_name        => 'default',
    },
  ];
}

1;

