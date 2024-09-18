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

 mbarba@ebi.ac.uk

=cut
package Bio::EnsEMBL::Pipeline::PipeConfig::BRC4_genome_compare_conf;

use strict;
use warnings;
use File::Spec::Functions qw/catfile catdir/;
use Data::Dumper;
use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Class::Inspector;
use File::Basename;

use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        ## inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

        ## General parameters
        'registry'      => $self->o('registry'),   
        'pipeline_name' => "brc4_genome_compare",
        'email'         => $self->o('ENV', 'USER').'@ebi.ac.uk',
        'output_dir'    => './output',
        'tmp_dir'    => './tmp',
      
        ## 'job_factory' parameters
        'species'     => [], 
        'division'    => [], 
        'antispecies' => [],
        'run_all'     => 0,	

        ## default LSF queue
	queue_name => 'standard',
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
            download_dir   => catdir($self->o('tmp_dir'), "download", '#species#'),
    };
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
                           'A->1' => 'ReportComparisons',
                         },
    },

	 { -logic_name     => 'Species_factory',
       -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
       -parameters     => {
                             species     => $self->o('species'),
                             division    => $self->o('division'),
                             antispecies => $self->o('antispecies'),
                             run_all     => $self->o('run_all'),
                          },
	    -analysis_capacity   => 1,
      -rc_name 	       => 'default',
      -max_retry_count => 0,
      -flow_into       => {
                           '2->A' => 'Files_makers',
                           'A->2' => 'Compare',
                         },
    },
 	
     { -logic_name     => 'Files_makers',
       -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -analysis_capacity  => -1,
       -rc_name 	   => 'default',
       -analysis_capacity => 1,
       -batch_size   => 100,
       -flow_into      => {
         '1' => [
           WHEN('not -e #fasta_dna_file#', 'Fasta_DNA'),
           WHEN('not -e #seq_region_file#', 'Seq_region'),
           'Get_accession',
         ]
       }
     },

    { -logic_name  => 'Fasta_DNA',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaDNA',
      -parameters => {
        dump_level => "toplevel",
      },
      -max_retry_count => 0,
      -analysis_capacity   => 20,
      -priority        => 5,
      -rc_name         => 'default',
      -flow_into  => {
        2 => '?accu_name=core_fasta_dna&accu_input_variable=fasta_file',
      },
    },

    {
      -logic_name  => 'Seq_region',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqRegionJson',
      -parameters     => {
        dump_level => 'toplevel',
      },
      -max_retry_count => 0,
      -analysis_capacity  => 20,
      -rc_name         => 'default',
      -flow_into  => {
        2 => '?accu_name=seq_region_json&accu_input_variable=metadata_json',
      },
    },

    # INSDC download
    { -logic_name  => 'Get_accession',
      -module      => 'Bio::EnsEMBL::Pipeline::Runnable::BRC4::GetMetaValue',
      -parameters => {
        provider => "assembly.provider_name",
        param_name => "accession",
        param_key => "assembly.accession",
      },
      -max_retry_count => 0,
      -analysis_capacity   => 1,
      -batch_size   => 100,
      -rc_name         => 'default',
      -flow_into  => {
        2 => [
           WHEN('not -e #insdc_dna_file#', 'INSDC_download'),
         ]
      },
    },
    {
      -logic_name     => 'INSDC_download',
      -module         => 'ensembl.brc4.runnable.download_assembly_data',
      -parameters => {
        accession => "#accession#",
        max_increment => 2,
      },
      -language => 'python3',
      -analysis_capacity => 1,
      -failed_job_tolerance => 100,
      -rc_name        => 'default',
      -flow_into  => {
        2 => [
          '?accu_name=insdc_fasta_dna&accu_input_variable=fasta_dna',
          '?accu_name=insdc_report&accu_input_variable=report',
          '?accu_name=accession&accu_input_variable=accession',
        ],
      },
    },

    # Compare
    {
      -logic_name => 'Compare',
      -module         => 'ensembl.brc4.runnable.compare_fasta',
      -parameters => {
        report => "#insdc_report#",
        fasta1 => "#insdc_fasta_dna#",
        fasta2 => "#core_fasta_dna#",
        seq_regions => "#seq_region_json#",
        comparison_name => "fasta_dna",
        accession => "#accession#",
      },
      -language => 'python3',
      -analysis_capacity => 5,
      -failed_job_tolerance => 0,
      -rc_name        => '8GB',
      -flow_into      => {
         '2' => '?accu_name=stats&accu_address={species}'
        }
    },

    # Report comparisons
    {
      -logic_name => 'ReportComparisons',
      -module         => 'ensembl.brc4.runnable.compare_report',
      -language => 'python3',
      -analysis_capacity => 1,
      -failed_job_tolerance => 0,
      -rc_name        => 'default',
    },
  ];
}

1;
