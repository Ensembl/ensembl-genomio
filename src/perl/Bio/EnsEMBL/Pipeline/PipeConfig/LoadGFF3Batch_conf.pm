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


=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::EG::PipeConfig::LoadGFF3Batch_conf

=head1 DESCRIPTION

Load valid GFF3 files into core databases.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::PipeConfig::LoadGFF3Batch_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Pipeline::PipeConfig::LoadGFF3_conf');


sub default_options {
  my ($self) = @_;
  return {
    %{$self->default_options_generic},
    
    pipeline_name => 'load_gff3_batch',
    
    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},
    
    # The code assumes that for all species there will be an appropriately
    # named subdirectory in the 'source_dir', and within that subdir the
    # GFF3 files will have the same name. Ditto other file types.
    gff3_tidy_filename     => $self->o('gff3_filename').'.tidied',
    fasta_filename         => undef,
    genbank_filename       => undef,
    protein_fasta_filename => undef,
    
    # Since this mode was designed to be used with patch builds,
    # default to not deleting any existing genes.
    delete_existing => 0,
    keep_logic_name => 1,
  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  
  return $self->pipeline_wide_parameters_generic;
}

sub pipeline_analyses {
  my ($self) = @_;
  
  return [
    {
      -logic_name      => 'SpeciesFactory',
      -module          => 'Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory',
      -input_ids       => [ {} ],
      -parameters      => {
                            species         => $self->o('species'),
                            antispecies     => $self->o('antispecies'),
                            division        => $self->o('division'),
                            run_all         => $self->o('run_all'),
                            meta_filters    => $self->o('meta_filters'),
                            chromosome_flow => 0,
                            regulation_flow => 0,
                            variation_flow  => 0,
                          },
      -max_retry_count => 1,
      -flow_into       => {
                            2 => ['ConstructFilenames'],
                          },
      -meadow_type     => 'LOCAL',
    },

    {
      -logic_name      => 'ConstructFilenames',
      -module          => 'Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ConstructFilenames',
      -max_retry_count => 0,
      -parameters      => {
                            source_dir             => $self->o('source_dir'),
                            gff3_filename          => $self->o('gff3_filename'),
                            gff3_tidy_filename     => $self->o('gff3_tidy_filename'),
                            fasta_filename         => $self->o('fasta_filename'),
                            genbank_filename       => $self->o('genbank_filename'),
                            protein_fasta_filename => $self->o('protein_fasta_filename'),
                          },
      -flow_into       => ['LoadGFF3Start'],
      -meadow_type     => 'LOCAL',
    },
    
    @{ $self->pipeline_analyses_generic}
  ];
}

1;
