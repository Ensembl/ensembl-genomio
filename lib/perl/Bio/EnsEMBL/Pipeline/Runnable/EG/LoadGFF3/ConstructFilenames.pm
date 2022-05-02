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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ConstructFilenames

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ConstructFilenames;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

use File::Spec::Functions qw(catdir);

sub write_output {
  my ($self) = @_;
  my $species                = $self->param_required('species');
  my $source_dir             = $self->param_required('source_dir');
  my $gff3_filename          = $self->param_required('gff3_filename');
  my $gff3_tidy_filename     = $self->param_required('gff3_tidy_filename');
  my $fasta_filename         = $self->param('fasta_filename');
  my $genbank_filename       = $self->param('genbank_filename');
  my $protein_fasta_filename = $self->param('protein_fasta_filename');
  
  my $output_ids = {
    gff3_file => catdir($source_dir, $species, $gff3_filename),
    gff3_tidy_file => catdir($source_dir, $species, $gff3_tidy_filename),
  };
  
  if (defined $fasta_filename) {
    $$output_ids{fasta_file} = catdir($source_dir, $species, $fasta_filename);
  } else {
    $$output_ids{fasta_file} = catdir($source_dir, $species, "genome.fa");
  }
  
  if (defined $genbank_filename) {
    $$output_ids{genbank_file} = catdir($source_dir, $species, $genbank_filename);
  }
  
  if (defined $protein_fasta_filename) {
    $$output_ids{protein_fasta_file} = catdir($source_dir, $species, $protein_fasta_filename);
  }
  
  $self->dataflow_output_id($output_ids, 1);
}

1;
