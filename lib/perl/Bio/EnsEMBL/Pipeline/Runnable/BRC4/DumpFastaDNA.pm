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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaDNA;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

use Bio::EnsEMBL::Utils::IO::FASTASerializer;

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);

sub param_defaults {
  my ($self) = @_;
  
  return {
    'chunk_factor'          => 1000,
    'line_width'            => 80,

    # DNA params
    # Alternative: toplevel
    'dump_level'            => 'seqlevel',
    'dump_cs_version'       => undef,
    'include_non_reference' => 1,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $fasta_type = 'dna';

  my $sub_dir = $self->create_dir('fasta_' . $fasta_type);
  my $fasta_file = catfile($sub_dir, $self->param_required("species") . '_' . $fasta_type . '.fa');
  $self->param('fasta_file', $fasta_file);
}

sub run {
  my ($self) = @_;

  my $fasta_file            = $self->param('fasta_file');
  my $chunk_factor          = $self->param('chunk_factor');
  my $line_width            = $self->param('line_width');

  # DNA params
  my $dump_level            = $self->param('dump_level');
  my $dump_cs_version       = $self->param('dump_cs_version');
  my $include_non_reference = $self->param('include_non_reference');

  my $header_function = sub {
    my $slice = shift;
    return $slice->seq_region_name;
  };
  open(my $fh, '>', $fasta_file) or $self->throw("Cannot open file $fasta_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new(
    $fh,
    $header_function,
    $chunk_factor,
    $line_width,
  );

  my $dba = $self->core_dba();
  my $sa = $dba->get_adaptor('Slice');
  my $slices = $sa->fetch_all($dump_level, $dump_cs_version, $include_non_reference);

  foreach my $slice (sort { $b->length <=> $a->length } @$slices) {
    $serializer->print_Seq($slice);
  }

  close($fh);
}

sub write_output {
  my ($self) = @_;
  $self->dataflow_output_id({ 'fasta_file' => $self->param('fasta_file')}, 2);
}

1;
