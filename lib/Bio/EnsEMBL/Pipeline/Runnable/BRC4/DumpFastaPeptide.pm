=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpFastaPeptide;

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

    # pep params
    'is_canonical' => undef,
    'allow_stop_codons' => 0,
  };
}

sub fetch_input {
  my ($self) = @_;
  
  $self->param( 'dba', $self->core_dba() );
  my $fasta_type = 'pep';

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
  my $is_canonical            = $self->param('is_canonical');
  my $allow_stop_codons       = $self->param('allow_stop_codons');

  open(my $fh, '>', $fasta_file) or $self->throw("Cannot open file $fasta_file: $!");
  my $serializer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new(
    $fh,
    undef,
    $chunk_factor,
    $line_width,
  );

  # Load transcripts with CDS
  my $dba = $self->core_dba();
  my $tra = $dba->get_adaptor('Transcript');
  my @coding_biotypes = qw(protein_coding pseudogene_with_CDS nontranslating_CDS);

  my @transcripts;
  for my $biotype (@coding_biotypes) {
    push @transcripts, @{ $tra->fetch_all_by_biotype($biotype) };
  }

  foreach my $transcript (sort { $a->stable_id cmp $b->stable_id } @transcripts) {
    if (defined $is_canonical) {
      next if $is_canonical != $transcript->is_canonical;
    }

    if ($transcript->translation) {
      my $seq_obj = $transcript->translate();
      $seq_obj->display_id($transcript->translation->stable_id);

      if ($seq_obj->seq() =~ /\*/ && !$allow_stop_codons) {
        $self->warning("Translation for transcript ".$transcript->stable_id." contains stop codons. Skipping.");
      } else {
        if ($seq_obj->seq() =~ /\*/) {
          $self->warning("Translation for transcript ".$transcript->stable_id." contains stop codons.");
        }
        $serializer->print_Seq($seq_obj);
      }
    }
  }

  close($fh);
}

sub write_output {
  my ($self) = @_;
  $self->dataflow_output_id({ 'fasta_file' => $self->param('fasta_file')}, 2);
}

1;
