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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::FixModels

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::FixModels;

use strict;
use warnings;
use feature 'say';

use List::Util qw(min);
use Path::Tiny qw(path);

use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use base ('Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base');

sub run {
  my ($self) = @_;
  my $logic_name         = $self->param_required('logic_name');
  my $protein_fasta_file = $self->param_required('protein_fasta_file');

  my $dba = $self->get_type_dba();
  
  $self->fix_models($dba, $logic_name, $protein_fasta_file);
  
  $self->set_protein_coding($dba, $logic_name);

  $dba->dbc->disconnect_if_idle();

  $self->dump_log();
}

sub fix_models {
  my ($self, $dba, $logic_name, $protein_fasta_file) = @_;
  
  # We attempt two fairly basic fixes, by shifting the model back and
  # forth to counter two common problems:
  # 1. Out-of-frame problems, caused by the wrong phase. These can be
  #    fixed by shuffling the translation start 1 or 2 bases.
  # 2. Incorrect gene positions, caused by a mismatch between the
  #    sequence the gene was annotated on, and the INSDC sequence.
  # These two cases are considered independently; if a gene is affected
  # by both, it won't be fixed.
  
  my $sa = $dba->get_adaptor('Slice');
  my $ta = $dba->get_adaptor('Transcript');
  
  my %protein = $self->load_fasta($protein_fasta_file);
  
  my $slices = $sa->fetch_all( 'toplevel', undef, 0, 1 );

  # split up a list of slices into smaller slices
  my $overlap    = 10000;
  my $max_length = 1e6;
  $slices     = split_Slices( $slices, $max_length, $overlap);
  my %done_transcripts;

  for my $slice (@$slices) {
    my $transcripts = $ta->fetch_all_by_Slice($slice, 0, $logic_name);

    foreach my $transcript (sort { $a->seq_region_start cmp $b->seq_region_start } @$transcripts) {
      next if $done_transcripts{ $transcript->dbID }++;
      if ($transcript->biotype eq 'nontranslating_CDS') {
        my $translation = $transcript->translation;

        if ($translation && exists $protein{$translation->stable_id}) {
          my $db_seq   = $translation->seq;
          my $file_seq = $protein{$translation->stable_id};
          $file_seq =~ s/\*$//;

          if ($db_seq ne $file_seq) {
            my $success = $self->shift_translation($dba, $transcript, $file_seq);

            if ($success) {
              $self->log_warning(
                'WARNING: Shifted translation (nontranslating_CDS) of ' . $transcript->stable_id);
            } else {
              $success = $self->shift_gene($dba, $transcript, $file_seq);

              if ($success) {
                $self->log_warning(
                  'WARNING: Shifted position (nontranslating_CDS) of ' . $transcript->stable_id);
              }
            }
          }
        }
      }
    }
  }

}

sub shift_translation {
  my ($self, $dba, $transcript, $target_seq) = @_;
  
  my $dbid = $transcript->translation->dbID;
  my $original_start = $transcript->translation->start;
  
  my @shift_by = (-2..2);
  my $success  = 0;
  
  foreach my $shift_by (@shift_by) {
    my $new_start = $original_start + $shift_by;
    if ($new_start > 0) {
      $self->update_translation_start($dba, $dbid, $new_start);
      $transcript = $transcript->adaptor->fetch_by_stable_id($transcript->stable_id);
      
      if ($transcript->translation->seq eq $target_seq) {
        $success = 1;
        last;
      }
    }
  }
  
  if (! $success) {
    $self->update_translation_start($dba, $dbid, $original_start);
    $transcript = $transcript->adaptor->fetch_by_stable_id($transcript->stable_id);
  }
  
  return $success;
}

sub update_translation_start {
  my ($self, $dba, $dbid, $start) = @_;
  
  my $sql = 'UPDATE translation SET seq_start = ? WHERE translation_id = ?;';
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($start, $dbid) or $self->throw("Failed to execute: $sql");
}

sub shift_gene {
  my ($self, $dba, $transcript, $target_seq) = @_;
  
  my $dbid = $transcript->get_Gene->dbID;
  my $original_start = $transcript->translation->start;
  
  # A large degree of negative shifting is allowed because sometimes
  # genes are on scaffolds that have leading Ns in the source data.
  # These are stripped out on submission to INSDC, effectively shifting
  # the gene model backwards. I don't know of any mechanism by which
  # models will be shifted forwards, so we don't try that.
  
  my $shift     = 0;
  my $max_shift = min($transcript->seq_region_start, 25);
  my $success   = 0;
  
  while ($shift-- > -$max_shift) {
    $self->update_gene_position($dba, $dbid, -1);
    $transcript = $transcript->adaptor->fetch_by_stable_id($transcript->stable_id);
    
    if ($transcript->translation->seq eq $target_seq) {
      $success = 1;
      last;
    }
  }
  
  if (! $success) {
    $self->update_gene_position($dba, $dbid, $max_shift);
    $transcript = $transcript->adaptor->fetch_by_stable_id($transcript->stable_id);
  }
  
  return $success;
}

sub update_gene_position {
  my ($self, $dba, $dbid, $shift_by) = @_;
  
  my $sql = "
    UPDATE
      gene g INNER JOIN
      transcript t USING (gene_id) INNER JOIN
      exon_transcript et USING (transcript_id) INNER JOIN
      exon e USING (exon_id)
    SET
      g.seq_region_start = cast(g.seq_region_start as signed) + cast($shift_by as signed),
      g.seq_region_end   = cast(g.seq_region_end as signed)   + cast($shift_by as signed),
      t.seq_region_start = cast(t.seq_region_start as signed) + cast($shift_by as signed),
      t.seq_region_end   = cast(t.seq_region_end as signed)   + cast($shift_by as signed),
      e.seq_region_start = cast(e.seq_region_start as signed) + cast($shift_by as signed),
      e.seq_region_end   = cast(e.seq_region_end as signed)   + cast($shift_by as signed)
    WHERE gene_id = ?;
  ";
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($dbid) or $self->log_throw("Failed to execute: $sql");
}

1;
