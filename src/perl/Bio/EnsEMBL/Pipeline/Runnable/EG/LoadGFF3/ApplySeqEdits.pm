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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ApplySeqEdits

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ApplySeqEdits;

use strict;
use warnings;
use feature 'say';

use Path::Tiny qw(path);

use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

use base ('Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base');

sub run {
  my ($self) = @_;
  my $logic_name         = $self->param_required('logic_name');
  my $genbank_file       = $self->param('genbank_file');
  my $protein_fasta_file = $self->param('protein_fasta_file');
  
  my $dba = $self->get_type_dba();
  
  if ($genbank_file) {
    $self->seq_edits_from_genbank($dba, $genbank_file);
  }
  
  if ($protein_fasta_file) {
    $self->seq_edits_from_protein($dba, $logic_name, $protein_fasta_file);
  }
  
  $self->set_protein_coding($dba, $logic_name);

  $dba->dbc->disconnect_if_idle();

  $self->dump_log();
}

sub seq_edits_from_genbank {
  my ($self, $dba, $genbank_file) = @_;
  
  my $ta  = $dba->get_adaptor('Transcript');
  my $aa  = $dba->get_adaptor("Attribute");
  
  my $genbank_path = path($genbank_file);
  my $genbank = $genbank_path->slurp;
  $genbank =~ s!\A//\s*!!m;
  $genbank =~ s!//\s*\Z!!m;
  my @genbank = split(m!//\s+!, $genbank);
  
  my ($edit_count, $total_count) = (0, 0);

  foreach my $record (@genbank) {
    if ($record =~ /##RefSeq-Attributes-START##.*(assembly gap|frameshifts)/ms) {
      
      # Probably will only have mRNA, but check just in case.
      my ($mol_type) = $record =~ /\s+\/mol_type="([^"]+)"/m;
      next if $mol_type ne 'mRNA';
      
      my ($accession) = $record =~ /^VERSION\s+(\S+)/m;
      my $transcript  = $ta->fetch_by_stable_id($accession);
      $self->log_throw("$accession in file, not in database") unless $transcript;
      #$self->log_warning($accession);
      
      my $cdna_seq = $self->extract_seq($record);
      
      # Check that the sequences don't already match before trying anything...
      if ($transcript->seq->seq ne $cdna_seq) {
        $self->log_warning('Attempting seq_edits for ' . $transcript->stable_id);
        $total_count++;
        
        # Start and end will only be defined if they fall within a seq_edit.
        my ($edits, $start, $end) = $self->extract_edits($record, $cdna_seq);
        
        my @attribs;
        
        if (defined $start) {
          my $start_attrib = $self->add_translation_start($transcript, $start);
          push @attribs, $start_attrib;
        }
        
        if (defined $end) {
          my $end_attrib = $self->add_translation_end($transcript, $end);
          push @attribs, $end_attrib;
        }
        
        foreach my $edit (@$edits) {
          my $seq_edit = $self->add_transcript_seq_edit($transcript, @$edit);
          $self->log_warning("Adding seq_edit for " . $transcript->stable_id . ": " . join(' ', @$edit));
          push @attribs, $seq_edit;
        }
        
        if ($transcript->seq->seq eq $cdna_seq) {
          $self->log_warning('Storing seq_edits for ' . $transcript->stable_id);
          $aa->store_on_Transcript($transcript, \@attribs);
          $edit_count++;
        } else {       
          $self->log_warning("Not storing seq_edits\nEdited seq:\n" . $transcript->seq->seq . "\nTarget seq:\n$cdna_seq\n");
          
          #my $seq_edit = $self->add_transcript_seq_edit($transcript, 1, $transcript->length, $cdna_seq);
          #$aa->store_on_Transcript($transcript, [$seq_edit]);
        }
      }
    }
  }
  
  $self->log_warning("$edit_count / $total_count seq_edits successful");
}

sub extract_edits {
  my ($self, $record, $cdna_seq) = @_;
  
  my %transcriptomes = ();
  my ($transcriptomes) = $record =~ /\s+transcript\s+sequences*\s+\((.*?)\)/ms;
  if ($transcriptomes) {
    my @transcriptomes = split(/,\s*|\s*and\s*/, $transcriptomes);
    %transcriptomes    = map { $_ => 1 } @transcriptomes;
  }
  
  my ($coords) = $record =~ /^PRIMARY[^\n]+\n(^\s+.+)/ms;
  $coords =~ s/\n^\S.+//ms;
  my @coords = split(/\n/, $coords);
  my @edits = ();
  
  # Need to detect if a coding region starts or ends in a seq_edit.
  my ($cds_start, $cds_end) = $record =~ /^\s+CDS\D+(\d+)\D+(\d+)\s*$/m;
  my ($transl_start, $transl_end);
  
  # Deletions are rather cryptic. Need to compare pairs of coords to look
  # for gaps of 1 or 2 bases, so need to keep track of last position.
  my ($last_seq_name, $last_seq_region_start, $last_seq_region_end);
  
  # The positions given in the file are all relative to the original
  # sequence. But in the Ensembl API, seq_edits are applied to the
  # result of any previous edits. So need to track what we've already
  # done, and modify positions accordingly.
  my $offset = 0;
  
  foreach my $coord (@coords) {
    my ($start, $end, $seq_name, $seq_region_start, $seq_region_end) =
      $coord =~ /^\s+(\d+)\-(\d+)\s+(\S+)\s+(\d+)\-(\d+)/;
    
    # Inserted sequence
    my $subseq;
    if (exists $transcriptomes{$seq_name}) {
      $subseq = substr($cdna_seq, $start - 1, $end-$start + 1);
      
      if ($cds_start >= $start && $cds_start <= $end) {
        $transl_start = $cds_start;
      }
      if ($cds_end >= $start && $cds_end <= $end) {
        $transl_end = $cds_end;
      }
      
    } elsif ($seq_name =~ /"(NN?)"/) {
      $subseq = $1;
    }
    if (defined $subseq) {
      push @edits, [$start - $offset, $start - $offset - 1, $subseq];
      $offset += length($subseq);
    }
    
    # Deleted sequence
    if (defined $last_seq_name && $last_seq_name eq $seq_name) {
      my $gap;
      if ($last_seq_region_end < $seq_region_start) {
        $gap = $seq_region_start - $last_seq_region_end - 1;
      } else {
        $gap = $last_seq_region_start - $seq_region_end - 1;
      }
      if ($gap <= 2) {
        push @edits, [$start - $offset, $start - $offset + $gap - 1, ''];
        $offset -= $gap
      }
    }
    
    $last_seq_name         = $seq_name;
    $last_seq_region_start = $seq_region_start;
    $last_seq_region_end   = $seq_region_end;
  }
  
  return (\@edits, $transl_start, $transl_end);
}

sub extract_seq {
  my ($self, $record) = @_;
  
  my ($seq) = $record =~ /^ORIGIN[^\n]+\n(.+)/ms;
  $seq =~ s/\n^\S.+//ms;
  $seq =~ s/\d+//gm;
  $seq =~ s/\s+//gm;
  
  return uc($seq);
}

sub seq_edits_from_protein {
  my ($self, $dba, $logic_name, $protein_fasta_file) = @_;
  
  my $sa = $dba->get_adaptor('Slice');
  my $ta = $dba->get_adaptor('Transcript');
  my $aa = $dba->get_adaptor("Attribute");
  
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
      my $translation = $transcript->translation;

      if ($translation) {
        # Get sequence from FASTA file if protein ID or transcript ID match
        my $file_seq;
        if ($protein{$translation->stable_id}) {
          $file_seq = $protein{$translation->stable_id};
        } elsif ($protein{$transcript->stable_id}) {
          $file_seq = $protein{$transcript->stable_id};
        } else {
          next;
        }
        $file_seq =~ s/\*$//;
        my $db_seq = $translation->seq;

        # Do not want to consider an amino acid derived from a partial
        # codon; but RefSeq do that, so need to lop off last amino acid
        # in that case.
        my $transcript_length = length($transcript->translateable_seq);
        my $partial_codon     = $transcript_length % 3;
        if ($partial_codon && length($file_seq) == length($db_seq) + 1) {
          $self->log_warning("WARNING: Loping off partial last codon for transcript "
              . $transcript->stable_id
              . " translation "
              . $translation->stable_id);
          $file_seq =~ s/.$//;
        }

        if ($db_seq ne $file_seq) {
          my @file_seq = split(//, $file_seq);
          my @atts;

          if (length($file_seq) == length($db_seq)) {
            while ((my $pos = index($db_seq, '*')) >= 0) {
              my $amino_acid = $file_seq[$pos];
              $db_seq =~ s/\*/$amino_acid/;

              $pos += 1;
              my $att = $self->add_translation_seq_edit($translation, $pos, $pos, $amino_acid);
              push @atts, $att;

              $self->log_warning(
                    "WARNING: replacing stop codon at $pos with $amino_acid  for transcript "
                  . $transcript->stable_id
                  . " translation "
                  . $translation->stable_id);
            }

            if ($db_seq eq $file_seq) {
              $aa->store_on_Translation($translation, \@atts);
            }
          } else {
            $self->log_warning(
              'WARNING: Protein sequence length mismatch for ' . $translation->stable_id);
          }
        }
      }
    }
  }

}

sub add_transcript_seq_edit {
  my ($self, $transcript, $start, $end, $seq) = @_;
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => '_rna_edit',
    -VALUE => "$start $end $seq",
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $transcript->add_Attributes($attribute);
  
  return $attribute;
}

sub add_translation_seq_edit {
  my ($self, $translation, $start, $end, $seq) = @_;
  
  my $code = 'amino_acid_sub';
  if ($seq eq 'U') {
    $code = '_selenocysteine';
  }
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => $code,
    -VALUE => "$start $end $seq",
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $translation->add_Attributes($attribute);
  
  return $attribute;
}

sub add_translation_start {
  my ($self, $transcript, $start) = @_;
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => '_transl_start',
    -VALUE => $start,
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $transcript->add_Attributes($attribute);
  
  return $attribute;
}

sub add_translation_end {
  my ($self, $transcript, $end) = @_;
  
  my $attribute = Bio::EnsEMBL::Attribute->new
  (
    -CODE  => '_transl_end',
    -VALUE => $end,
  );
  
  # Need to add the attribute in order to force recalculation
  # of the sequence for the benefit of any subsequent edits.
  $transcript->add_Attributes($attribute);
  
  return $attribute;
}

1;
