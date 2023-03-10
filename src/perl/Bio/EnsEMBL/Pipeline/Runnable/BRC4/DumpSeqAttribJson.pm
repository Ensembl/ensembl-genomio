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


package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpSeqAttribJson;

use strict;
use warnings;

use JSON;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpJsonBase');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    metadata_name => 'seq_attrib',
    allowed_transcript_attribs => {
      _rna_edit     => "sequence_alteration",
      Frameshift    => "frameshift",
      _transl_start => "coding_start",
      _transl_end   => "coding_end",
    },
    allowed_translation_attribs => {
      amino_acid_sub  => "sequence_alteration",
      initial_met     => "sequence_alteration",
      _selenocysteine => "selenocysteine",
    },
  };
}

sub prepare_data {
  my ($self) = @_;

  my $dba = $self->core_dba();

  # Get genes
  my @features;
  my $ta = $dba->get_adaptor('Transcript');
  my $pa = $dba->get_adaptor('Translation');

  my @items;
  push @items, @{$ta->fetch_all()};
  push @items, @{$pa->fetch_all()};
  
  my @seq_attribs;

  # Get transcript attribs
  my $type = 'transcript';
  foreach my $transcript (@{$ta->fetch_all()}) {
    my $attribs = $self->get_transcript_seq_attribs($transcript);
    push @seq_attribs, @$attribs;
  }
  
  # Get translation attribs
  $type = 'translation';
  foreach my $translation (@{$pa->fetch_all()}) {
    my $attribs = $self->get_translation_seq_attribs($translation);
    push @seq_attribs, @$attribs;
  }

  $dba->dbc()->disconnect_if_idle();

  # Sort for easier file comparison
  @seq_attribs = sort { $a->{object_id} cmp $b->{object_id} } @seq_attribs;

  return \@seq_attribs;
}


sub get_transcript_seq_attribs {
  my ($self, $transcript) = @_;

  my $allowed_attribs = $self->param('allowed_transcript_attribs');

  # Get all attributes at once
  my $attribs = $transcript->get_all_Attributes();

  my @selected_attribs;
  my $object_type = 'transcript';
  for my $attrib (@$attribs) {
    my $code = $attrib->code;
    if ($allowed_attribs->{$code}) {
      $code = $allowed_attribs->{$code};
      my %attrib = (
        object_type => $object_type,
        object_id   => $transcript->stable_id,
        seq_attrib_type => $code,
      );
      my $attrib_values;

      # Get attrib specific values
      if ($code eq 'sequence_alteration') {
        $attrib_values = $self->format_edit($attrib);
      } elsif ($code eq 'coding_start' or
               $code eq 'coding_end') {
        $attrib_values = $self->format_position($attrib);
      } elsif ($code eq 'frameshift') {
        $attrib_values = $self->format_frameshift($attrib);
      }
      die("Could not get attrib values for $code, ".$transcript->stable_id) if not $attrib_values;

      # Merge attrib values
      %attrib = (%attrib, %$attrib_values);
      push @selected_attribs, \%attrib;
    }
  }
  return \@selected_attribs;
}

sub get_translation_seq_attribs {
  my ($self, $translation) = @_;

  my $allowed_attribs = $self->param('allowed_translation_attribs');

  # Get all attributes at once
  my $attribs = $translation->get_all_Attributes();

  my @selected_attribs;
  my $object_type = 'translation';
  for my $attrib (@$attribs) {
    my $code = $attrib->code;
    if ($allowed_attribs->{$code}) {
      $code = $allowed_attribs->{$code};
      my %attrib = (
        object_type => $object_type,
        object_id   => $translation->stable_id,
        seq_attrib_type => $code,
      );
      my $attrib_values;

      # Get attrib specific values
      if ($code eq 'sequence_alteration') {
        $attrib_values = $self->format_edit($attrib);
      } elsif ($code eq 'selenocysteine') {
        $attrib_values = $self->format_position($attrib);
      }
      die("Could not get attrib values for $code, ".$translation->stable_id) if not $attrib_values;

      # Merge attrib values
      %attrib = (%attrib, %$attrib_values);
      push @selected_attribs, \%attrib;
    }
  }
  return \@selected_attribs;
}

sub format_edit {
  my ($self, $attrib) = @_;

  my ($start, $end, $seq) = split / /, $attrib->value;
  ($start, $end) = $self->ensembl_to_interbase($start, $end);

  my %values = (
    start => int($start),
    end => int($end),
    sequence => $seq,
  );
  return \%values;
}

sub ensembl_to_interbase {
  my ($self, $start, $end) = @_;

  # The only difference between Ensembl current system (1-based) and
  # an interbase system (0-based) is that Ensembl start is +1
  # The result is more legible
  $start--;

  return ($start, $end);
}

sub format_position {
  my ($self, $attrib) = @_;

  my %values = (
    position => int($attrib->value),
  );

  return \%values;
}

sub format_frameshift {
  my ($self, $attrib) = @_;

  my %values = (
    intron_number => int($attrib->value),
  );

  return \%values;
}

1;
