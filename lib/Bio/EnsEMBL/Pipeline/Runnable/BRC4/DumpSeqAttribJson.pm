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
       my %attrib = (
          object_type => $object_type,
          object_id   => $transcript->stable_id,
          seq_attrib_type => $allowed_attribs->{$code},
        );
      my $attrib_values;

      # Get attrib specific values
      if ($code eq '_rna_edit') {
        $attrib_values = $self->format_edit($attrib);
      } elsif ($code eq '_transl_start' or
               $code eq '_transl_end') {
        $attrib_values = $self->format_position($attrib);
      } elsif ($code eq 'Frameshift') {
        $attrib_values = $self->format_frameshift($attrib);
      }
      die("Could not get attrib values for $code, ".$transcript->stable_id) if not $attrib_values;

      # Merge attrib values
      %attrib = (%attrib, %$attrib_values);
      push \@selected_attribs, \%attrib;
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
       my %attrib = (
          object_type => $object_type,
          object_id   => $translation->stable_id,
          seq_attrib_type => $allowed_attribs->{$code},
        );
      my $attrib_values;

      # Get attrib specific values
      if ($code eq 'amino_acid_sub') {
        $attrib_values = $self->format_edit($attrib);
      } elsif ($code eq '_selenocysteine') {
        $attrib_values = $self->format_position($attrib);
      }
      die("Could not get attrib values for $code, ".$translation->stable_id) if not $attrib_values;

      # Merge attrib values
      %attrib = (%attrib, %$attrib_values);
      push \@selected_attribs, \%attrib;
    }
  }
  return \@selected_attribs;
}

sub format_edit {
  my ($self, $attrib) = @_;

  my ($start, $end, $seq) = split / /, $attrib->value;
  my %values = (
    start => int($start),
    end => int($end),
    sequence => $seq,
  );
  return \%values;
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
