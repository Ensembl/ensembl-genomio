=head1 LICENSE

Copyright [1999-2019] EMBL-European Bioinformatics Institute
and Wellcome Trust Sanger Institute

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

Bio::EnsEMBL::Pipeline::Runnable::EG:LoadGFF3::LoadGFF3

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::LoadGFF3;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base');

use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::PredictionExon;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

my $biotype_map = {
  gene => {},
  transcript => {
    lnc_RNA => 'lncRNA',
  },
};

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    gene_types      => ['gene', 'pseudogene', 'miRNA_gene', 'ncRNA_gene',
                        'rRNA_gene', 'snoRNA_gene', 'snRNA_gene', 'tRNA_gene',
                        'transposable_element'],
    mrna_types      => ['mRNA', 'transcript', 'pseudogenic_transcript',
                        'pseudogenic_rRNA', 'pseudogenic_tRNA',
                        'ncRNA', 'lincRNA', 'lncRNA', 'miRNA', 'pre_miRNA',
                        'RNase_MRP_RNA', 'RNAse_P_RNA', 'rRNA', 'snoRNA',
                        'snRNA', 'sRNA', 'SRP_RNA', 'tRNA', 'scRNA', 'guide_RNA',
                        'telomerase_RNA', 'antisense_RNA',
                        'transposable_element'],
    exon_types      => ['exon', 'pseudogenic_exon'],
    cds_types       => ['CDS'],
    utr_types       => ['five_prime_UTR', 'three_prime_UTR'],
    ignore_types    => ['misc_RNA', 'RNA',
                        'match', 'match_part',
                        'sequence_feature',
                        'cDNA_match', 'nucleotide_match', 'protein_match',
                        'polypeptide', 'protein',
                        'chromosome', 'supercontig', 'contig',
                        'region', 'biological_region',
                        'regulatory_region', 'repeat_region'],
    types_complete  => 1,
    use_name_field  => undef,
    polypeptides    => 1,
    min_intron_size => undef,
    nontranslating  => 'nontranslating_CDS',
    load_pseudogene_with_CDS => 0,
    prediction      => 0,
    gene_source     => 'Ensembl',
    stable_ids      => {},
  };
}

sub run {
  my ($self) = @_;
  my $gff3_file  = $self->param_required('gff3_file');
  my $fasta_file = $self->param_required('fasta_file');
  
  # Load the genes into an in-memory BioPerl database.
  my $db = $self->load_db($gff3_file, $fasta_file);
  $self->check_db($db);
  
  # Convert in-memory genes into EnsEMBL objects, and save to database.
  $self->load_genes($db);

  $self->dump_log();
}

sub load_db {
  my ($self, $gff_file, $fasta_file) = @_;
  
  my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
  
  my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store => $db, -ignore_seqregion => 1, -chunk_size => 10000000);
  
  $loader->load($gff_file, $fasta_file);
  
  return $db;
}

sub check_db {
  my ($self, $db) = @_;
  my @gene_types     = @{ $self->param_required('gene_types')   };
  my @mrna_types     = @{ $self->param_required('mrna_types')   };
  my @exon_types     = @{ $self->param_required('exon_types')   };
  my @cds_types      = @{ $self->param_required('cds_types')    };
  my @utr_types      = @{ $self->param_required('utr_types')    };
  my @ignore_types   = @{ $self->param_required('ignore_types') };
  my $types_complete = $self->param_required('types_complete');
  
  my %all_types =
    map {lc($_) => 1}
    (@gene_types, @mrna_types, @exon_types, @cds_types, @utr_types, @ignore_types);
  
  my @type_list = $db->types;
  my @types = map {$_->method} @type_list;
  my %types = map {$_ => 1} @types;
  
  my $info = "Loaded features to in-memory database:\n";
  
  my @unrecognised_types;
  foreach my $type (sort keys %types) {
    my @features = $db->get_features_by_type($type);
    $info .= "$type: ".scalar(@features)."\n";
    if (! exists $all_types{$type}) {
      push @unrecognised_types, $type;
    }
  }
  $self->log_warning($info);
  
  if (scalar(@unrecognised_types)) {
    my $msg = "Unrecognised types in GFF3 file: ".join(', ', @unrecognised_types);
    if ($types_complete) {
      $self->log_throw($msg);
    } else {
      $self->log_warning($msg);
    }
  }
}

sub load_genes {
  my ($self, $db) = @_;
  my $logic_name = $self->param_required('logic_name');
  my @gene_types = @{ $self->param_required('gene_types') };
  
  my $dba = $self->url2dba($self->param_required('db_url'));

  # Fetch slices and their synonyms into a lookup hash.
  my %slices = $self->fetch_slices($dba);
  
  # 'gff_genes' are the BioPerl gene objects.
  my @gff_genes = $db->get_features_by_type(@gene_types);
  
  # Only proceed if all genes are on scaffolds which exist in the db.
  $self->check_seq_ids(\%slices, \@gff_genes);
  
  # All genes and transcripts in the file are linked to a single analysis.
  my $aa       = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  if (! defined $analysis) {
    $self->log_throw("Analysis '$logic_name' does not exist in the database.");
  }
  
  # Instantiate adaptors for storing the data.
  my $ga  = $dba->get_adaptor('Gene');
  my $pta = $dba->get_adaptor('PredictionTranscript');
  
  foreach my $gff_gene (@gff_genes) {
    my $slice = $slices{$gff_gene->seq_id};
    
    my $gene = $self->new_gene($gff_gene, $slice, $analysis);
    
    $self->add_transcripts($db, $ga, $pta, $gff_gene, $gene);
    $self->set_pseudogene($ga, $gene);
  }

  $dba->dbc->disconnect_if_idle();
}

sub set_pseudogene {
  my ($self, $ga, $gene) = @_;

  # Only set the gene as pseudogene with CDS if all transcripts are pseudogenes and at least one has a CDS
  my $transcripts = $gene->get_all_Transcripts;
  my $num_tr = scalar @$transcripts;
  my $num_pseudo = scalar (grep { $_->biotype eq 'pseudogene' } @$transcripts);
  my $num_pseudo_CDS = scalar (grep { $_->biotype eq 'pseudogene_with_CDS' } @$transcripts);

  if ($num_tr eq ($num_pseudo + $num_pseudo_CDS)) {
    if ($num_pseudo_CDS > 0) {
      $self->log_warn("Set ".$gene->stable_id." as a pseudogene_with_CDS because all transcripts are pseudogenic and $num_pseudo_CDS/$num_tr are pseudogene_with_CDS\n");
      $gene->biotype('pseudogene_with_CDS');
      $ga->update($gene);
    } else {
      $self->log_warn("Set ".$gene->stable_id." as a pseudogene (without CDS) because all transcripts are pseudogenic\n");
      $gene->biotype('pseudogene');
      $ga->update($gene);
    }
  }
}

sub check_seq_ids {
  my ($self, $slices, $gff_objects) = @_;
  
  my %unrecognised;
  foreach my $gff_object (@$gff_objects) {
    if (! exists $$slices{$gff_object->seq_id}) {
      $unrecognised{$gff_object->seq_id}++;
    }
  }
  
  if (scalar(keys %unrecognised)) {
    my $msg = "Unrecognised sequences in GFF3 file: ".
              join(', ', sort keys %unrecognised);
    $self->log_throw($msg);
  }
}

sub add_transcripts {
  my ($self, $db, $ga, $pta, $gff_gene, $gene) = @_;
  my @mrna_types = @{ $self->param_required('mrna_types') };
  my $prediction = $self->param_required('prediction');
  
  my @gff_transcripts = $gff_gene->get_SeqFeatures(@mrna_types);
  if (scalar(@gff_transcripts) == 0) {
    my $transcript = $self->add_pseudogenic_transcript($db, $gff_gene, $gene);
    $gene->add_Transcript($transcript);
    $ga->store($gene);
  } else {
    if ($prediction) {
      foreach my $gff_transcript (@gff_transcripts) {
        my $transcript = $self->add_predicted_transcript($db, $gff_transcript, $gene);
        $pta->store($transcript);
      }
    } else {
      foreach my $gff_transcript (@gff_transcripts) {
        my $transcript = $self->add_transcript($db, $gff_transcript, $gene);

        if ($transcript->biotype eq 'protein_coding'
            or $transcript->biotype eq 'pseudogene_with_CDS') {
          $self->add_translation($db, $gff_transcript, $gene, $transcript);
        }
          
        $gene->add_Transcript($transcript);
      }
      
      $self->set_nontranslating_gene($gene);
    
      $ga->store($gene);
    }
  }
}

sub add_transcript {
  my ($self, $db, $gff_transcript, $gene) = @_;
  my @exon_types      = @{ $self->param_required('exon_types') };
  my $min_intron_size = $self->param('min_intron_size');
  
  my $transcript = $self->new_transcript($db, $gff_transcript, $gene);
  $transcript->stable_id($gene->stable_id) unless $transcript->stable_id;
  
  my @gff_exons = $gff_transcript->get_SeqFeatures(@exon_types);
  my $gff_exons;
  if (scalar(@gff_exons) == 0) {
    $gff_exons = $self->infer_exons($db, $gff_transcript);
  } else {
    $gff_exons = $self->merge_exons(\@gff_exons, $min_intron_size);
  }
  
  $self->add_exons($gff_exons, $transcript, 0);
  
  return $transcript;
}

sub add_predicted_transcript {
  my ($self, $db, $gff_transcript, $gene) = @_;
  my @exon_types      = @{ $self->param_required('exon_types') };
  my $min_intron_size = $self->param('min_intron_size');
  
  my $transcript = $self->new_predicted_transcript($gff_transcript, $gene);
  
  my @gff_exons = $gff_transcript->get_SeqFeatures(@exon_types);
  my $gff_exons;
  if (scalar(@gff_exons) == 0) {
    $gff_exons = $self->infer_exons($db, $gff_transcript);
  } else {
    $gff_exons = $self->merge_exons(\@gff_exons, $min_intron_size);
  }
  $self->add_exons($gff_exons, $transcript, 1);
  
  return $transcript;
}

sub add_pseudogenic_transcript {
  my ($self, $db, $gff_gene, $gene) = @_;
  my @exon_types      = @{ $self->param_required('exon_types') };
  my $min_intron_size = $self->param('min_intron_size');
  
  my $gff_transcript = $self->infer_transcript($db, $gff_gene);
  
  my $transcript = $self->new_transcript($db, $gff_transcript, $gene);
  $transcript->stable_id($gene->stable_id) unless $transcript->stable_id;
  
  my @gff_exons = $gff_gene->get_SeqFeatures(@exon_types);
  my $gff_exons;
  if (scalar(@gff_exons) == 0) {
    $gff_exons = $self->infer_exons($db, $gff_transcript);
  } else {
    $gff_exons = $self->merge_exons(\@gff_exons, $min_intron_size);
  }
  $self->add_exons($gff_exons, $transcript, 0);
  
  return $transcript;
}

sub add_exons {
  my ($self, $gff_exons, $transcript, $prediction) = @_;
  
  my $n = 1;
  
  foreach my $gff_exon (sort sort_coding @$gff_exons) {
    my $exon_id = $transcript->stable_id.'-E'.$n++;
    my $exon;
    if ($prediction) {
      $exon = $self->new_predicted_exon($gff_exon, $transcript);
    } else {
      $exon = $self->new_exon($gff_exon, $transcript, $exon_id);
    }
    $transcript->add_Exon($exon);
  }
}

sub add_translation {
  my ($self, $db, $gff_transcript, $gene, $transcript) = @_;
  my $polypeptides = $self->param_required('polypeptides');

  my ($translation_id, $gff_object, $genomic_start, $genomic_end);

  if ($polypeptides) {
    ($translation_id, $gff_object, $genomic_start, $genomic_end) = $self->get_polypeptide($db, $transcript->stable_id);
  }

  if (! defined $translation_id) {
    ($translation_id, $gff_object, $genomic_start, $genomic_end) = $self->infer_translation($gff_transcript, $transcript);
  }

  if (defined $genomic_start && defined $genomic_end && ($genomic_end - $genomic_start) > 1 ) {
    my ($start_exon, $end_exon, $seq_start, $seq_end) = $self->translation_coordinates($transcript, $genomic_start, $genomic_end);

    my $translation = $self->new_translation($translation_id, $start_exon, $end_exon, $seq_start, $seq_end);

    $self->add_xrefs($gff_object, $transcript->analysis, $translation, 'translation');

    $transcript->translation($translation);

    $self->set_exon_phase($transcript);

    if (!$transcript->translate()) {
      $self->log_warn("no translation for transcript ", $transcript->stable_id, "\n");
      $self->set_nontranslating_transcript($transcript);
    } else {
      my $seq = $transcript->translate()->seq;

      my $biotype = $transcript->biotype;

      if ($biotype eq 'pseudogene_with_CDS') {
        # Truncate the sequence if it has stop codons
        # Only keep whatever is before the first stop codon
        # (it might be an empty string, in which case the translation is empty)
        if ($seq and $seq =~ /\*/) {
          $self->log_warn("Pseudogene_with_CDS has stop codons: truncating (", $transcript->stable_id, ")\n");
          my ($seq) = split(/\*/, $seq);
        }

        # Change to normal pseudogene if there is no sequence
        if (not $seq or $seq  eq '') {
          $self->log_warn("changing biotype from pseudogene_with_CDS to pseudogene for", $transcript->stable_id, "\n");
          $transcript->biotype("pseudogene");
          $transcript->translation(undef);
        }
      } elsif (not $seq or $seq eq '' or $seq =~ /\*/) {
        $self->log_warn("no translation seq or one with the stop for transcript ", $transcript->stable_id,
                        ": $seq from ", $transcript->translateable_seq(),"\n");
        $self->set_nontranslating_transcript($transcript);
      }
    }
  } else {
    $self->log_warn("no translation (genomic start/end) for transcript ", $transcript->stable_id, "\n");
    $self->set_nontranslating_transcript($transcript);
  }
}

sub get_stable_id {
  my ($self, $gff_object) = @_;
  my $use_name_field = $self->param('use_name_field');
  my $stable_ids     = $self->param_required('stable_ids');
  
  my $stable_id = $gff_object->load_id;
  if (defined $use_name_field && $use_name_field eq 'stable_id') {
    $stable_id = $gff_object->name if $gff_object->name;
    
    # Because the provided Name won't necessarily be unique,
    # need to ensure uniqueness ourselves...
    if (exists $$stable_ids{$stable_id}) {
      $stable_id .= '_'.$$stable_ids{$stable_id}++;
      $self->log_warning("Added suffix to make stable ID unique: $stable_id");
    } else {
      $$stable_ids{$stable_id} = 1;
    }
  }
  
  return $stable_id;
}

sub get_cds {
  my ($self, $gff_transcript) = @_;
  my @cds_types = @{ $self->param_required('cds_types') };
  
  my @cds = ();
  my @gff_cds = sort sort_genomic $gff_transcript->get_SeqFeatures(@cds_types);
  foreach my $gff_cds (@gff_cds) {
    if (defined $gff_cds->segments && scalar($gff_cds->segments) > 0) {
      push @cds, sort sort_genomic $gff_cds->segments;
    } else {
      if ($gff_cds->start != $gff_cds->end) {
        push @cds, $gff_cds;
      }
    }
  }
  
  return @cds;
}

sub get_cds_id {
  my ($self, $gff_transcript) = @_;
  my @cds_types = @{ $self->param_required('cds_types') };
  
  my ($translation_id, $gff_object);
  my @gff_cds = sort sort_genomic $gff_transcript->get_SeqFeatures(@cds_types);
  if (scalar(@gff_cds) == 1) {
    $gff_object = $gff_cds[0];
    $translation_id = $self->get_stable_id($gff_object);
  }
  return ($translation_id, $gff_object);
}

sub get_utr {
  my ($self, $gff_transcript) = @_;
  my @utr_types = @{ $self->param_required('utr_types') };
  
  my @gff_utr = sort sort_genomic $gff_transcript->get_SeqFeatures(@utr_types);
  
  return @gff_utr;
}

sub get_polypeptide {
  my ($self, $db, $transcript_id) = @_;
  
  my ($translation_id, $polypeptide, $genomic_start, $genomic_end);
  
  my @polypeptides = $db->get_features_by_attribute(Derives_from => $transcript_id);
  if (scalar(@polypeptides)) {
    if (scalar(@polypeptides) > 1) {
      $self->log_warning("More than one polypeptide defined for $transcript_id. ".
                     "Only the longest will be processed.");
      my @sorted = sort { length($b->seq->seq) <=> length($a->seq->seq) } @polypeptides;
      $polypeptide = $sorted[0];
    } else {
      $polypeptide = $polypeptides[0];
    }
    
    $translation_id = $self->get_stable_id($polypeptide);
    $genomic_start  = $polypeptide->start;
    $genomic_end    = $polypeptide->end;
  }
  
  return ($translation_id, $polypeptide, $genomic_start, $genomic_end);
}

sub infer_transcript {
  my ($self, $db, $gff_gene) = @_;
  
  # If we have no mRNA/transcript, then we can assume that the transcript
  # has the same coordinates as the gene.
  my $gff_transcript = $db->new_feature(
    -primary_tag => 'transcript',
    -seq_id      => $gff_gene->seq_id,
    -start       => $gff_gene->start,
    -end         => $gff_gene->end,
    -strand      => $gff_gene->strand,
  );
  $gff_gene->add_SeqFeature($gff_transcript);
  
  return $gff_transcript;
}

sub infer_exons {
  my ($self, $db, $gff_transcript) = @_;
  
  # If we have no exons and no CDSs, then we can assume that the exon
  # has the same coordinates as the transcript. If we have CDSs, then
  # create exons with the same co-ordinates as them. If we additionally
  # have UTRs, then add new exons or adjust the boundaries of the exons
  # accordingly.
  my @gff_exons;
  my @gff_cds = $self->get_cds($gff_transcript);
  my @gff_utr = $self->get_utr($gff_transcript);
  
  if (scalar(@gff_cds) == 0) {
    my $exon = $db->new_feature(
      -primary_tag => 'exon',
      -seq_id      => $gff_transcript->seq_id,
      -start       => $gff_transcript->start,
      -end         => $gff_transcript->end,
      -strand      => $gff_transcript->strand,
    );
    $gff_transcript->add_SeqFeature($exon);
    push @gff_exons, $exon;
  } else {
    foreach my $gff_cds (sort sort_genomic @gff_cds) {
      my $exon = $db->new_feature(
        -primary_tag => 'exon',
        -seq_id      => $gff_cds->seq_id,
        -start       => $gff_cds->start,
        -end         => $gff_cds->end,
        -strand      => $gff_cds->strand,
        -phase       => $gff_cds->phase,
      );
      $gff_transcript->add_SeqFeature($exon);
      push @gff_exons, $exon;
    }
  }
  
  my (@before, @after);
  foreach my $gff_utr (sort sort_genomic @gff_utr) {
    # Do this in an unsophisticated, exhautive/long-winded way to keep the logic simple.
    # Rely on everything having been sorted into genomic order, so don't worry about strand.
    if ($gff_exons[0]->start - $gff_utr->end == 1) {
      # UTR tacked onto front of CDS
      $gff_exons[0]->start($gff_utr->start);
    } elsif ($gff_exons[0]->start - $gff_utr->end > 1) {
      # UTR is an exon in front of CDS
      push(@before, $gff_utr);
    } elsif ($gff_utr->start - $gff_exons[-1]->end == 1) {
      # UTR tacked onto end of CDS
      $gff_exons[-1]->end($gff_utr->end);
    } elsif ($gff_utr->start - $gff_exons[-1]->end > 1) {
      # UTR is an exon after CDS
      push(@after, $gff_utr);
    }
  }
  unshift(@gff_exons, @before);
  push(@gff_exons, @after);
  
  return \@gff_exons;
}

sub infer_translation {
  my ($self, $gff_transcript, $transcript) = @_;
  
  my ($genomic_start, $genomic_end);
  
  my @gff_cds = $self->get_cds($gff_transcript);
  my ($translation_id, $gff_object) = $self->get_cds_id($gff_transcript);
  $translation_id = $transcript->stable_id unless defined $translation_id;
  
  if (scalar(@gff_cds) > 0) {
    $genomic_start = $gff_cds[0]->start;
    $genomic_end   = $gff_cds[-1]->end;
    
    if ($transcript->strand == -1) {
      $genomic_end   -= $gff_cds[-1]->phase if defined $gff_cds[-1]->phase;
    } else {
      $genomic_start += $gff_cds[0]->phase if defined $gff_cds[0]->phase;
    }
  }
  
  return ($translation_id, $gff_object, $genomic_start, $genomic_end);
}

sub merge_exons {
  my ($self, $gff_exons, $min_intron_size) = @_;
  
  # If one exon is wholly within another ignore it; if there is
  # overlap at the 3' end, merge the exon; because we sort from 5' to 3',
  # there is no need to consider 5' overlaps.
  # If min_intron_size is defined, then close any gaps between
  # exons that are less than that by merging exons.
  my @corrected;
  my $previous;
  
  $min_intron_size = 0 unless defined $min_intron_size;
  
  foreach my $gff_exon (sort sort_coding @$gff_exons) {
    if ($gff_exon->strand == -1) {
      if (defined $previous) {
        if ( $gff_exon->end >= ($previous - $min_intron_size) ) {
          # The current exon either overlaps, or shares a sufficiently
          # small intron, with the previous exon, so adjust the previous
          # exon accordingly, as long as the current exon isn't wholly
          # within the previous exon.
          if ($gff_exon->start < $previous) {
            $corrected[-1]->start($gff_exon->start);
            $previous = $gff_exon->start;
          }
        } else {
          # Current exon does not need to be merged.
          push @corrected, $gff_exon;
          $previous = $gff_exon->start;
        }
      } else {
        # First exon, nothing to do.
        push @corrected, $gff_exon;
        $previous = $gff_exon->start;
      }
      
    } else {
      if (defined $previous) {
        if ( $gff_exon->start <= ($previous + $min_intron_size) ) {
          # The current exon either overlaps, or shares a sufficiently
          # small intron, with the previous exon, so adjust the previous
          # exon accordingly, as long as the current exon isn't wholly
          # within the previous exon.
          if ($gff_exon->end > $previous) {
            $corrected[-1]->end($gff_exon->end);
            $previous = $gff_exon->end;
          }
        } else {
          # Current exon does not need to be merged.
          push @corrected, $gff_exon;
          $previous = $gff_exon->end;
        }
      } else {
        # First exon, nothing to do.
        push @corrected, $gff_exon;
        $previous = $gff_exon->end;
      }
    }
  }
  
  return \@corrected;
}

sub translation_coordinates {
  my ($self, $transcript, $genomic_start, $genomic_end) = @_;
  
  my ($start_exon, $end_exon, $seq_start, $seq_end);
  
  my @exons = @{ $transcript->get_all_Exons };
  
  $start_exon = $exons[0];
  $end_exon   = $exons[-1];
  $seq_start  = 1;
  $seq_end    = $exons[-1]->length;
  
  foreach my $exon (@exons) {
    if ($genomic_start >= $exon->start && $genomic_start <= $exon->end) {
      if ($exon->strand == -1) {
        $end_exon = $exon;
        $seq_end = $exon->end - $genomic_start + 1;
      } else {
        $start_exon = $exon;
        $seq_start = $genomic_start - $exon->start + 1;
      }
    }
    if ($genomic_end >= $exon->start && $genomic_end <= $exon->end) {
      if ($exon->strand == -1) {
        $start_exon = $exon;
        $seq_start = $exon->end - $genomic_end + 1;
      } else {
        $end_exon = $exon;
        $seq_end = $genomic_end - $exon->start + 1;
      }
    }
  }
  
  return ($start_exon, $end_exon, $seq_start, $seq_end);
}

sub set_exon_phase {
  my ($self, $transcript) = @_;
  
  my ($phase, $end_phase) = (undef, undef);
  my $exons = $transcript->get_all_Exons;
  my $translation = $transcript->translation;
  
  # The phase and end_phase have defaults of -1, so only need to change
  # these when dealing with coding regions. The exons are automatically
  # returned in 5' -> 3' order.
  my $start_exon = 1;
  foreach my $exon (@{$exons}) {
    if (defined $exon->coding_region_start($transcript)) {
      if ($start_exon) {
        if ($translation->start == 1) {
          $phase = 0;
          $end_phase = ($exon->length) % 3;
        } else {
          $phase = -1;
          my $offset = $translation->start - 1;
          $end_phase = ($exon->length - $offset) % 3;
        }
        $start_exon = 0;
      } else {
        $end_phase = ($phase + $exon->length) % 3;
      }
      $exon->phase($phase);
      $exon->end_phase($end_phase);
      $phase = $end_phase;
    }
  }
  
  # End phase of -1 is conditional on there being a 3' UTR.
  if ($translation->end < $translation->end_Exon->length) {
    $translation->end_Exon->end_phase(-1);
  }
}

sub set_nontranslating_transcript {
  my ($self, $transcript) = @_;
  my $nontranslating = $self->param_required('nontranslating');
  
  $transcript->biotype($nontranslating);

  $self->log_warn("setting transcript ", $transcript->stable_id, " nontranslating\n");
  
  if ($nontranslating ne 'nontranslating_CDS') {
    $transcript->translation(undef);
  }
}

sub set_nontranslating_gene {
  my ($self, $gene) = @_;
  my $nontranslating = $self->param_required('nontranslating');
  
  if ($nontranslating eq 'pseudogene') {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      $transcript->biotype('pseudogene');
    }
  } else {
    # If any of its transcripts are protein_coding,
    # we want the gene to be so too.
    my $protein_coding_transcript = 0;
    my $nontranslating_transcript = 0;
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      if ($transcript->biotype eq 'protein_coding') {
        $protein_coding_transcript = 1;
      }
      if ($transcript->biotype eq 'nontranslating_CDS') {
        $nontranslating_transcript = 1;
      }
    }
    if ($protein_coding_transcript) {
      $gene->biotype('protein_coding');
    } elsif ($nontranslating_transcript) {
      $self->log_warn("Set ".$gene->stable_id." biotype to nontranslating_CDS (from set_nontranslating_gene, there are 'nontranslating_CDS' transcripts)\n");
      $gene->biotype('nontranslating_CDS');
    }
  }
}

sub map_biotype_gene {
  my ($self, $biotype, $gene) = @_;
  
  my $biotype_attr = $self->get_feature_biotype($gene);
  $biotype = $biotype_attr if $biotype_attr;
  
  if (exists $biotype_map->{gene}->{$biotype}) {
    return $biotype_map->{gene}->{$biotype};
  } else {
    return $biotype;
  }
}

sub map_biotype_transcript {
  my ($self, $biotype, $transcript) = @_;
  
  my $biotype_attr = $self->get_feature_biotype($transcript);
  $biotype = $biotype_attr if $biotype_attr;
  
  if (exists $biotype_map->{transcript}->{$biotype}) {
    return $biotype_map->{transcript}->{$biotype};
  } else {
    return $biotype;
  }
}

sub get_feature_biotype {
  my ($self, $feature) = @_;
  
  # Use biotype from GFF3
  my $biotype = $feature->type;
  
  # Get biotype attribute
  my %attr = $feature->attributes();
  my $biotype_attr = $attr{'biotype'} ? $attr{biotype}[0] : undef;
  
  # Replace biotype with the more precise biotype attribute
  if ($biotype_attr) {
    my $known_biotypes = $self->get_db_biotypes();
    
    if ($known_biotypes->{$biotype_attr}) {
      return $biotype_attr;
    }
  }
}

sub get_db_biotypes {
  my ($self) = @_;
  
  if (not $self->{_known_biotypes}) {
    my $dba = $self->url2dba($self->param_required('db_url'));
    my $ba = $dba->get_adaptor('Biotype');
    my %biotypes = map { $_->name => 1 } @{$ba->fetch_all()};
    $self->{_known_biotypes} = \%biotypes;
  }
  
  return $self->{_known_biotypes};
}

sub new_gene {
  my ($self, $gff_gene, $slice, $analysis) = @_;
  my $source = $self->param_required('gene_source');
  
  my $stable_id = $self->get_stable_id($gff_gene);
  
  my %atts = $gff_gene->attributes;
  
  my $biotype;
  if (exists $atts{'pseudo'}) {
    $biotype = 'pseudogene';
  } elsif ($gff_gene->type =~ /^gene:*/) {
    $biotype = 'protein_coding';
  } elsif ($gff_gene->type =~ /^(\w+):*/) {
    ($biotype) = $1;
  }
  $biotype = $self->map_biotype_gene($biotype, $gff_gene);
  
  my $gene = Bio::EnsEMBL::Gene->new(
    -stable_id     => $stable_id,
    -biotype       => $biotype,
    -source        => $source,
    -slice         => $slice,
    -start         => $gff_gene->start,
    -end           => $gff_gene->end,
    -strand        => $gff_gene->strand,
    -created_date  => time,
    -modified_date => time,
    -analysis      => $analysis,
  );
  $gene->version(undef);
  
  $self->add_xrefs($gff_gene, $analysis, $gene, 'gene');
  
  return $gene;
}

sub new_transcript {
  my ($self, $db, $gff_transcript, $gene) = @_;
  
  my $stable_id = $self->get_stable_id($gff_transcript);

  # Decide the biotype of the transcript
  my $biotype;
  my $gene_type = $gene->biotype;
  my $transcript_type = $gff_transcript->type;
  $transcript_type =~ s/:.+$//;
  
  # Pseudogene
  if ($gene_type eq 'pseudogene' or $transcript_type =~ /^pseudogenic_/) {
    my $tr_type = $transcript_type;
    $tr_type =~ s/^pseudogenic_//;
    
    # Protein_coding pseudogene: CDS or not?
    if ($tr_type eq "transcript" or $tr_type eq "mRNA") {
      # Check if there is a translation
      my ($translation_id) = $self->get_cds_id($gff_transcript);
      my $translatable = defined $translation_id;

      if ($translatable and $self->param('load_pseudogene_with_CDS')) {
        $self->log_warn("Pseudogene has CDSs: $stable_id\n");
        $biotype = 'pseudogene_with_CDS';
      } else {
        $self->log_warn("Pseudogene has no CDS: $stable_id ($transcript_type, $gene_type)\n");
        $biotype = 'pseudogene';
      }
      
    # Pseudogenic_tRNA
    } elsif ($tr_type eq "tRNA") {
        $self->log_warn("Pseudogenic tRNA: $stable_id\n");
        $biotype = 'tRNA_pseudogene';
        $gene->biotype($biotype);

    # Pseudogenic_rRNA
    } elsif ($tr_type eq "pseudogenic_rRNA") {
        $self->log_warn("Pseudogenic rRNA: $stable_id\n");
        $biotype = 'rRNA_pseudogene';
        $gene->biotype($biotype);
    
    # Other kinds of pseudogenes?
    } else {
      die("Unrecognized pseudogene biotype: $transcript_type (gene: $gene_type)");
    }
    
  # transposable_element
  } elsif ($gene_type eq 'transposable_element') {
    $biotype = "transposable_element";
    $stable_id = $gene->stable_id . '-RA' if (!$stable_id && $gene->stable_id);
    $self->log_warn("transposable_element: $stable_id\n");
  
  # Protein coding
  } elsif ($transcript_type eq 'mRNA' or $transcript_type eq "transcript") {
    $biotype = "protein_coding";
    $self->log_warn("protein_codig: $stable_id\n");
  
  # Non protein coding
  } else {
    $biotype = $transcript_type;
    $biotype = $self->map_biotype_transcript($biotype, $gff_transcript);
    $self->log_warn("non-protein_codig: $stable_id\n");

    # NB: we may need a control of what biotypes are known or not
    $self->log_warn("seting biotype to $biotype for non-codig: $stable_id (gene ", $gene->stable_id, ")\n");
    $gene->biotype($biotype);
  }

  my $transcript = Bio::EnsEMBL::Transcript->new(
    -stable_id     => $stable_id,
    -biotype       => $biotype,
    -source        => $gene->source,
    -slice         => $gene->slice,
    -start         => $gff_transcript->start,
    -end           => $gff_transcript->end,
    -strand        => $gff_transcript->strand,
    -created_date  => time,
    -modified_date => time,
    -analysis      => $gene->analysis,
  );
  $transcript->version(undef);
  
  $self->add_xrefs($gff_transcript, $transcript->analysis, $transcript, 'transcript');
  
  return $transcript;
}

sub new_predicted_transcript {
  my ($self, $gff_transcript, $gene) = @_;
  
  my $display_label = $self->get_stable_id($gff_transcript);
  
  my $transcript = Bio::EnsEMBL::PredictionTranscript->new(
    -display_label => $display_label,
    -slice         => $gene->slice,
    -start         => $gff_transcript->start,
    -end           => $gff_transcript->end,
    -strand        => $gff_transcript->strand,
    -analysis      => $gene->analysis,
  );
  
  return $transcript;
}

sub new_exon {
  my ($self, $gff_exon, $transcript, $exon_id) = @_;
  
  my $exon = Bio::EnsEMBL::Exon->new(
    -stable_id     => $exon_id,
    -slice         => $transcript->slice,
    -start         => $gff_exon->start,
    -end           => $gff_exon->end,
    -strand        => $gff_exon->strand,
    -phase         => -1,
    -end_phase     => -1,
    -created_date  => time,
    -modified_date => time,
  );
  $exon->version(undef);
  
  return $exon;
}

sub new_predicted_exon {
  my ($self, $gff_exon, $transcript) = @_;
  
  my %atts    = $gff_exon->attributes;
  my $att     = $atts{'p_value'};
  my $p_value = $$att[0];
  
  my $exon = Bio::EnsEMBL::PredictionExon->new(
    -slice   => $transcript->slice,
    -start   => $gff_exon->start,
    -end     => $gff_exon->end,
    -strand  => $gff_exon->strand,
    -phase   => -1,
    -score   => $gff_exon->score,
    -p_value => $p_value,
  );
  
  return $exon;
}

sub new_translation {
  my ($self, $translation_id, $start_exon, $end_exon, $seq_start, $seq_end) = @_;
  
  my $translation = Bio::EnsEMBL::Translation->new(
    -stable_id     => $translation_id,
    -start_exon    => $start_exon,
    -end_exon      => $end_exon,
    -seq_start     => $seq_start,
    -seq_end       => $seq_end,
    -created_date  => time,
    -modified_date => time,
  );
  $translation->version(undef);
  
  return $translation;
}

sub add_xrefs {
  my ($self, $gff_object, $analysis, $object, $object_type) = @_;
  my $use_name_field = $self->param('use_name_field');
  
  if (defined $use_name_field && $use_name_field eq 'xref') {
    my $name = $gff_object->name;
    if ($name) {
      my $edb_name = $self->param_required("xref_$object_type\_external_db");
      
      my $xref = $self->new_xref($edb_name, $name, $name, $analysis);
      
      $object->add_DBEntry($xref);
    }
  }
}

sub new_xref {
  my ($self, $external_db_name, $primary_id, $display_id, $analysis) = @_;
  
  my $xref = Bio::EnsEMBL::DBEntry->new
  (
    -dbname      => $external_db_name,
    -primary_id  => $primary_id,
    -display_id  => $display_id,
  );
  $xref->analysis($analysis);
  
  return $xref;
}

sub sort_coding {  
  if ($a->strand == -1) {
    return $b->start <=> $a->start;
  } else {
    return $a->start <=> $b->start;
  }
}

sub sort_genomic {  
  return $a->start <=> $b->start;
}

1;
