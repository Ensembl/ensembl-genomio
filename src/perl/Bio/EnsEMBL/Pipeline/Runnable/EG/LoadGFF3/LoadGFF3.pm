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

Bio::EnsEMBL::Pipeline::Runnable::EG:LoadGFF3::LoadGFF3

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::LoadGFF3;

use strict;
use warnings;
use feature 'say';
use List::Util qw/sum0/;

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
  gene => {
    lnc_RNA => 'lncRNA',
    #misc_RNA => 'ncRNA',
    ncRNA_gene => 'ncRNA',
    other      => undef, # let the 'sub new_transcript' set the proper one
    transcript => 'misc_RNA',
  },
  transcript => {
    lnc_RNA => 'lncRNA',
    mRNA => 'ncRNA',
    #misc_RNA => 'ncRNA',
    transcript => 'misc_RNA',
  },
};

# use to force transfer biotypes from transcripts to genes, see `sub new_transcript`
my $transfer_biotype_tr2gene = {};

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},

    # Lists of the types that we expect to see in the GFF3 file
    gene_types      => ['gene', 'pseudogene', 'miRNA_gene', 'ncRNA_gene',
                        'rRNA_gene', 'snoRNA_gene', 'snRNA_gene', 'tRNA_gene',
                        'transposable_element'],
    mrna_types      => ['mRNA', 'transcript', 'misc_RNA', 'RNA',
                        'pseudogenic_transcript', 'pseudogenic_rRNA', 'pseudogenic_tRNA',
                        'ncRNA', 'lincRNA', 'miRNA', 'pre_miRNA',
                        'lncRNA', 'lnc_RNA', 'piRNA',
                        'RNase_MRP_RNA', 'RNAse_P_RNA', 'rRNA', 'snoRNA',
                        'snRNA', 'sRNA', 'SRP_RNA', 'tRNA', 'scRNA', 'guide_RNA',
                        'tmRNA', 'scaRNA',
                        'telomerase_RNA', 'antisense_RNA',
                        'transposable_element',
                        'TR_V_gene','TR_C_gene',
                        'IG_V_gene','IG_C_gene',
                        'ribozyme',
			'antisense_RNA', 'miRNA_primary_transcript'],
    exon_types      => ['exon', 'pseudogenic_exon'],
    cds_types       => ['CDS'],
    utr_types       => ['five_prime_UTR', 'three_prime_UTR'],
    ignore_types    => ['match', 'match_part',
                        'sequence_feature',
                        'cDNA_match', 'nucleotide_match', 'protein_match',
                        'polypeptide', 'protein',
                        'chromosome', 'supercontig', 'contig',
                        'region', 'biological_region',
                        'regulatory_region', 'repeat_region',
                        'long_terminal_repeat', 'STS',
                        'D_loop','origin_of_replication',
			'non_canonical_five_prime_splice_site','non_canonical_three_prime_splice_site'],

    # By default, it is assumed that the above type lists are exhaustive.
    # If there is a type in the GFF3 that is not listed, an error will be
    # thrown, unless 'types_complete' = 0.
    types_complete  => 1,

    # By default, load the GFF3 "ID" fields as stable_ids, and ignore "Name"
    # fields. If they exist, can load them as stable IDs instead, with the
    # value 'stable_id'; or load them as xrefs by setting to 'xref'.
    use_name_field  => undef,

    # If there are polypeptide rows in the GFF3, defined by 'Derives_from'
    # relationships, those will be used to determine the translation
    # (rather than inferring from CDS), unless 'polypeptides' = 0.
    # N.B. if on, could lead to models with the missing stop codon
    polypeptides    => 1,

    # Some sources have 1- or 2-base introns
    # defined to deal with readthrough stop codons. But their sequence
    # files contradict this, and include those intronic bases. The NCBI
    # .gbff files then define a 1- or 2-base insertion of 'N's, which
    # fixes everything up (and which this pipeline can handle).
    # So, the pipeline can merge exons that are separated by a small intron.
    min_intron_size => undef,

    # Set the biotype for transcripts that produce invalid translations. We
    # can treat them as pseudogenes by setting 'nontranslating' => "pseudogene".
    # The default is "nontranslating_CDS", and in this case the translation
    # is still added to the core db, on the assumption that remedial action
    # will fix it (i.e. via the ApplySeqEdits module).
    nontranslating  => 'nontranslating_CDS',

    load_pseudogene_with_CDS => 0,

    prediction      => 0,

    gene_source     => 'Ensembl',
    stable_ids      => {},

    # use common prefix as the stable_id for CDS (within the same multifeature)
    find_multifeature_commmon_name => 1,
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
  
  # append source types from $biotype_map
  my @biotype_map_source_types = ();
  for my $feature_type (keys %$biotype_map) {
    push @biotype_map_source_types, keys %{$biotype_map->{$feature_type}};
  }
  @biotype_map_source_types = keys %{{ map {$_ => 1} @biotype_map_source_types }};

  my %all_types =
    map {lc($_) => 1}
    (@gene_types, @mrna_types, @exon_types, @cds_types, @utr_types, @ignore_types, @biotype_map_source_types);
  
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

  my $dba = $self->get_type_dba();

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
    $self->update_gene_coords_based_on_transcripts_on_circular_sr($gene, $ga);
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
      $self->log_warning("Set " . $gene->stable_id . " as a pseudogene_with_CDS because all transcripts are pseudogenic and $num_pseudo_CDS/$num_tr are pseudogene_with_CDS");
      $gene->biotype('pseudogene_with_CDS');
      $ga->update($gene);
    } else {
      $self->log_warning("Set " . $gene->stable_id . " as a pseudogene (without CDS) because all transcripts are pseudogenic");
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
            or $transcript->biotype eq 'IG_V_gene'
            or $transcript->biotype eq 'TR_V_gene'
            or $transcript->biotype eq 'IG_C_gene'
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

sub update_gene_coords_based_on_transcripts_on_circular_sr {
  my ($self, $gene, $ga) = @_;

  # fixing only on circular seq_regions
  return if (!$gene->slice->is_circular);
  
  # we can assume that we work with the top-level seq_regions,
  #   thus Gene::start/end and Gene::seq_region_start/end are equivalent
  my ($gene_start, $gene_end) = ( $gene->start, $gene->end );
  $gene->recalculate_coordinates();
  if ($gene_start != $gene->start || $gene_end != $gene->end) {
    # ideally, just call
    #   $ga->update_coords($gene);
    # unfortunately this doesn't work for circular seq_regions
    # when update_coords is called it uses Gene::seq_region_start/end, but
    # Gene::seq_region_start/end uses BaseFeatureAdaptor::_seq_region_boundary_from_db
    # failback mechanism for circular regions and ignore Feature:: coords being updated
    # (see https://github.com/Ensembl/ensembl/blob/release/111/modules/Bio/EnsEMBL/Feature.pm#L1067)

    # going vanila SQL update way
    $self->update_gene_seq_region_coords_from_coords($gene, $ga);
  }
}


sub update_gene_seq_region_coords_from_coords {
  # function to update gene coordinates in db based on the object state
  # ideally,
  #   $ga->update_coords($gene);
  # should be called, but unfortunately this doesn't work for _circular_ seq_regions
  #   when update_coords is called it uses Gene::seq_region_start/end, but
  #   Gene::seq_region_start/end uses BaseFeatureAdaptor::_seq_region_boundary_from_db
  #   failback mechanism for circular regions and ignore Feature:: coords being updated
  #   (see https://github.com/Ensembl/ensembl/blob/release/111/modules/Bio/EnsEMBL/Feature.pm#L1067)
  # for _circular_ regions going vanila SQL UPDATE way

  my ($self, $gene, $ga) = @_;

  return if !$gene;
  
  # if not circular use the default mechanism
  if (!$gene->slice->is_circular) {
    $ga->update_coords($gene);
    return;
  }

  # for circular
  #   (copied from Bio::EnsEMBL::DBSQL::GeneAdaptor::update_coords)
  my $update_sql = qq(
    UPDATE gene
       SET seq_region_start = ?,
           seq_region_end = ?
       WHERE gene_id = ?
    );
  
  # we can assume that we work with the top-level seq_regions,
  #   thus Gene::start/end and Gene::seq_region_start/end are equivalent
  my $sth = $ga->prepare($update_sql);
  $sth->bind_param(1, $gene->start); # ~ seq_region_start
  $sth->bind_param(2, $gene->end); # ~ seq_region_end
  $sth->bind_param(3, $gene->dbID);
  $sth->execute();
}

sub add_exons {
  my ($self, $gff_exons, $transcript, $prediction) = @_;
  
  # force exon ranking for strange cases like transplacing or for coords > seq_region length (on circular)
  my @outliers = grep { $_ ->start > $transcript->slice->length } @$gff_exons;
  my @strands = map {$_->strand} @$gff_exons;
  my $different_strands = scalar(keys %{{ map { $_=>1 } @strands }});
  my $force_ranking = ($different_strands > 1 || scalar(@outliers) > 0) ? 1 : 0;

  # Use the strandness majority to decide in which order the exons of transspliced genes are ranked
  my $strandness = sum0 @strands;
  my $sorting = ($strandness > 0)? \&sort_genomic : \&sort_genomic_desc;

  my $n = 1;
  foreach my $gff_exon (sort $sorting @$gff_exons) {
    my $exon_id = $transcript->stable_id.'-E'.$n++;
    my $exon;
    if ($prediction) {
      $exon = $self->new_predicted_exon($gff_exon, $transcript);
    } else {
      $exon = $self->new_exon($gff_exon, $transcript, $exon_id);
    }
    $transcript->add_Exon($exon, ($force_ranking ? $n-1 : undef));
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
      $self->log_warning("WARNING: no translation for transcript " . $transcript->stable_id);
      $self->set_nontranslating_transcript($transcript);
    } else {
      my $seq = $transcript->translate()->seq;

      my $biotype = $transcript->biotype;

      if ($biotype eq 'pseudogene_with_CDS') {
        # Truncate the sequence if it has stop codons
        # Only keep whatever is before the first stop codon
        # (it might be an empty string, in which case the translation is empty)
        if ($seq and $seq =~ /\*/) {
          $self->log_warning("WARNING: Pseudogene_with_CDS has stop codons: truncating (" . $transcript->stable_id . ")");
          my ($seq) = split(/\*/, $seq);
        }

        # Change to normal pseudogene if there is no sequence
        if (not $seq or $seq  eq '') {
          $self->log_warning("WARNING: Changing biotype from pseudogene_with_CDS to pseudogene for " . $transcript->stable_id);
          $transcript->biotype("pseudogene");
          $transcript->translation(undef);
        }
      } elsif (not $seq or $seq eq '' or $seq =~ /\*/) {
        $self->log_warning("WARNING: No translation seq or one with the stop for transcript " . $transcript->stable_id .
                        ": $seq from " . $transcript->translateable_seq());
        $self->set_nontranslating_transcript($transcript);
      }
    }
  } else {
    $self->log_warning("WARNING: No translation (genomic start/end) for transcript " . $transcript->stable_id);
    $self->set_nontranslating_transcript($transcript);
  }
}


sub common_prefix {
  my ($self, $name, @names) = @_;
  return $name if (!@names);

  ($name, @names) = keys %{{ map { $_ => 1 } ($name, @names) }};
  return $name if (!@names);

  my $common = $name;
  foreach my $other (@names) {
    my $i = 0;
    for ($i = 0; $i < length($common) && $i < length($other); $i++) {
      last if (substr($common, $i, 1) ne substr($other, $i, 1));
    }
    $common = substr($common, 0, $i);
  }
  return $common;
}

sub get_stable_id {
  my ($self, $gff_object, @rest) = @_;
  my $use_name_field = $self->param('use_name_field');
  my $stable_ids     = $self->param_required('stable_ids');

  my $stable_id = $gff_object->load_id;
  my @all_ids = grep { defined $_ && $_ } map { $_->load_id } ($gff_object, @rest);

  if ($self->param('find_multifeature_commmon_name') && @all_ids) {
    my $new_stable_id = $self->common_prefix(@all_ids);
    if ($new_stable_id ne $stable_id) {
      $new_stable_id =~ s/[-\.]$//; # remove trailing hyphen or dot
      if ($new_stable_id) {
        $self->log_warning("using common prefix $new_stable_id as multifeature stable_id for $stable_id");
        # we assume uniquness among multifeature (CDS) prefices: same prefix -- for same CDS only
        $stable_id = $new_stable_id;
      }
    }
  }

  # but ignore multifeatures if $use_name_field
  if (defined $use_name_field && $use_name_field eq 'stable_id') {
    $stable_id = $gff_object->name if $gff_object->name;
    
    # Because the provided Name won't necessarily be unique,
    # need to ensure uniqueness ourselves...
    if (exists $$stable_ids{$stable_id}) {
      $stable_id .= '_'.$$stable_ids{$stable_id}++;
      $self->log_warning("WARNING: Added suffix to make stable ID unique: $stable_id");
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

  my $tr_stable_id = $self->get_stable_id($gff_transcript);
  
  my ($translation_id, $gff_object);
  my @gff_cds = sort sort_genomic $gff_transcript->get_SeqFeatures(@cds_types);
  my $cds_count = scalar(@gff_cds);
  if ($cds_count > 0) {
    $gff_object = $gff_cds[0];
    $translation_id = $self->get_stable_id(@gff_cds);
    $self->log_warning("get_cds_id: cds_count $cds_count > 1 for transcript $tr_stable_id, "
	                          . "using first translation_id " . ($translation_id // '')) if ($cds_count > 1);
  } else {
    $self->log_warning("get_cds_id: no gff_cds for transcript $tr_stable_id");
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
      $self->log_warning("WARNING: More than one polypeptide defined for $transcript_id. ".
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
  # produce ($genomic_start, $genomic_end) to be used by translation_coordinates
  # to infer offsets within start / end exons
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
    } elsif ($transcript->strand == 1) {
      $genomic_start += $gff_cds[0]->phase if defined $gff_cds[0]->phase;
    } else { # unknow strand for trans spliced
      # (partly duplicating self::translation_coordinates)
      # pick first "[0]" / last "[-1]" CDSs parts from unadjusted (genomically sorted) list of CDS
      
      # but because we are to iterate through adjusted exon (sort_coding) we need to adjust (circularise) CDS coordinates
      my $start_cds = [ $self->exon_coords($gff_cds[0], $transcript, "CDS"), "start" ]; # start, end, strand, label
      my $end_cds = [ $self->exon_coords($gff_cds[-1], $transcript, "CDS"), "end" ];

      # find out which exon has first or last CDS part
      # and adjust phase if there's a need
      my @exons = @{ $transcript->get_all_Exons };
      for my $exon (@exons) {
        # check if there are any CDSs within the current exon
        my @filtered_cds = grep { $exon->start <= $_->[0] && $_->[1] <= $exon->end } ($start_cds, $end_cds);

        # if we have one -- translation starts (or ends)
        if (@filtered_cds) {
          # update start/end based on label (assume exon and CDS have the same strand)
          my $label = ($exon->strand == 1)? $filtered_cds[0]->[3] : $filtered_cds[-1]->[3];
           
          # alter genomic_(start|end) only for "closest" CDSs
          if ($label eq "start" && $exon->strand == 1) {
            $genomic_start += $gff_cds[0]->phase if defined $gff_cds[0]->phase;
          }
          if ($label eq "end" && $exon->strand == -1) {
            $genomic_end   -= $gff_cds[-1]->phase if defined $gff_cds[-1]->phase;
          }
        }
      }
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
  
  # fix seq_start, seq_end for circular seq_regions
  my $slice_len = $transcript->slice->length();
  my $is_circular = $transcript->slice->is_circular();
  my $transcript_id = $transcript->stable_id;

  my ($genomic_start_raw, $genomic_end_raw) = ($genomic_start, $genomic_end);
  $genomic_start = $self->circularise_coord($genomic_start, $slice_len, $is_circular, "fixing translation_coordinates genomic_start for $transcript_id");
  $genomic_end = $self->circularise_coord($genomic_end, $slice_len, $is_circular, "fixing translation_coordinates genomic_end for $transcript_id");

  foreach my $exon (@exons) {
    if ($exon->start <= $genomic_start && $genomic_start <= $exon->end) {
      # genomic_start (fixed coords) is within exon
      if ($exon->strand == -1) {
        $end_exon = $exon;
        $seq_end = $exon->end - $genomic_start + 1;
      } else {
        $start_exon = $exon;
        $seq_start = $genomic_start - $exon->start + 1;
      }
    }
    # genomic_end (fixed coords) is within exon
    if ($exon->start <= $genomic_end && $genomic_end <= $exon->end) {
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

  my $coding_region_start = $transcript->coding_region_start();
  my $coding_region_end = $transcript->coding_region_end();
  return if (!defined $coding_region_start || !defined $coding_region_end);

  # warn if Transcript::coding_region_start invariant
  #   "By convention, the coding_region_start is always lower than the value returned by the coding_end method"
  # is broken
  # https://github.com/Ensembl/ensembl/blob/ae2dd9f7392c152c1aa07fa70c7eca4416fc1171/modules/Bio/EnsEMBL/Transcript.pm#L1072
  if ($coding_region_end < $coding_region_start) {
    $self->log_warning("Set exon phase for transcript " . $transcript->stable_id . ": (coding_region_start <= coding_region_end) ($coding_region_start <= $coding_region_end) invariant is broken. Can be a transspliced one and/or on a circular region");
  }
  
  # The phase and end_phase have defaults of -1, so only need to change
  # these when dealing with coding regions. The exons are automatically
  # returned in 5' -> 3' order.
  my $fixed_start_phase = 0;
  my $translatable = 0;
  foreach my $exon (@{$exons}) {
    # find if translation starts / ends within this exon
    my $started = $exon->start <= $coding_region_start && $coding_region_start <= $exon->end;
    my $ended = $exon->start <= $coding_region_end && $coding_region_end <= $exon->end;
    
    # with weird circular stuff not sure what we see first
    $translatable ^= 1 if ($started || $ended);
    if ($translatable || $started || $ended) {
      if (!$fixed_start_phase) {
        if ($translation->start == 1) {
          $phase = 0;
          $end_phase = ($exon->length) % 3;
        } else {
          $phase = -1;
          my $offset = $translation->start - 1;
          $end_phase = ($exon->length - $offset) % 3;
        } # translation start else
        $fixed_start_phase = 1;
      } else {
        # if already fixed start phase
        $end_phase = ($phase + $exon->length) % 3;
      }
      # for any exon
      $exon->phase($phase);
      $exon->end_phase($end_phase);
      $phase = $end_phase;
    } # $translatable

    # mark everything else as not translatable
    $translatable = 0 if ($started && $ended); # same start & end exon
  } # foreach exon
  
  # End phase of -1 is conditional on there being a 3' UTR.
  if ($translation->end < $translation->end_Exon->length) {
    $translation->end_Exon->end_phase(-1);
  }
}

sub set_nontranslating_transcript {
  my ($self, $transcript) = @_;
  my $nontranslating = $self->param_required('nontranslating');
  
  $transcript->biotype($nontranslating);

  $self->log_warning("Setting transcript " . $transcript->stable_id . " as nontranslating");
  
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
    my $igv_transcript = 0;
    my $trv_transcript = 0;
    my $igc_transcript = 0;
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      if ($transcript->biotype eq 'protein_coding') {
        $protein_coding_transcript = 1;
      } elsif ($transcript->biotype eq 'nontranslating_CDS') {
        $nontranslating_transcript = 1;
      } elsif ($transcript->biotype eq 'IG_V_gene') {
        $igv_transcript = 1;
      } elsif ($transcript->biotype eq 'TR_V_gene') {
        $trv_transcript = 1;
      } elsif ($transcript->biotype eq 'IG_C_gene') {
        $igc_transcript = 1;
      }
    }

    if ($protein_coding_transcript) {
      $gene->biotype('protein_coding');
    } elsif ($igv_transcript) {
      $gene->biotype('IG_V_gene');
    } elsif ($trv_transcript) {
      $gene->biotype('TR_V_gene');
    } elsif ($igc_transcript) {
      $gene->biotype('IG_C_gene');
    } elsif ($nontranslating_transcript) {
      $self->log_warning("Set " . $gene->stable_id . " biotype to nontranslating_CDS (from set_nontranslating_gene, there are 'nontranslating_CDS' transcripts)");
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
  my $biotype_attr = $attr{'biotype'} ? $attr{biotype}->[0] : undef;
  
  # Replace biotype with the more precise biotype attribute
  if ($biotype_attr) {
    $biotype_attr = lc($biotype_attr);
    my $known_biotypes = $self->get_db_biotypes();
    my $mapped_biotype = $known_biotypes->{$biotype_attr};
    
    return $mapped_biotype if $mapped_biotype;
  }
}

sub get_db_biotypes {
  my ($self) = @_;
  
  if (not $self->{_known_biotypes}) {
    my $dba = $self->get_type_dba();
    my $ba = $dba->get_adaptor('Biotype');
    my %biotypes = map { lc($_->name) => $_->name } @{$ba->fetch_all()};

    # add biotypes from map
    for my $feature_type (keys %$biotype_map) {
      for my $bt (keys %{$biotype_map->{$feature_type}}) {
        my $bt_lc = lc($bt);
        # biotypes from db has higher priorities
        $biotypes{$bt_lc} = $biotype_map->{$feature_type}->{$bt} if (!exists $biotypes{$bt_lc});
      }
    }

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
  
  $self->log_warning("Adding new gene: stable_id $stable_id biotype $biotype");
  
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

  # NB: do not to convert biotypes at this stage, as it will change mRNA to the ncRNA

  my ($translation_id, $gff_object) = $self->get_cds_id($gff_transcript);
  my $translatable = defined $gff_object;

  $self->log_warning("Preparing new transcript [raw data]: stable_id $stable_id biotype $transcript_type translatable $translatable gene stable_id " . $gene->stable_id . " gene_type $gene_type");
  
  # Pseudogene
  if ($gene_type eq 'pseudogene' or $transcript_type =~ /^pseudogenic_/) {
    my $tr_type = $transcript_type;
    $tr_type =~ s/^pseudogenic_//;
    
    # Protein_coding pseudogene: CDS or not?
    if ($tr_type eq "transcript" or $tr_type eq "mRNA") {
      # Check if there is a translation

      if ($translatable and $self->param('load_pseudogene_with_CDS')) {
        $self->log_warning("Pseudogene has CDSs: $stable_id");
        $biotype = 'pseudogene_with_CDS';
      } else {
        $self->log_warning("Pseudogene has no CDS: $stable_id ($transcript_type, $gene_type)");
        $biotype = 'pseudogene';
      }
      
    # Pseudogenic_tRNA
    } elsif ($tr_type eq "tRNA") {
        $self->log_warning("Pseudogenic tRNA: $stable_id");
        $biotype = 'tRNA_pseudogene';
        $gene->biotype($biotype);

    # Pseudogenic_rRNA
    } elsif ($tr_type eq "rRNA") {
        $self->log_warning("Pseudogenic rRNA: $stable_id");
        $biotype = 'rRNA_pseudogene';
        $gene->biotype($biotype);
    
    # Other kinds of pseudogenes?
    } else {
      $self->log_throw("Unrecognized pseudogene biotype: $transcript_type (gene: $gene_type)");
    }
    
  # transposable_element
  } elsif ($gene_type eq 'transposable_element') {
    $biotype = "transposable_element";
    $stable_id = $gene->stable_id . '-RA' if (!$stable_id && $gene->stable_id);
    $self->log_warning("Transposable_element: $stable_id");
  
  # Protein coding: if there are CDSs, don't rely on the biotype
  } elsif ($translatable) {
    if ($transcript_type eq "IG_V_gene") {
      $biotype = "IG_V_gene";
      $self->log_warning("IG V gene: $stable_id");
    } elsif ($transcript_type eq "TR_V_gene") {
      $biotype = "TR_V_gene";
      $self->log_warning("TR V gene: $stable_id");
    } elsif ($transcript_type eq "IG_C_gene") {
      $biotype = "IG_C_gene";
      $self->log_warning("IG C gene: $stable_id");
    } else {
      $biotype = "protein_coding";
      $self->log_warning("Protein_coding: $stable_id");
    }
  } else { # Non protein coding 
    $biotype = $transcript_type;
    $biotype = $self->map_biotype_transcript($biotype, $gff_transcript);
    $self->log_warning("Non-protein_coding: $stable_id");

    # NB: we may need a control of what biotypes are known or not
    # there's a healthcheck: genes have at least one transcript with a matching biotype group
    #    see: ensembl-datacheck/lib/Bio/EnsEMBL/DataCheck/Checks/GeneBiotypes.pm
    if (exists $transfer_biotype_tr2gene->{$biotype} || !$gene_type || lc($gene_type) eq "protein_coding"){
      $self->log_warning("Setting biotype to $biotype for non-coding: $stable_id (gene " . $gene->stable_id . ")");
      $gene->biotype($biotype);
    } else {
      $self->log_warning("Not updating biotype to $biotype for non-coding: $stable_id (gene " . $gene->stable_id . " gene type " . $gene_type . ")");
    }
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
  
  $self->log_warning("Adding new transcript: stable_id $stable_id biotype $biotype gene stable_id " . $gene->stable_id . " gene biotype " . $gene->biotype);

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


sub circularise_coord {
  my ($self, $coord, $length, $circular, $msg) = @_;
  return $coord if ($coord <= $length);

  if ($coord > $length && !$circular) {
    $self->log_throw("coord behind the non-circular slice end: $msg");
  } else {
    # moving to 0-based coords, mod len and moving back to 1-based coords
    $coord = int( ($coord - 1) % $length + 1 );
    $self->log_warning("coord behind the non-circular slice end, circularising (mod $length, $coord): $msg");
  }
  return $coord;
}


sub exon_coords {
  my ($self, $gff_exon, $transcript, $exon_id) = @_;
  # gets coords from $gff_exon, wraps them for circular seq_regions
  #   throws error for non-circular seq_regions behind the seq_region length
  #   returns ($exon_start, $exon_end, $exon_strand)

  my $exon_start = $gff_exon->start;
  my $exon_end = $gff_exon->end;
  my $exon_strand = $gff_exon->strand;

  my $is_circular = $transcript->slice->is_circular(); 
  my $slice_len = $transcript->slice->length(); 
  my $slice_name = $transcript->slice->name(); 

  my $msg = "exon $exon_id (start $exon_start, end $exon_end, strand $exon_strand) on slice $slice_name";
  $self->log_throw("not a valid strand for $msg") if ($exon_strand != 1 && $exon_strand != -1);

  # circularise only if seq_region start is not within exon (exon[ | ]exon)
  if ($exon_start > $slice_len) {
    $exon_start = $self->circularise_coord($exon_start, $slice_len, $is_circular, $msg);
    $exon_end = $self->circularise_coord($exon_end, $slice_len, $is_circular, $msg);
  }

  return ($exon_start, $exon_end, $exon_strand);
}


sub new_exon {
  my ($self, $gff_exon, $transcript, $exon_id) = @_;

  my ($exon_start, $exon_end, $exon_strand) = $self->exon_coords($gff_exon, $transcript, $exon_id);
  
  my $exon = Bio::EnsEMBL::Exon->new(
    -stable_id     => $exon_id,
    -slice         => $transcript->slice,
    -start         => $exon_start,
    -end           => $exon_end,
    -strand        => $exon_strand,
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

  my ($exon_start, $exon_end, $exon_strand) = $self->exon_coords($gff_exon, $transcript, "predicted");
  
  my $exon = Bio::EnsEMBL::PredictionExon->new(
    -slice   => $transcript->slice,
    -start   => $exon_start,
    -end     => $exon_end,
    -strand  => $exon_strand,
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

sub sort_genomic_desc {
  return $b->start <=> $a->start;
}

1;
