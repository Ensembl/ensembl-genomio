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

Bio::EnsEMBL::Pipeline::Runnable::EG:LoadGFF3::EmailReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ReportSeqEdits;

use strict;
use warnings;

use File::Basename qw(fileparse);
use File::Path qw(make_path);

use base qw(Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base);

sub run {
  my ($self) = @_;

  my $logic_name         = $self->param_required('logic_name');
  my $protein_fasta_file = $self->param('protein_fasta_file');

  my $biotype_report_filename = $self->param('biotype_report_filename');
  my $seq_edit_tt_report_filename = $self->param('seq_edit_tt_report_filename');
  my $seq_edit_tn_report_filename = $self->param('seq_edit_tn_report_filename');
  my $protein_seq_report_filename = $self->param('protein_seq_report_filename');
  my $protein_seq_fixes_filename = $self->param('protein_seq_fixes_filename');
  
  my $dba = $self->url2dba($self->param_required('db_url'));
  my $dbh = $dba->dbc->db_handle;
  
  my $biotype_report_data = $self->biotype_report($dbh, $logic_name);
  my $seq_edit_tt_report_data = $self->seq_edit_tt_report($dbh, $logic_name);
  my $seq_edit_tn_report_data = $self->seq_edit_tn_report($dbh, $logic_name);

  $self->dump_report($biotype_report_data, $biotype_report_filename);
  $self->dump_report($seq_edit_tt_report_data, $seq_edit_tt_report_filename);
  $self->dump_report($seq_edit_tn_report_data, $seq_edit_tn_report_filename);

  if ($protein_fasta_file && -e $protein_fasta_file) {
    my ($seq_report, $fixes) = $self->protein_seq_report($dba, $logic_name, $protein_fasta_file);
    $self->dump_report($seq_report, $protein_seq_report_filename);
    $self->dump_report($fixes, $protein_seq_fixes_filename);
  }

  $dba->dbc->disconnect_if_idle();
}

sub biotype_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      g.biotype AS gene_biotype, 
      t.biotype AS transcript_biotype, 
      COUNT(DISTINCT g.stable_id) AS count_of_genes, 
      COUNT(DISTINCT t.stable_id) AS count_of_transcripts 
    FROM 
      gene g INNER JOIN 
      transcript t USING (gene_id) INNER JOIN 
      analysis a ON g.analysis_id = a.analysis_id 
    WHERE 
      a.logic_name = ? 
    GROUP BY 
      g.biotype, 
      t.biotype 
    ORDER BY 
      g.biotype, 
      t.biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = 'Imported genes summarised by biotype:';
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub seq_edit_tt_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      t.biotype, 
      COUNT(DISTINCT t.stable_id) AS count_of_transcripts, 
      COUNT(*) AS count_of_seq_edits 
    FROM 
      analysis a INNER JOIN 
      transcript t USING (analysis_id) INNER JOIN
      transcript_attrib ta USING (transcript_id) INNER JOIN
      attrib_type at USING (attrib_type_id)
    WHERE 
      a.logic_name = ? AND
      at.code = '_rna_edit'
    GROUP BY 
      t.biotype 
    ORDER BY 
      t.biotype 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = 'Imported genes requiring transcript-level sequence edits:';
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub seq_edit_tn_report {
  my ($self, $dbh, $logic_name) = @_;
  
  my $sql = "
    SELECT 
      tt.biotype,
      at.name,
      COUNT(DISTINCT tn.stable_id) AS count_of_translations, 
      COUNT(*) AS count_of_seq_edits 
    FROM 
      analysis a INNER JOIN 
      transcript tt USING (analysis_id) INNER JOIN
      translation tn USING (transcript_id) INNER JOIN
      translation_attrib ta USING (translation_id) INNER JOIN
      attrib_type at USING (attrib_type_id)
    WHERE 
      a.logic_name = ? AND
      at.code IN ('_selenocysteine', 'amino_acid_sub')
    GROUP BY 
      tt.biotype,
      at.name 
    ORDER BY 
      tt.biotype,
      at.name 
  ;";
  
  my $sth = $dbh->prepare($sql);
  $sth->execute($logic_name);
  
  my $title = 'Imported genes requiring translation-level sequence edits:';
  my $columns = $sth->{NAME};
  my $results = $sth->fetchall_arrayref();
  
  return $self->format_table($title, $columns, $results);
}

sub protein_seq_report {
  my ($self, $dba, $logic_name, $protein_fasta_file) = @_;
  
  my $ta          = $dba->get_adaptor('Transcript');
  my $transcripts = $ta->fetch_all_by_logic_name($logic_name);  
  my %protein     = $self->load_fasta($protein_fasta_file);
  
  my @results;
  my @fixes;
  my %translations;
  
  foreach my $transcript (@$transcripts) {
    my $tt_id       = $transcript->stable_id;
    my $translation = $transcript->translation;
    
    if ($translation) {
      my $tn_id  = $translation->stable_id;
      my $db_seq = $translation->seq;
      
      if (exists $protein{$tn_id}) {
        $translations{$tn_id} = 1;
        
        my $file_seq = $protein{$tn_id};
        $file_seq =~ s/(\*|\-)$//;
        
        if ($db_seq ne $file_seq) {
          my $length = length($db_seq);

          # Initial methionine only?
          (my $met_first_db_seq = $file_seq) =~ s/^./M/;
          if ($met_first_db_seq eq $file_seq) {
            push @results, [$tt_id, $tn_id, 'First amino acid in db changed to Met', $db_seq, $file_seq];
            push @fixes, "INSERT INTO translation_attrib SELECT translation_id, 170, '1 1 M' FROM translation WHERE stable_id = '$tn_id';";
          } else {
            push @results, [$tt_id, $tn_id, 'Sequences do not match', $db_seq, $file_seq];
            push @fixes, "INSERT INTO translation_attrib SELECT translation_id, 144, '1 $length $file_seq' FROM translation WHERE stable_id = '$tn_id';";
          }
          push @fixes, "UPDATE gene INNER JOIN transcript USING (gene_id) SET gene.biotype = 'protein_coding', transcript.biotype = 'protein_coding' WHERE transcript.stable_id = '$tt_id';";
        }
        
      } else {
        push @results, [$tt_id, $tn_id, 'Not in file', $db_seq, ''];
      }
    }
  }
  
  foreach my $id (keys %protein) {
    if (! exists $translations{$id}) {
      push @results, ['', $id, 'Not in db', '', $protein{$id}];
    }
  }
  
  my $title = 'Protein sequence discrepancies:';
  my @columns = qw(
    transcript_stable_id
    translation_stable_id
    message
    db_seq
    file_seq
  );
  
  my $report = "$title\n".join("\t", @columns)."\n";
  
  foreach my $result (@results) {
    $report .= join("\t", @$result)."\n";
  }
  
  my $fixes = join("\n", @fixes)."\n";
  
  return ($report, $fixes);
}

sub dump_report {
  my ($self, $data, $path) = @_;
  return unless defined $path;

  my ($filename, $dir, undef) = fileparse($path);
  if (!-e $dir) {
    make_path($dir) or $self->throw("Failed to create directory '$dir'");
  }

  open( my $fh, ">", "$path")
    or die "Can't open > $path: $!";
  print $fh $data;
  close($fh)
}

# from eg-pipeline:  Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport
sub format_table {
  my ($self, $title, $columns, $results) = @_;

  my @lengths;
  foreach (@$columns) {
    push @lengths, length($_) + 2;
  }

  foreach (@$results) {
    for (my $i=0; $i < scalar(@$_); $i++) {
      if (defined $$_[$i]) {
        my $len = length($$_[$i]) + 2;
        $lengths[$i] = $len if $len > $lengths[$i];
      }
    }
  }

  my $table = "\n$title\n";
  $table .= '+'.join('+', map {'-' x $_ } @lengths).'+'."\n";

  for (my $i=0; $i < scalar(@lengths); $i++) {
    my $column = $$columns[$i];
    my $padding = $lengths[$i] - length($column) - 2;
    $table .= '| '.$column.(' ' x $padding).' ';
  }

  $table .= '|'."\n".'+'.join('+', map {'-' x $_ } @lengths).'+'."\n";

  foreach (@$results) {
    for (my $i=0; $i < scalar(@lengths); $i++) {
      my $value = $$_[$i] || '';
      my $padding = $lengths[$i] - length($value) - 2;
      $table .= '| '.$value.(' ' x $padding).' ';
    }
    $table .= '|'."\n"
  }

  $table .= '+'.join('+', map {'-' x $_ } @lengths).'+'."\n";

  return $table;
}

1;
