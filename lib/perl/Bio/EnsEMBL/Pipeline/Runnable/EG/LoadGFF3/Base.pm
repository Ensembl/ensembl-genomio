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

Bio::EnsEMBL::Pipeline::Runnable::EG:LoadGFF3::Base

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base;

use strict;
use warnings;

use File::Basename qw(fileparse);
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Hive::Process');
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Utils::URL qw/ parse /;


sub param_defaults {
  my ($self) = @_;
  return {
    db_type => 'core',
  };
}

sub fetch_input {
  my $self = shift @_;
  
  if (defined $self->param('escape_branch') and 
      $self->input_job->retry_count >= $self->input_job->analysis->max_retry_count) 
  {
    $self->dataflow_output_id($self->input_id, $self->param('escape_branch'));
    $self->input_job->autoflow(0);
    $self->complete_early("Failure probably due to memory limit, retrying with a higher limit.");
  }
  
  # Need explicit disconnects, so that connections are freed up
  # while the analysis is running.
  my $dba = $self->url2dba($self->param_required('db_url'));
  $dba->dbc && $dba->dbc->disconnect_if_idle();
  $self->dbc && $self->dbc->disconnect_if_idle();
}

sub fetch_slices {
  my ($self, $dba) = @_;
  
  my %slices;
  
  my $slice_adaptor = $dba->get_adaptor("Slice");
  my $synonym_adaptor = $dba->get_adaptor("SeqRegionSynonym");
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    my $seq_region_name = $slice->seq_region_name;
    $slices{$seq_region_name} = $slice;
    my $synonyms = $synonym_adaptor->get_synonyms($slice->get_seq_region_id);
    foreach my $synonym (@$synonyms) {
      $slices{$synonym->name} = $slices{$seq_region_name};
    }
  }
  
  return %slices;
}

sub load_fasta {
  my ($self, $fasta_file) = @_;
  my %fasta;
  
  open(FASTA, $fasta_file);
  
  my $id;
  while (my $row = <FASTA>) {
    if ($row =~ /^>(\S+)/) {
      $id = $1;
    } else {
      $row =~ s/\s//gm;
      $fasta{$id} .= uc($row);
    }
  }
  
  close(FASTA);
  
  return %fasta;
}

sub set_protein_coding {
  my ($self, $dba, $logic_name) = @_;
  
  my $ga = $dba->get_adaptor('Gene');
  my $ta = $dba->get_adaptor('Transcript');
  
  my $genes = $ga->fetch_all_by_logic_name($logic_name);
  
  foreach my $gene (@$genes) {
    next if $gene->biotype ne 'protein_coding' and $gene->biotype ne 'nontranslating_CDS';
    my $nontranslating_transcript = 0;
    my $protein_coding_transcript = 0;
    
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      next if $transcript->biotype ne 'protein_coding' and $gene->biotype ne 'nontranslating_CDS';
      if ($transcript->translation) {
        if ($transcript->translation->seq =~ /\*/) {
          $transcript->biotype('nontranslating_CDS');
          $nontranslating_transcript++;
        } else {
          $transcript->biotype('protein_coding');
          $protein_coding_transcript++;
        }
        $ta->update($transcript);
      }
    }
    
    if ($protein_coding_transcript) {
      $gene->biotype('protein_coding');
      $ga->update($gene);
    } elsif ($nontranslating_transcript) {
      $gene->biotype('nontranslating_CDS');
      $ga->update($gene);
    }
  }
}

sub url2dba {
  my ($self, $url) = @_;

  my $dbp = Bio::EnsEMBL::Hive::Utils::URL::parse($url);
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => $dbp->{host},
    -port   => $dbp->{port},
    -user   => $dbp->{user},
    -pass   => $dbp->{pass},
    -dbname => $dbp->{dbname},
  );
  return $dba;
}

# logging methods

sub log {
  my ($self, @msg) = @_;
  $self->{_local_log} = [] if (!defined $self->{_local_log});
  push @{$self->{_local_log}}, join(" ", @msg);
}

sub log_warning {
  my ($self, @msg) = @_;
  $self->log(@msg);
  $self->warning(@msg);
}

sub log_warn {
  my ($self, @msg) = @_;
  $self->log(@msg);
  warn(@msg);
}

sub log_throw {
  my ($self, @msg) = @_;
  $self->log(@msg, "dying...");
  $self->dump_log();
  $self->throw(@msg);
}

sub dump_log {
  my ($self) = @_;

  my $path = $self->param("log");
  return unless defined $path;

  my ($filename, $dir, undef) = fileparse($path);
  if (!-e $dir) {
    make_path($dir) or $self->throw("Failed to create directory '$dir'");
  }

  open( my $fh, ">", "$path")
    or die "Can't open > $path: $!";
  print $fh join("\n", @{ $self->{_local_log} // [] });
  close($fh)
}

1;
