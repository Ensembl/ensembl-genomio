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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AddSynonyms

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AddSynonyms;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base');

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    synonym_external_db => 'ensembl_internal_synonym',
  };
}

sub run {
  my ($self) = @_;
  my $db_type    = $self->param_required('db_type');
  my $fasta_file = $self->param_required('fasta_file');
  
  my $dba = $self->get_DBAdaptor($db_type);
  
  # Fetch slices and their synonyms into a lookup hash.
  my %slices = $self->fetch_slices($dba);
  
  # Read in the fasta file.
  my %fasta = $self->load_fasta($fasta_file);
  
  $self->slice_synonyms($dba, \%slices, \%fasta)
}

sub slice_synonyms {
  my ($self, $dba, $slices, $fasta) = @_;
  
  # To save lots of obviously wrong comparisons,
  # only compare regions of the same length.
  my %slice_length;
  foreach my $slice_name (keys %$slices) {
    my $length = $$slices{$slice_name}->length;
    push @{$slice_length{$length}}, $slice_name;
  }
  
  foreach my $fasta_id (keys %$fasta) {
    if (! exists $$slices{$fasta_id}) {
      $self->warning("Warning: Fasta sequence $fasta_id does not exist in database.");
      
      # Leading Ns are sometimes in the source data, but
      # will have been stripped out by INSDC.
      $$fasta{$fasta_id} =~ s/^N*//;
      
      my $length = length($$fasta{$fasta_id});
      my %sr_ids;
      
      foreach my $slice_name (@{$slice_length{$length}}) {
        if ($$fasta{$fasta_id} eq $$slices{$slice_name}->seq) {
          $sr_ids{$$slices{$slice_name}->get_seq_region_id}++;
        }
      }
      
      my @sr_ids = keys %sr_ids;
      if (scalar(@sr_ids) == 1) {
        my $seq_region_id = $sr_ids[0];
        $self->warning("Adding $fasta_id as a synonym where seq_region_id = $seq_region_id.");
        $self->add_synonym_to_db($dba, $seq_region_id, $fasta_id);
      } elsif (scalar(@sr_ids) > 1) {
        $self->warning("Multiple seq_region_ids for $fasta_id: ".join(", ", @sr_ids));
      } else {
        $self->warning("No synonyms for $fasta_id.");
      }
    }
  }
}

sub add_synonym_to_db {
  my ($self, $dba, $seq_region_id, $synonym) = @_;
  my $synonym_external_db = $self->param_required('synonym_external_db');
  
  my $external_db_sql =
    "SELECT external_db_id FROM external_db WHERE db_name = ?;";
  
  my $synonym_sql =
    "INSERT IGNORE INTO seq_region_synonym ".
      "(seq_region_id, synonym, external_db_id) VALUES (?, ?, ?);";
  
  my $sth = $dba->dbc->prepare($external_db_sql);
  $sth->execute($synonym_external_db);
  my @results = $sth->fetchrow_array;
  
  if (scalar(@results)) {
    my $external_db_id = $results[0];
    $sth = $dba->dbc->prepare($synonym_sql);
    $sth->execute($seq_region_id, $synonym, $external_db_id);
  } else {
    die("Unrecognised external_db '$synonym_external_db'.");
  }
}

1;
