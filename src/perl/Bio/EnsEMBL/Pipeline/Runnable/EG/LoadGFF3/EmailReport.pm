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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::EmailReport

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::EmailReport;

use strict;
use warnings;

use File::Basename qw(fileparse);
use File::Path qw(make_path);

use base qw(
  Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EmailReport
  Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::ReportSeqEdits
);

sub fetch_input {
  my ($self) = @_;
  my $species            = $self->param_required('species');
  my $db_type            = $self->param_required('db_type');
  my $logic_name         = $self->param_required('logic_name');
  my $protein_fasta_file = $self->param('protein_fasta_file');

  my $biotype_report_filename = $self->param('biotype_report_filename');
  my $seq_edit_tt_report_filename = $self->param('seq_edit_tt_report_filename');
  my $seq_edit_tn_report_filename = $self->param('seq_edit_tn_report_filename');
  my $protein_seq_report_filename = $self->param('protein_seq_report_filename');
  my $protein_seq_fixes_filename = $self->param('protein_seq_fixes_filename');
  
  my $dba = $self->get_DBAdaptor($db_type);
  my $dbh = $dba->dbc->db_handle;
  
  my $report = '';
  my $biotype_report_data = $self->biotype_report($dbh, $logic_name);
  $report .= $biotype_report_data; 
  my $seq_edit_tt_report_data = $self->seq_edit_tt_report($dbh, $logic_name);
  $report .= $seq_edit_tt_report_data;  
  my $seq_edit_tn_report_data = $self->seq_edit_tn_report($dbh, $logic_name);
  $report .= $seq_edit_tn_report_data;  

  $self->dump_report($biotype_report_data, $biotype_report_filename);
  $self->dump_report($seq_edit_tt_report_data, $seq_edit_tt_report_filename);
  $self->dump_report($seq_edit_tn_report_data, $seq_edit_tn_report_filename);

  if ($protein_fasta_file && -e $protein_fasta_file) {
    my ($seq_report, $fixes) = $self->protein_seq_report($dba, $logic_name, $protein_fasta_file);
    $report .= $seq_report;
    $report .= $fixes;
    $self->dump_report($seq_report, $protein_seq_report_filename);
    $self->dump_report($fixes, $protein_seq_fixes_filename);
  }
  
  $self->param('text', $report);
}

1;
