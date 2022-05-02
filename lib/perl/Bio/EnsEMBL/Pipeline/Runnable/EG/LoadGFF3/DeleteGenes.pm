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

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::DeleteGenes

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::DeleteGenes;

use strict;
use warnings;

use base ('Bio::EnsEMBL::EGPipeline::Common::RunnableDB::Base');

sub run {
  my ($self) = @_;
  my $db_type    = $self->param_required('db_type');
  my $logic_name = $self->param_required('logic_name');
  
  my $dba = $self->get_DBAdaptor($db_type);
  my $ga  = $dba->get_adaptor('Gene');
  
  my $genes = $ga->fetch_all_by_logic_name($logic_name);
  
  foreach my $gene (@$genes) {
    $ga->remove($gene);
  }
}

1;
