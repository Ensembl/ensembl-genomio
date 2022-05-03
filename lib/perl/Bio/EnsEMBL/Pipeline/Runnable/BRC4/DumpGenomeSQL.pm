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


package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpGenomeSQL;

use strict;
use warnings;

use JSON;
use File::Spec::Functions qw(catdir catfile);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::DatabaseDumper');

sub fetch_input {
  my ($self) = @_;

  my $output_dir = $self->param_required('sql_dir');
  my $species = $self->param_required("species");

  # Get meta table
  my $dba = $self->core_dba();
  my $ma = $dba->get_adaptor('MetaContainer');
  
  # Get component
  my ($component) = @{ $ma->list_value_by_key("BRC4.component") };
  my ($organism) = @{ $ma->list_value_by_key("BRC4.organism_abbrev") };

  $dba->dbc()->disconnect_if_idle();

  if ($component) {
    $output_dir = catdir($output_dir, $component);
  }
  
  # Dump file
  my $dump_file = $organism ? $organism . ".sql.gz" : $species . "sql.gz";
  my $dump_path = catfile($output_dir, $dump_file);
  $self->param('output_file', $dump_path);
  
  my $options = [
    "--max_allowed_packet 1G",
  ];
  $self->param("dump_options", $options);
  $self->SUPER::fetch_input();
}

1;
