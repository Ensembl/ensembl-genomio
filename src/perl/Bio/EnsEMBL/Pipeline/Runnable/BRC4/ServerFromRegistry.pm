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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::ServerFromRegistry;

use strict;
use warnings;

use URI;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');


sub run {
  my ($self) = @_;

  my $server_url = "";
  my $db_name = "";
  my $server_params = $self->server_params_from_registry();

  use Data::Dumper;
  print(Dumper($server_params) . "\n");

  $self->dataflow_output_id($server_params, 2);
}

sub server_params_from_registry {
  my ($self) = @_;

  my $dbc = $self->core_dba()->dbc;

  my $credentials = $dbc->username;
  if ($dbc->password) {
    $credentials .= ":" . $dbc->password;
  }
  my $url = sprintf("mysql://%s@%s:%d/", $credentials, $dbc->host, $dbc->port);
  my $database = $dbc->dbname;

  my $server_params = {
    dbsrv_url => $url,
    db_name => $database
  };

  return $server_params;
}

1;
