
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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::GetMetaValue;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub run {
    my ($self) = @_;

    my $key      = $self->param_required("param_key");
    my $name     = $self->param_required("param_name");
    my $provider = $self->param_required("provider");
    my $dba      = $self->core_dba();
    my $ma       = $dba->get_adaptor('MetaContainer');

    my $prov  = $ma->single_value_by_key($provider);
    my $value = $ma->single_value_by_key($key);
    print("Get value for $key: $value\n");

    if ( $prov eq "RefSeq" ) {
        $value =~ s/GCA/GCF/;
    }

    $self->dataflow_output_id( { $name => $value }, 2 );
}

1;
