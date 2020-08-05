package Bio::EnsEMBL::Pipeline::Runnable::BRC4::GetMetaValue;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub run {
  my ($self) = @_;
  
  my $key = $self->param_required("param_key");
  my $name = $self->param_required("param_name");

  my $dba = $self->core_dba();
  my $ma = $dba->get_adaptor('MetaContainer');
  
  my $value = $ma->single_value_by_key($key);
  
  print("Get value for $key: $value\n");
  
  $self->dataflow_output_id({ $name => $value }, 2);
}

1;
