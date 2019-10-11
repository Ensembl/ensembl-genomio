package Bio::EnsEMBL::Pipeline::Runnable::BRC4::FilterGFF3;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub run {
  my ($self) = @_;

  my $out_file = $self->param_required('out_file');

  # Only continue with the gff3 that contains the whole gene set
  if ($out_file =~ /\.chr_patch_hapl_scaff\.gff3\.gz/
      or $out_file =~ /\.chr\.gff3\.gz/) {
    unlink $out_file;
  } else {
    $self->dataflow_output_id({ "out_file" => $out_file }, 2);
  }
  return;
}


1;
