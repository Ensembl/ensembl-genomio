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


package Bio::EnsEMBL::Pipeline::Runnable::BRC4::FilterGFF3;

use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec::Functions qw(catfile);

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub run {
  my ($self) = @_;

  my $out_file = $self->param_required('out_file');
  my $species = $self->param_required("species");

  my $dir = dirname($out_file);
  my $final_gff_file = catfile($dir, $species . ".gff3.gz");

  # Only continue with the gff3 that contains the whole gene set
  if ($out_file =~ /\.chr_patch_hapl_scaff\.gff3\.gz/
      or $out_file =~ /\.chr\.gff3\.gz/) {
    unlink $out_file;
  } else {
    $self->dataflow_output_id({
        "filtered_gff_file" => $out_file,
        "final_gff_file" => $final_gff_file,
      }, 2);
  }
  return;
}

1;
