package Bio::EnsEMBL::Pipeline::Runnable::BRC4::Manifest;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');
use JSON;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use File::Copy qw(move);

sub run {
  my ($self) = @_;

  my $manifest = $self->param_required('manifest');
  my $output_dir = $self->param_required('output_dir');
  my $species = $self->param_required('species');

  my $dir = catdir($output_dir, $species);
  make_path($dir);

  # First move all files to the species dir
  for my $name (keys %$manifest) {
    move($manifest->{$name}, $dir);
  }

  # Then create a manifest file
  my $manifest_path = catfile($dir, 'manifest.json');
  
  open my $json, '>', $manifest_path or
    die "Could not open $manifest_path for writing";
  print $json encode_json($manifest);
  close $json;

  return;
}

1;
