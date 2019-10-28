package Bio::EnsEMBL::Pipeline::Runnable::BRC4::Manifest;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');
use JSON;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use File::Copy qw(cp);
use File::Basename qw(basename);
use Digest::MD5 qw(md5_hex);
use v5.14;

sub run {
  my ($self) = @_;

  my $manifest = $self->param_required('manifest');
  my $output_dir = $self->param_required('output_dir');
  my $species = $self->param_required('species');

  my $dir = catdir($output_dir, $species);
  make_path($dir);

  # First move all files to the species dir
  my %final_manifest;
  for my $name (keys %$manifest) {

    if (ref($manifest->{$name}) eq 'HASH') {
      foreach my $subname (sort keys %{ $manifest->{$name} }) {
        cp($manifest->{$name}->{$subname}, $dir);
        my $file = basename($manifest->{$name}->{$subname});

        open my $fh, '<', catfile($dir, $file);
        my $md5sum = Digest::MD5->new->addfile($fh)->hexdigest;
        close $fh;

        $final_manifest{$name}{$subname} = {
          file => $file,
          md5sum => $md5sum,
        };
      }
    } else {

      cp($manifest->{$name}, $dir);
      my $file = basename($manifest->{$name});

      open my $fh, '<', catfile($dir, $file) or die("$dir/$file: $!");
      my $md5sum = Digest::MD5->new->addfile($fh)->hexdigest;
      close $fh;

      $final_manifest{$name} = {
        file => $file,
        md5sum => $md5sum,
      };
    }
  }

  # Then create a manifest file
  my $manifest_path = catfile($dir, 'manifest.json');
  
  open my $json, '>', $manifest_path or
    die "Could not open $manifest_path for writing";
  print $json encode_json(\%final_manifest);
  close $json;

  $self->dataflow_output_id({ "manifest" => $manifest_path }, 2);

  return;
}

1;
