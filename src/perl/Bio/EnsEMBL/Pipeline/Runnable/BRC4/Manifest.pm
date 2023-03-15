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


package Bio::EnsEMBL::Pipeline::Runnable::BRC4::Manifest;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');
use JSON;
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use File::Copy qw(cp);
use File::Basename qw(basename dirname);
use Digest::MD5 qw(md5_hex);
use v5.14;

sub run {
  my ($self) = @_;

  my $species = $self->param_required('species');
  my $manifest = $self->param('manifest');
  if (not manifest_has_files($manifest)) {
    say("No manifest data for $species");
    return;
  }
  my $output_dir = $self->param_required('output_dir');
  my $species_component = $self->get_meta_value('BRC4.component');
  my $species_abbrev = $self->get_meta_value('BRC4.organism_abbrev');
  my $species_name = $species_abbrev // $species;

  my $dir = catdir($output_dir, $species_component, $species_name);
  make_path($dir);

  # First move all files to the species dir
  my %final_manifest;
  for my $name (keys %$manifest) {

    if (ref($manifest->{$name}) eq 'HASH') {
      foreach my $subname (keys %{ $manifest->{$name} }) {
        my $file = prepare_file($manifest->{$name}->{$subname}, $dir, $species, $species_abbrev);
        $final_manifest{$name}{$subname} = $file if $file;
      }
    } else {
      my $file = prepare_file($manifest->{$name}, $dir, $species, $species_abbrev);
      $final_manifest{$name} = $file if $file;
    }
  }

  # Then create a manifest file
  my $manifest_path = catfile($dir, 'manifest.json');
  
  # Print pretty JSON
  my $json = JSON->new;
  open my $jsonfh, '>', $manifest_path or
  die "Could not open $manifest_path for writing";
  print $jsonfh $json->pretty->encode(\%final_manifest);
  close $jsonfh;

  $self->dataflow_output_id({ "manifest" => $manifest_path }, 2);

  return;
}

sub manifest_has_files {
  my ($manifest) = @_;

  my $has_files = 0;
  for my $name (keys %$manifest) {

    if (ref($manifest->{$name}) eq 'HASH') {
      foreach my $subname (keys %{ $manifest->{$name} }) {
        my $file = $manifest->{$name}->{$subname};
        $has_files++ if -s $file;
      }
    } else {
      my $file = $manifest->{$name};
      $has_files++ if -s $file;
    }
  }
  
  return $has_files;
}

sub get_meta_value {
  my ($self, $key) = @_;

  my $dba = $self->core_dba();
  my $ma = $dba->get_adaptor('MetaContainer');
  my ($value) = @{ $ma->list_value_by_key($key) };

  return $value;
}

sub prepare_file {
  my ($old_path, $dir, $species, $species_abbrev) = @_;

  return if not -s $old_path;

  my $new_path = new_file_path($old_path, $dir, $species, $species_abbrev);
  cp($old_path, $new_path);
  my $md5sum = md5sum($new_path);
  my $file_name = basename($new_path);

  return {
    file => $file_name,
    md5sum => $md5sum,
  };
}

sub new_file_path {
  my ($file_path, $dir, $from, $to) = @_;

  my $file_name = basename($file_path);

  # Rename file?
  if ($to) {
    $file_name =~ s/$from/$to/;
  }
  my $new_file_path = catfile($dir, $file_name);
  return $new_file_path;
}

sub md5sum {
  my ($file_path) = @_;

  open my $fh, '<', $file_path;
  my $md5sum = Digest::MD5->new->addfile($fh)->hexdigest;
  close $fh;

  return $md5sum;
}

1;
