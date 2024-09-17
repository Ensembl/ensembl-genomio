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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DbUrlToParams;

use strict;
use warnings;

use URI;

use base ('Bio::EnsEMBL::Hive::Process');


sub param_defaults {
  my ($self) = @_;

  return {
    db_urls => {},
  };
}

sub run {
  my ($self) = @_;

  my $db_urls = $self->param('db_urls') || {};
  return if (!$db_urls);

  my $out = {};
  for my $pfx (keys %$db_urls) {
    next if (!exists $db_urls->{$pfx});

    my $res = db_params_from_url ($pfx && "${pfx}_" || "", $db_urls->{$pfx});
    map { $out->{$_} = $res->{$_} } grep { defined $res->{$_} } keys %$res;
  }

  $self->dataflow_output_id($out, 1) if ($out && %$out);
}

sub db_params_from_url {
  my ($pfx, $url) = @_;

  my ($proto, @no_proto) = split (/:\/\//, $url);
  # faking scheme for URI to work
  my $uri = URI->new(join('://', ('ftp', @no_proto)))->canonical;

  my $uri_str = $uri->as_string;
  $uri_str =~ s,^ftp://,${proto}://,;

  my $path = $uri->path;
  if (defined $path) {
    $path =~ s,/+,/,g;
    $path =~ s,^/,,;
    $path =~ s,/$,,;
  }

  my %pre = (
    url  => $uri_str,
    host => $uri->host,
    port => $uri->port,
    user => $uri->user,
    pass => $uri->password,
    dbname => $path,
  );

  return { map { $pfx . $_ => $pre{$_} } grep { defined $pre{$_} } keys %pre };
}

1;
