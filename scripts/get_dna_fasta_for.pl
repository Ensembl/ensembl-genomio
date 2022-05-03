#!/usr/bin/env perl
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


=pod

=head1 NAME

  get_dna_fasta_for.pl

=head1 SYNOPSIS

  Get toplevel dna sequences for the list of ids in stdin

=head1 DESCRIPTION

  Get toplevel dna sequences for the list of ids in stdin.
  Provide list of seq_region ids as input (/dev/stdin). 

=head1 ARGUMENTS

  perl get_dna_fast_for.pl
         -dbname
         -host
         -port
         -user
         -pass
         -help

=head1 EXAMPLE

  cat spec.gff | awk -F "\t" '$0 !~ /^#/ && $3 != "seq_region" {print $1}' |
    sort | uniq |
    perl ./get_dna_fasta4.pl \
      -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> \
      -dbname <core_db>

=cut

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $port, $user, $pass, $dbname);
my $help = 0;

&GetOptions(
  'host=s'      => \$host,
  'port=s'      => \$port,
  'user=s'      => \$user,
  'pass=s'      => \$pass,
  'dbname=s'    => \$dbname,
  'help|?'      => \$help,
) or pod2usage(-message => "use -help", -verbose => 1);

pod2usage(-verbose => 2) if $help;

# get ids of interest
my %ids = ();
while (<STDIN>) {
  chomp;
  next if m/^\s*$/;
  $ids{$_} = 0;
}

my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $sa = $core_db->get_SliceAdaptor();
for my $slice (@{ $sa->fetch_all('toplevel') }) {
  my $name = $slice->seq_region_name;
  my @names = ($name, @{ $slice->get_all_synonyms() });

  if (%ids) {
    @names = grep { exists $ids{$_} } @names;
    warn (join("\t", "MULTIPLE SYNONYMS USED", @names), "\n") if (scalar(@names) > 1);
    next if (scalar(@names) < 1);
    map {$ids{$_} = 1} @names; 
  }

  print ">", join(" ", @names), "\n";
  print $slice->seq, "\n";
} 

for my $id (keys %ids) {
  warn "NO TOPLEVEL SEQ REGION\t$id\n" if (!$ids{$id});
}
