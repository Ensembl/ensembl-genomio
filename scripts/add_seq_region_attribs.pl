#!/usr/env perl
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


use v5.14.00;
use strict;
use warnings;
use Carp;
use autodie qw(:all);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 
use File::JSON::Slurper qw(read_json);

use Bio::EnsEMBL::Registry;


my @attribs_to_update = qw(
  EBI_seq_region_name
  BRC4_seq_region_name
);
my @syns_to_update = qw(
  INSDC
);

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';

my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $species = $opt{species};

# Get all seq_regions
my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice');
my $attribute_adaptor = $registry->get_adaptor($species, 'core', 'attribute');

my $json_seq = read_json($opt{json});
my %seqs = map { $_->{name} => $_ } @$json_seq;

for my $seq_name (keys %seqs) {
  my $slice = $slice_adaptor->fetch_by_toplevel_location($seq_name);
  if ($slice) {
    my $seq = $seqs{$seq_name};

    # Update attribs
    for my $attrib_name (@attribs_to_update) {
      my $attrib_value = $seq->{$attrib_name};
      if ($attrib_value) {
        change_attrib($slice, $attrib_name, $attrib_value);
      } else {
        say "No value to update for $attrib_name in " . $slice->seq_region_name;
      }
    }

    # Update syns too
    # TODO

  } else {
    say "$seq_name\tNo slices";
  }
}

sub change_attrib {
  my ($slice, $name, $value) = @_;
  
  # Get current attrib value
  my ($cur_attr) = @{ $slice->get_all_Attributes($name) };
  if ($cur_attr) {
    my $cur_value = $cur_attr->value;
    if ($value eq $cur_value) {
      #say "Current value for attrib $name is the same ($value)";
    } else {
      say "Update the value for attrib $name from $cur_value to $value";
      
    }
  } else {
    say "Insert the new value for attrib $name ($value)";
    insert_attrib($slice, $name, $value);
  }
}

sub insert_attrib {
  my ($slice, $name, $value) = @_;

  my $attr = Bio::EnsEMBL::Attribute->new(
    -CODE => $name,
    -VALUE => $value,
  );
  my @attrs = ($attr);

  if ($opt{apply}) {
    $attribute_adaptor->store_on_Slice($slice, \@attrs);
  }
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    Update seq_region attribs and synonyms from a seq_region json from the INSDC record
    Adds:
      - BRC4_seq_region_name
      - EBI_seq_region_name
      - INSDC accession

    --registry <path>
    --species <str>
    --json <path>
    
    --apply           : Actually apply the changes to the db (default: just report)
    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "registry=s",
    "species=s",
    "json=s",
    "apply",
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage("Seq region json needed") if not $opt{json};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}
__END__

__END__

