#!/usr/env perl
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

use Bio::EnsEMBL::Registry;
use Try::Tiny;

use Data::Dumper;

###############################################################################
# MAIN
main();

sub main {
  # Get command line args
  my %opt = %{ opt_check() };

  $logger->info("Load registry");
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($opt{registry}, 1);
  
  check_core($registry, $opt{species});
}

sub check_core {
  my ($registry, $species) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  $logger->info("Look for genes");
  
  my %count;
  for my $gene (@{$ga->fetch_all()}) {
    my $desc = $gene->description;
    my ($is_xref, $is_empty);
    
    if ($desc) {
      $is_xref = ($desc =~ /\[Source:/);
    } else {
      $is_empty = 1;
    }
    
    $count{empty}++ if $is_empty;
    $count{xref}++ if $is_xref;
  }
  
  for my $stat (sort keys %count) {
    print("$stat: $count{$stat}\n");
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
    Check genes and transcript descriptions in a core

    --registry <path> : Ensembl registry for the core database
    --species <str>   : production_name of one species
    
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
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

