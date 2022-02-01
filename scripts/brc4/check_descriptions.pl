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

  my @species = ($opt{species}) || @{$registry->get_all_species()};
  for my $species (@species) {
    say("\tSpecies: $species");
    check_core($registry, $species);
  }
}

sub check_core {
  my ($registry, $species) = @_;

  my $counts = check_genes($registry, $species);
}

sub check_genes {
  my ($registry, $species) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  
  my %count = (
    empty_transferable => 0,
    empty_not_transferable => 0,
    gene_full => 0,
    gene_and_transcript_empty => 0,
    gene_xref => 0,
    gene_empty => 0,
    gene_total => 0,
  );

  use Data::Dumper;
  
  GENE: for my $gene (@{$ga->fetch_all()}) {
    $count{gene_total}++;
    my $g_status = check_description_status($gene->description);
    if ($g_status->{full}) {
      $count{gene_full}++;
      next GENE;
    } elsif ($g_status->{xref}) {
      $count{gene_xref}++;
    } else {
      $count{gene_empty}++;
    }
    
    my @transcripts = @{$gene->get_all_Transcripts()};

    my $t_desc = "";
    for my $transc (@transcripts) {
      my $t_status = check_description_status($transc->description);
      
      if (not $t_desc) {
        $t_desc = $transc->description;
      } elsif ($transc->description and $t_desc ne $transc->description) {
        $count{empty_not_transferable}++;
        next GENE;
      }
    }
    
    if ($t_desc) {
      $count{empty_transcript_transferable}++;
    } else {
      $count{transcript_empty}++;
    }
    
    # TODO: translation
  }
  
  for my $stat (sort keys %count) {
    say("\t\t$stat: $count{$stat}");
  }
  
  return %count;
}

sub check_description_status {
  my ($desc) = @_;

  my %status = (
    full => 0,
    empty => 0,
    xref => 0,
  );
  
  if ($desc) {
    if ($desc =~ /\[Source:/) {
      $status{xref} = 1;
    } else {
      $status{full} = 1;
    }
  } else {
    $status{empty} = 1;
  }
  
  return \%status;
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

    Optional:
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
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

