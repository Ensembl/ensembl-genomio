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

  say("#build\torgabbrev\tspecies\tgenes\ttransferable\tputative\tnot_transferable");
  my @species = ($opt{species}) || @{$registry->get_all_species()};
  for my $species (@species) {
    my $count = check_genes($registry, $species);
    my $build = get_build($registry, $species);
    my $org = get_organism_abbrev($registry, $species);
    
    my $total = $count->{gene_total};
    my $transfer = $count->{empty_transferable};
    my $not_transfer = $count->{empty_untransferable};
    my $putative = $count->{empty_putative};
    if ($transfer) {
      say("$build\t$org\t$species\t$total\t$transfer\t$putative\t$not_transfer");
    }
  }
}

sub get_build {
  my ($registry, $species) = @_;

  my $ga = $registry->get_adaptor($species, "core", "gene");
  my $dbname = $ga->dbc->dbname;
  if ($dbname =~ /_(\d+)_\d+_\d+$/) {
    return $1;
  }
}

sub get_organism_abbrev {
  my ($registry, $species) = @_;

  my $meta = $registry->get_adaptor($species, "core", "MetaContainer");
  my ($value) = @{ $meta->list_value_by_key("BRC4.organism_abbrev") };

  return $value;
}

sub check_genes {
  my ($registry, $species) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  my $dbname = $ga->dbc->dbname;
  $logger->info("Database:\t$dbname");
  $logger->info("Species:\t$species");
  
  my %count = (
    empty_putative => 0,
    empty_untransferable => 0,
    empty_untransferable => 0,
    gene_and_transcript_empty => 0,
    gene_full => 0,
    gene_xref => 0,
    gene_empty => 0,
    gene_total => 0,
  );
  
  # To check some names are not highly repeated (e.g. hypothetical protein)
  my %tname;
  
  GENE: for my $gene (@{$ga->fetch_all()}) {
    my $stable_id = $gene->stable_id;
    
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
      if (not $t_desc) {
        $t_desc = $transc->description;
      } elsif ($transc->description and $t_desc ne $transc->description) {
        $logger->info("Gene $stable_id has several transcript descriptions:\n\t'$t_desc'\n\t'".$transc->description ."'");
        $count{empty_untransferable}++;
        next GENE;
      }
    }
    
    if ($t_desc) {
      $tname{$t_desc}++;
      my $t_status = check_description_status($t_desc);
      
      if ($t_status->{putative}) {
        $count{empty_putative}++;
      } else {
        $count{empty_transferable}++;
      }
    } else {
      $count{gene_and_transcript_empty}++;
    }
    
    # TODO: translation
  }
  
  check_repeated_names(\%tname);
  
  for my $stat (sort keys %count) {
    $logger->info(sprintf("\t%7d %s", $count{$stat}, $stat));
  }
  
  return \%count;
}

sub check_repeated_names {
  my ($name) = @_;
  
  my $sus_threshold = 20;
  my @sus = grep { $name->{$_} > $sus_threshold } keys %$name;
  for my $desc (sort { $name->{$b} <=> $name->{$a} } @sus) {
    my $status = check_description_status($desc);
    if (not $status->{putative}) {
      $logger->info("The transcript description '$desc' is repeated $name->{$desc} times");
    }
  }
}

sub check_description_status {
  my ($desc) = @_;

  my %status = (
    full => 0,
    empty => 0,
    xref => 0,
    putative => 0,
  );
  
  if ($desc) {
    if ($desc =~ /\[Source:/) {
      $status{xref} = 1;
    } elsif ($desc =~ /(conserved)? *hypothetical protein(, conserved)?/i) {
      $status{putative} = 1;
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

