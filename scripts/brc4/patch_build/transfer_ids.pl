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
use File::Path qw(make_path);
use List::MoreUtils qw(uniq);

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

# First, load the mapping
my $mapping = load_mapping($opt{mapping})
my @old_gene_ids = values(%mapping);

my $species = $opt{species};

# Get gene ids and transcript ids from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);

my $old_genes = get_genes_data($registry, $species, \@old_gene_ids);
$logger->info(scalar(@$old_genes) . " genes from old database");

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);

# Transfer the IDs to the new genes using the mapping
transfer_gene_ids($registry, $species, $mapping)

###############################################################################
sub load_mapping {
  my ($mapping_file) = @_;

  my %mapping;
  open my $map_fh, '<', $mapping_file;
  while (my $line = readline $map_fh) {
    chomp $line;
    my ($old_id, $new_id) = split("\t", $line);
    $mapping{$new_id} = $old_id;
  }

  return \%mapping;
}

sub get_genes_data {
  my ($registry, $species, $gene_ids) = @_;
  
  my %genes;
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  for my $gene_id (@$gene_ids) {
    my $gene = $ga->fetch_by_stable_id($gene_id);
    my $gene_id = $feat->stable_id;

    # Get transcripts
    my @trs;
    for my $tr (@{$gene->get_all_Transcripts()}) {
      my $tr_data = {
        id => $tr->stable_id,
        fingerprint => tr_fingerprint($tr),
        seq_checksum => {},
      };
      my $prot = $tr->translation;
      if ($prot) {
        $tr_data{protein_id} = $prot->stable_id;
        my $seq = $prot->seq();
        $tr_data{seq_checksum}{$seq} = 1;
      }
      push @trs, $tr_data;
    }
    my %gene_data = (
      id => $gene_id,
      transcripts => \@trs,
    );
    $genes{$gene_id} = $gene_data;
  }
  
  return \%feats;
}

sub tr_fingerprint {
  my ($tr) = @_;

  
}

sub update_descriptions {
  my ($registry, $species, $old_genes, $update) = @_;
  
  my $update_count = 0;
  my $empty_count = 0;
  my $new_count = 0;
  
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  for my $gene (@{$ga->fetch_all}) {
    my $id = $gene->stable_id;
    my $description = $gene->description;

    if (not defined $description) {
      my $old_gene = $old_genes->{$id};
      if (not $old_gene) {
        $new_count++;
        next;
      }
      my $old_description = $old_gene->{description};

      if (defined $old_description) {
        my $new_description = $old_description;
        $logger->debug("Transfer gene $id description: $new_description");
        $update_count++;

        if ($update) {
          $gene->description($new_description);
          $ga->update($gene);
        }
      } else {
        $empty_count++;
      }
    }
  }
  
  $logger->info("$update_count gene descriptions transferred");
  $logger->info("$empty_count genes without description remain");
  $logger->info("$new_count new genes, without description");
  $logger->info("(Use --write to update the descriptions in the database)") if $update_count > 0 and not $update;
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
    Transfer the gene versions and descriptions to a patched database

    --old_registry <path> : Ensembl registry
    --new_registry <path> : Ensembl registry
    --species <str>       : production_name of one species
    --mapping <path>      : Old gene ids to new gene ids mapping file
    
    --update          : Do the actual changes (default is no changes to the database)
    
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
    "old_registry=s",
    "new_registry=s",
    "species=s",
    "mapping=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};
  usage("Mapping needed") if not $opt{mapping};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

