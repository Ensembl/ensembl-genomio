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
my $mapping = load_mapping($opt{mapping});
my @old_gene_ids = values(%$mapping);

my $species = $opt{species};

# Get gene ids and transcript ids from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);

my $old_genes = get_genes_data($registry, $species, \@old_gene_ids);
$logger->info(scalar(%$old_genes) . " genes from old database");

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);

# Transfer the IDs to the new genes using the mapping
my $missed_transcripts = transfer_gene_ids($registry, $species, $mapping, $old_genes, $opt{update});

# Print out the missed transcripts
print_missed($missed_transcripts, $opt{out_transcripts});

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
    my $gene_id = $gene->stable_id;

    # Get transcripts
    my %tr_data;
    for my $tr (@{$gene->get_all_Transcripts()}) {
      my $tr_fingerprint = tr_fingerprint($tr);
      $tr_data{$tr_fingerprint} = $tr;
    }
    $genes{$gene_id} = \%tr_data;
  }
  
  return \%genes;
}

sub tr_fingerprint {
  my ($tr) = @_;

  my @exon_fingerprints;
  my $exons = $tr->get_all_Exons();
  for my $exon (@$exons) {
    my @exon_coords = (
      $exon->seq_region_name,
      $exon->seq_region_start,
      $exon->seq_region_end,
      $exon->strand,
    );
    my $exon_fingerprint = join(":", @exon_coords);
    push @exon_fingerprints, $exon;
  }
  my $fingerprint = join("__", @exon_fingerprints);

  # print("Fingerprint for " . $tr->stable_id . " : $fingerprint\n");

  return $fingerprint;
}

sub transfer_gene_ids {
  my ($registry, $species, $mapping, $old_genes, $update) = @_;

  my %stats = (
    gene_id => 0,
    total => 0,
    not_found => 0,
  );
  my $ga = $registry->get_adaptor($species, "core", 'gene');

  my @missed_transcripts;
  for my $new_id (sort keys %$mapping) {
    my $old_id = $mapping->{$new_id};
    my $gene = $ga->fetch_by_stable_id($new_id);
    if (not $gene) {
      $stats{not_found}++;
      next;
    }
    $gene->stable_id($old_id);
    $stats{gene_id}++;
    $stats{total}++;

    # Update transcripts as well?
    my @missed_tr = transfer_transcripts($gene, $old_genes);
    push @missed_transcripts, @missed_tr;
    $ga->update($gene) if $update;
  }

  for my $name (sort keys %stats) {
    my $count = $stats{$name};
    print("Transfered $name = $count\n") if $count;
  }
  if ($stats{total}) {
    print("(Add --update to actually make the transfers)\n") if not $update;
  }

  return \@missed_transcripts;
}

sub transfer_transcripts {
  my ($gene, $old_genes, $stats) = @_;

  my $old_gene = $old_genes->{$gene->stable_id};
  die("No old gene to transfer transcript from") if not $old_gene;

  my @missed_trs;
  my $transcripts = $gene->get_all_Transcripts();
  for my $transcript (@$transcripts) {
    my $cur_fingerprint = tr_fingerprint($transcript);

    my $old_tr = $old_gene->{$cur_fingerprint};
    if (not $old_tr) {
      push @missed_trs, [$gene->stable_id, $transcript->stable_id];
    }
    next if not $old_tr;

    # Same transcript fingerprint = same ID (both transcript and translation)
    $transcript->stable_id($old_tr->stable_id);
    $stats->{transcript}++;

    my $translation = $transcript->translation;
    my $old_translation = $old_tr->translation;
    next if not $translation or not $old_translation;

    $translation->stable_id($old_translation->stable_id);
    $stats->{translation}++;
  }
  return @missed_trs;
}

sub print_missed {
  my ($list, $out_file) = @_;

  open my $out_fh, ">", $out_file;
  for my $ids (@$list) {
    print $out_fh join("\t", @$ids) . "\n";
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
    Transfer the gene versions and descriptions to a patched database

    --old_registry <path> : Ensembl registry
    --new_registry <path> : Ensembl registry
    --species <str>       : production_name of one species
    --mapping <path>      : Old gene ids to new gene ids mapping file
    --out_transcripts <path: List of transcript IDs not replaced by an older ID
    
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
    "out_transcripts=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};
  usage("Mapping needed") if not $opt{mapping};
  usage("Output transcript needed") if not $opt{out_transcripts};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

