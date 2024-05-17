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

use Bio::EnsEMBL::Registry;
use Try::Tiny;
use File::Path qw(make_path);
use List::MoreUtils qw(uniq);

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

# First, load the mapping of changed genes
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
my ($missed_trs, $mapped_trs) = transfer_gene_ids($registry, $species, $mapping, $old_genes, $opt{update});

# Print out the missed transcripts
print_list($missed_trs, $opt{missed_transcripts});
print_list($mapped_trs, $opt{mapped_transcripts});

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

#Get all the genes and transcripts associated through the adaptor, its a key value pair
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
      my $translation = $tr->translation;
      $tr_data{$tr_fingerprint} = {
        "stable_id" => $tr->stable_id,
        "translation_id" => $translation->stable_id,
      };
    }
    $genes{$gene_id} = \%tr_data;
  }
  
  return \%genes;
}

#Retrieve all the exons coords and pack them into "fingerprint"
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

#Gene ids are updated mapping stable_id, the old gene_id is stored as stable_id in the new db
sub transfer_gene_ids {
  my ($registry, $species, $mapping, $old_genes, $update) = @_;

  my %stats = (
    gene_id => 0,
    total => 0,
    not_found => 0,
  );
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  my $tra = $registry->get_adaptor($species, "core", 'transcript');

  my @missed_transcripts;
  my @mapped_transcripts;
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
    my ($missed_trs, $mapped_trs) = transfer_transcripts($gene, $old_genes, \%stats, $tra, $update);
    push @missed_transcripts, @$missed_trs;
    push @mapped_transcripts, @$mapped_trs;
    $ga->update($gene) if $update;
  }

  for my $name (sort keys %stats) {
    my $count = $stats{$name};
    print("Transfered $name = $count\n") if $count;
  }
  if ($stats{total}) {
    print("(Add --update to actually make the transfers)\n") if not $update;
  }

  return (\@missed_transcripts, \@mapped_transcripts);
}

# Transcripts are transferred mapping the stable_id
sub transfer_transcripts {
  my ($gene, $old_genes, $stats, $tra, $update) = @_;

  my $old_gene = $old_genes->{$gene->stable_id};
  die("No old gene to transfer transcript from") if not $old_gene;

  my @missed_trs;
  my @mapped_trs;
  my $transcripts = $gene->get_all_Transcripts();
  my @old_transcripts = values(%$old_gene);

  # Special case: 1-to-1 gene with 1 transcript = automatic transfer
  if (@$transcripts == 1 and @old_transcripts == 1) {
    $logger->debug("1-to-1 transcript transfer of ID");
    my $transcript = $transcripts->[0];
    my $old_transcript = $old_transcripts[0];
    transfer_transcript_id($old_transcript, $transcript, $stats, $tra, $update);
    push @mapped_trs, [$gene->stable_id, $old_transcript->{"stable_id"}, $old_transcript->{"translation_id"}];
  } else {
    # Otherwise: 
    for my $transcript (@$transcripts) {
      my $cur_fingerprint = tr_fingerprint($transcript);

      my $old_transcript = $old_gene->{$cur_fingerprint};
      if (not $old_transcript) {
        $logger->debug("Missed id for " . $transcript->stable_id);
        push @missed_trs, [$gene->stable_id, $transcript->stable_id];
        $stats->{missed_transcripts}++;
        next;
      }
      transfer_transcript_id($old_transcript, $transcript, $stats, $tra, $update);
      push @mapped_trs, [$gene->stable_id, $old_transcript->{"stable_id"}, $old_transcript->{"translation_id"}];
    }
  }
  return (\@missed_trs, \@mapped_trs);
}

sub transfer_transcript_id {
  my ($old_transcript, $transcript, $stats, $tra, $update) = @_;
  $logger->debug("Update id for " . $transcript->stable_id . " to " . $old_transcript->{"stable_id"});
  
  my $translation = $transcript->translation;
  my $old_translation_id = $old_transcript->{"translation_id"};

  # Stop if we don't have translation in both old and new transcripts
  if (not $old_translation_id and $translation) {
      die("Old transcript had not translation, but new one has one " . $transcript->stable_id);
  } elsif ($old_translation_id and not $translation) {
      die("Old transcript had a translation, but new one doesn't " . $transcript->stable_id);
  }

  # Actual transcript update
  $transcript->stable_id($old_transcript->{"stable_id"});
  $tra->update($transcript) if $update;
  $stats->{transcript}++;

  # Also update translations
  $logger->debug("Update id for " . $translation->stable_id . " to " . $old_translation_id);
  update_translation_id($tra, $old_translation_id, $translation->stable_id) if $update;
  $translation->stable_id($old_transcript->{"translation_id"});
  $stats->{translation}++;
}

sub update_translation_id {
  my ($tra, $old_id, $new_id) = @_;

  # There is no translationAdaptor update, so make it manually
  my $sth = $tra->prepare("UPDATE translation SET stable_id=? WHERE stable_id=?");
  $sth->bind_param(1, $old_id);
  $sth->bind_param(2, $new_id);
  $sth->execute();
}

sub print_list {
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
    --missed_transcripts <path> : List of transcript IDs not replaced by an older ID
    --mapped_transcripts <path> : List of transcript IDs replaced by an older ID
    
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
    "missed_transcripts=s",
    "mapped_transcripts=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};
  usage("Mapping needed") if not $opt{mapping};
  usage("Missed transcripts output needed") if not $opt{missed_transcripts};
  usage("Mapped transcripts output needed") if not $opt{mapped_transcripts};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

