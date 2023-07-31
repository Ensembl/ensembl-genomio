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

my @current_xrefs = (
  'BRC4_Community_Annotation',
  'RefSeq_gene_name',
  'PUBMED',
  'EntrezGene',
);
my %alias_xrefs = (
  VB_Community_Annotation => 'BRC4_Community_Annotation',
);

my %ok_xrefs = map { $_ => $_ } @current_xrefs;
%ok_xrefs = (%ok_xrefs, %alias_xrefs);

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

my $species = $opt{species};

# Features to transfer versions from
my @features = qw(gene);

# Get genes and transcripts from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);
my $ma = $registry->get_adaptor($species, "core", 'MetaContainer');

my $old_metadata = get_meta($ma);
my %old_data;
for my $feat (@features) {
  $old_data{$feat} = get_feat_data($registry, $species, $feat);
  $logger->info(scalar(keys %{$old_data{$feat}}) . " ${feat}s from old database");
}

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);
$ma = $registry->get_adaptor($species, "core", 'MetaContainer');
my $new_metadata = get_meta($ma);

my %new_data;
for my $feat (@features) {
  $new_data{$feat} = get_feat_data($registry, $species, $feat);
  $logger->info(scalar(keys %{$new_data{$feat}}) . " ${feat}s from new database");
}

# Transfer descriptions
if ($opt{descriptions}) {
  $logger->info("Gene descriptions transfer:");
  update_descriptions($registry, $species, $old_data{gene}, $opt{write});
}

# Transfer xrefs
if ($opt{xrefs}) {
  $logger->info("Gene xrefs transfer:");
  update_xrefs($registry, $species, $old_data{gene}, $opt{write});
}

if ($opt{events}) {
  # Get events from features differences
  my ($old_ids, $new_ids) = diff_events($old_data{gene}, $new_data{gene});

  # Load all events from the file
  my $del_events = load_deletes($opt{deletes}, $old_ids);
  my $file_events = load_events($opt{events}, $old_ids, $new_ids);

  my %events = (%$del_events, %$file_events);

  for my $event_type (sort keys %events) {
    my @feat_events = @{$events{$event_type}};
    $logger->info("Event: $event_type = " . scalar(@feat_events));
  }

  # Compile metadata for mapping
  my $metadata = {
    old_db_name => $old_metadata->{db_name},
    new_db_name => $new_metadata->{db_name},
    old_release => $old_metadata->{release},
    new_release => $new_metadata->{release},
    old_assembly => $old_metadata->{assembly},
    new_assembly => $new_metadata->{assembly},
  };

  my $ga = $registry->get_adaptor($species, "core", 'gene');
  my $dbc = $ga->dbc;
  add_events(\%events, $old_data{gene}, $dbc, $metadata, $opt{write});
}

# Transfer versions
if ($opt{versions}) {
  $logger->info("Gene versions transfer:");
  # Update the version for the features we want to transfer and update
  for my $feat (@features) {
    update_versions($registry, $species, $feat, $old_data{$feat}, $opt{write});
  }

  # Blanket version replacement for transcripts, translations, exons
  if ($opt{write}) {
    $logger->info("Transcripts, translations and exons version update");
    my $ga = $registry->get_adaptor($species, "core", 'gene');
    my $dbc = $ga->dbc;

    # Transfer the versions from genes to transcripts
    my $transcript_query = "UPDATE transcript LEFT JOIN gene USING(gene_id) SET transcript.version = gene.version";
    $dbc->do($transcript_query);

    # Also transfer the versions from transcripts to translations
    my $translation_query = "UPDATE translation LEFT JOIN transcript USING(transcript_id) SET translation.version = transcript.version";
    $dbc->do($translation_query);

    # And set the exons version to 1
    my $exon_query = "UPDATE exon SET version = 1";
    $dbc->do($exon_query);
  } else {
    $logger->info("Transcripts, translations and exons init versions were not updated (use --write to do so)");
  }
}

###############################################################################
sub load_deletes {
  my ($deletes_file, $old_ids) = @_;

  if (not $deletes_file) {
    return {};
  }

  my %events = (
    deleted => [],
  );
  open(my $deletes_fh, "<", $deletes_file);
  while (my $line = readline($deletes_fh)) {
    chomp $line;
    my $id = $line;

    if ($old_ids->{$id}) {
      push @{$events{deleted}}, {from => [$id], to => []};
      delete $old_ids->{$id};
    } else {
      $logger->warn("Warning: '$id' in the deletes file, but not deleted in the new core");
    }
  }

  return \%events;
}

sub load_events {
  my ($events_file, $old_ids, $new_ids) = @_;

  if (not $events_file) {
    return {};
  }

  my %events = (
    change => [],
    new => [],
    split => [],
    merge => [],
  );
  my %merge_to = ();
  my %split_from = ();
  open(my $events_fh, "<", $events_file);
  while (my $line = readline($events_fh)) {
    chomp $line;
    my ($id1, $event_name, $id2) = split("\t", $line);

    # New gene
    if ($id1 and not $id2) {
      push @{$events{new}}, {from => [], to => [$id1]};
      if ($new_ids->{$id1}) {
        delete $new_ids->{$id1};
      }
    # Changed gene
    } elsif ($id1 eq $id2) {
      push @{$events{change}}, {from => [$id1], to => [$id1]};
    # Merge or split
    } elsif ($event_name eq 'merge_gene' or $event_name eq 'split_gene') {
      if (not $merge_to{$id1}) {
        $merge_to{$id1} = [];
      }
      if (not $split_from{$id2}) {
        $split_from{$id2} = [];
      }
      push @{$merge_to{$id1}}, $id2;
      push @{$split_from{$id2}}, $id1;
    }
    else {
      $logger->warn("Unsupported event '$event_name'? $line");
    }
  }

  # Ensure all event ids are unique in their groups
  for my $merge_id (keys %merge_to) {
    my @ids = uniq @{ $merge_to{$merge_id} };
    $merge_to{$merge_id} = \@ids;
  }
  for my $split_id (keys %split_from) {
    my @ids = uniq @{ $split_from{$split_id} };
    $split_from{$split_id} = \@ids;
  }

  # Check for merge and splits involving the same ids (multi event)
  my $multi_events = check_multi_events(\%split_from, \%merge_to);

  for my $event_name (keys %$multi_events) {
    my @named_events = @{ $multi_events->{$event_name} };
    for my $named_event (@named_events) {
      for my $from_id (@{ $named_event->{from} }) {
        if ($old_ids->{$from_id}) {
          delete $old_ids->{$from_id};
        }
      }
      for my $to_id (@{ $named_event->{to} }) {
        if ($new_ids->{$to_id}) {
          delete $new_ids->{$to_id};
        }
      }
    }
  }
  %events = (%events, %$multi_events);

  if (%$new_ids) {
    $logger->warn(scalar(%$new_ids) . " new ids not in the event file: " . join("; ", sort keys %$new_ids));
  }

  if (%$old_ids) {
    $logger->warn(scalar(%$old_ids) . " old ids not in the event file: " . join("; ", sort keys %$old_ids));
  }

  close($events_fh);

  return \%events;
}

sub check_multi_events {
  my ( $from, $to ) = @_;

  my %events = (
    merge => [],
    split => [],
    multi => [],
  );
  for my $from_id ( sort keys %$from ) {
    next if not $from->{$from_id};
    my @to_ids = @{ $from->{$from_id} };
    my %group  = ( from => [$from_id], to => [@to_ids] );

    my $extended = 1;
    while ($extended) {
      $extended = 0;
      my %from_ids = map { $_ => 1 } @{$group{from}};
      my %to_ids = map { $_ => 1 } @{$group{to}};

      # Expand the group in the to ids
      for my $to_id (sort keys %to_ids) {
        if ( exists $to->{$to_id} ) {
          my @to_from_ids = @{ $to->{$to_id} };

          # Add to the from list?
          my @new_from_ids;
          for my $to_from_id (@to_from_ids) {
            if (not $from_ids{$to_from_id}) {
              push @new_from_ids, $to_from_id;
            }
          }
          if (@new_from_ids) {
            push @{ $group{from} }, @new_from_ids;
            $extended = 1;
          }
        }
      }
      # Expand the group in the from ids
      for my $from_id (sort keys %from_ids) {
        if ( exists $from->{$from_id} ) {
          my @from_to_ids = @{ $from->{$from_id} };

          # Add to the to list?
          my @new_to_ids;
          for my $from_to_id (@from_to_ids) {
            if (not $to_ids{$from_to_id}) {
              push @new_to_ids, $from_to_id;
            }
          }
          if (@new_to_ids) {
            push @{ $group{to} }, @new_to_ids;
            $extended = 1;
          }
        }
      }
    }

    # We have a group! Clean up and save it
    for my $from_id (@{ $group{from} }) {
      delete $from->{$from_id};
    }
    for my $to_id (@{ $group{to} }) {
      delete $to->{$to_id};
    }
    
    # Save the group, guess the event name
    if (@{ $group{from} } > 1) {
      if (@{ $group{to} } > 1) {
        push @{$events{multi}}, \%group;
      } else {
        push @{$events{merge}}, \%group;
      }
    } else {
      push @{$events{split}}, \%group;
    }
  }

  return \%events;
}

sub diff_events {
  my ($old, $new) = @_;

  my %events = (
    new => [],
    deleted => [],
  );

  my %old_ids = map { $_ => 1 } keys %$old;
  my %new_ids = map { $_ => 1 } keys %$new;

  for my $old_id (keys %old_ids) {
    if ($new_ids{$old_id}) {
      delete $old_ids{$old_id};
      delete $new_ids{$old_id};
    }
  }

  return \%old_ids, \%new_ids;
}

sub get_feat_data {
  my ($registry, $species, $feature) = @_;
  
  my %feats;
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version // 1;
    my $description;
    if ($feature eq 'gene') {
      $description = $feat->description;
    }
    my $xrefs = $feat->get_all_DBEntries();
    $feats{$id} = { version => $version, description => $description, xrefs => $xrefs };
  }
  
  return \%feats;
}

sub update_versions {
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  
  my $update_count = 0;
  my $new_count = 0;
  
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;

    if (not defined $version) {
      my $old_feat = $old_feats->{$id};
      
      if (not $old_feat) {
        $new_count++;
        next;
      }
      
      my $old_version = $old_feat->{version};
      my $new_version = 1;

      if (defined $old_version) {
        $new_version = $old_version + 1;
        $logger->debug("Updated $feature $id must be upgraded from $old_version to $new_version");
        $update_count++;
      } else {
        $logger->debug("New $feature $id must be initialized to 1");
        $new_count++;
      }

      if ($update) {
        $feat->version($new_version);
        $fa->update($feat);
      }
    }
  }
  
  $logger->info("$update_count $feature version updated");
  $logger->info("$new_count $feature version to be initialized");
  $logger->info("(Use --write to make the changes to the database))") if $update_count + $new_count > 0 and not $update;
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

sub update_xrefs {
  # Transfer the xrefs from old gene entries to new gene entries, if those do not exist
  my ($registry, $species, $old_genes, $update) = @_;
  
  my $update_count = 0;
  my $new_count = 0;
  my %no_transfer = ();
  my %yes_transfer = ();
  my $total_transfer = 0;

  my $aa = $registry->get_adaptor($species, "core", 'analysis');
  
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  my $xa = $registry->get_adaptor($species, "core", 'DBentry');
  for my $gene (@{$ga->fetch_all}) {
    my $id = $gene->stable_id;

    my $old_gene = $old_genes->{$id};
    if (not $old_gene) {
      $new_count++;
      next;
    }
    my $xrefs = $gene->get_all_DBEntries();
    my %xref_dict = map { $_->dbname => $_ } @$xrefs;

    my $old_xrefs = $old_gene->{xrefs};

    for my $xref (@$old_xrefs) {
      my $dbname = $xref->dbname;
      # We only include dbnames that we expect
      if (not exists $ok_xrefs{$dbname}) {
        $logger->debug("NO TRANSFER for gene $id xref:\t$dbname\twith ID " . $xref->primary_id);
        $no_transfer{$dbname}++;
        next;
      }
      # Rename the dbname in case we need to transfer old xrefs which had their name changed
      $dbname = $ok_xrefs{$dbname};
      $xref->dbname($dbname);
      if (not exists $xref_dict{$dbname}) {
        $yes_transfer{$dbname}++;
        $logger->debug("Transfer gene $id xref: $dbname with ID " . $xref->primary_id);
        $update_count++;
        if ($update) {
          # Ensure we use an up to date analysis
          my $analysis_name = $xref->analysis->logic_name;
          my $analysis = $aa->fetch_by_logic_name($analysis_name);
          $xref->analysis($analysis);
          $xa->store($xref, $gene->dbID, 'Gene', 1); # ignore release
          $total_transfer++;
        }
      }
    }
  }
  
  $logger->info("$update_count gene xrefs transferred");
  $logger->info("$new_count new genes, without xref to transfer");
  for my $dbname (sort keys %yes_transfer) {
    my $count = $yes_transfer{$dbname};
    $logger->info("Transfered: $count\t$dbname");
  }
  for my $dbname (sort keys %no_transfer) {
    my $count = $no_transfer{$dbname};
    $logger->info("NOT transfered: $count from external_db '$dbname'");
  }
  $logger->info("$total_transfer written xref transfers") if $total_transfer > 0;
  if ($update_count > 0) {
    if (not $update) {
      $logger->info("(Use --write to update the xrefs in the database)");
    }
  } else {
    $logger->info("(No xrefs to update in the database)");
  }
}

sub add_events {
  my ($events, $data, $dbc, $metadata, $write) = @_;

  my $feat = 'gene';
  my $mapping_id = add_mapping_session($dbc, $metadata, $write);

  for my $event_name (keys %$events) {
    my @events = @{$events->{$event_name}};
    my $nevents = scalar @events;
    my $msg = "Storing $nevents events $event_name...";
    $msg = "Fake $msg" unless $write;
    $logger->info($msg);
    for my $event (@events) {
      my @from = @{$event->{from}};
      my @to = @{$event->{to}};

      if ($event_name eq 'new') {
        insert_event($dbc, $write, [$mapping_id, $feat, undef, undef, $to[0], 1]);
      }
      elsif ($event_name eq 'deleted') {
        my $id = $from[0];
        my $version = $data->{$id}->{version} // 1;
        insert_event($dbc, $write, [$mapping_id, $feat, $id, $version, undef, undef]);
      }
      elsif ($event_name eq 'change') {
        my $id = $from[0];
        my $old_version = $data->{$id}->{version} // 1;
        my $new_version = $old_version + 1;
        insert_event($dbc, $write, [$mapping_id, $feat, $id, $old_version, $id, $new_version]);
      }
      elsif ($event_name eq 'split' or $event_name eq 'merge' or $event_name eq 'multi') {
        for my $merge_id (@from) {
          my $old_version = $data->{$merge_id}->{version} // 1;
          for my $split_id (@to) {
            insert_event($dbc, $write, [$mapping_id, $feat, $merge_id, $old_version, $split_id, 1]);
          }
          insert_event($dbc, $write, [$mapping_id, $feat, $merge_id, $old_version, undef, undef]);
        }
      } else {
        $logger->warn("Unsupported event: $event_name");
      }
    }
  }
}

sub insert_event {
  my ($dbc, $write, $values) = @_;
  my $sql = "INSERT INTO stable_id_event(mapping_session_id, type, old_stable_id, old_version, new_stable_id, new_version) VALUES(?,?,?,?,?,?)";
  my $sth = $dbc->prepare($sql);
  my $msg = "Insert values: " . join(", ", map { $_ // 'undef' } @$values);
  if ($write) {
    $logger->debug($msg);
    $sth->execute(@$values);
  } else {
    $logger->debug("Fake $msg");
  }
}

sub get_meta {
  my ($ma) = @_;

  # Get production name
  my $prod_name = get_meta_value($ma, 'species.production_name');

  # Get db name
  my $db_name = $ma->dbc->dbname;

  # Get db release
  my $release = 0;
  if ($db_name =~ /_core_(\d+)_\d+_\d+$/) {
    $release = $1;
  } else {
    die("Can't get release from db name: $db_name");
  }

  # Remove prefix if any
  if ($db_name =~ /^.+_(${prod_name}.+$)/) {
    $db_name = $1;
  }

  # Get db assembly
  my $assembly = get_meta_value($ma, 'genebuild.version');

  my %meta = (
    db_name => $db_name,
    assembly => $assembly,
    release => $release
  );

  return \%meta;  
}

sub get_meta_value {
  my ($ma, $key) = @_;

  my ($value) = @{ $ma->list_value_by_key($key) };

  return $value;
}

sub add_mapping_session {
  my ($dbc, $meta, $write) = @_;

  my $mapping_id = 1;
  if ($write) {
    $logger->info("Storing new mapping session...");
    my @fields = qw(old_db_name new_db_name old_release new_release old_assembly new_assembly);
    my $sql = "INSERT INTO mapping_session(".join(", ", @fields).", created) VALUES(?,?,?,?,?,?,NOW())";
    my $sth = $dbc->prepare($sql);
    my @values = map { $meta->{$_} } @fields;
    $sth->execute(@values);
    my $dbh = $dbc->db_handle();
    $mapping_id = $dbh->last_insert_id;
  } else {
    $logger->info("Fake storing new mapping session (set to 1).");
    use Data::Dumper;
    $logger->debug(Dumper($meta));
  }

  return $mapping_id;
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
    --species <str>   : production_name of one species

    You can do 2 things:
    - Update the history
    - Transfer the descriptions and version
    (you can do both at the same time)

    History (use both events and deletes at the same time to make them part of the same session):
    --events <path>   : Path to an events file, to update the history and versions
    --deletes <path>  : Path to a list of deleted genes (to use with the events file)

    NB: only run the history update once (otherwise you will have duplicates).
    
    Transfer:
    --descriptions : Transfer the gene descriptions
    --versions     : Transfer and increment the gene versions, and init the others
    --xrefs        : Transfer the xrefs associated with genes
    
    Use this to make actual changes:
    --write           : Do the actual changes (default is no changes to the database)
    
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
    "descriptions",
    "versions",
    "xrefs",
    "events=s",
    "deletes=s",
    "write",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

