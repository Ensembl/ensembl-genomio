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

my %alias_xrefs = (
  VB_Community_Annotation => 'BRC4_Community_Annotation',
);
my %skip_analysis = map { $_ => 1 } (
  "xrefuniprot",
  "xrefchecksum",
  "xref_sprot_blastp",
  "xref_trembl_blastp",
  "xrefuniparc",
  "gouniprot",
  "interpro2go",
);

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

my $species = $opt{species};

# Features to transfer versions from
my @features = qw(gene transcript translation);

# Get genes and transcripts from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);
my $ma = $registry->get_adaptor($species, "core", 'MetaContainer');

#Get the metadata from the database
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
  for my $feat (@features) {
    $logger->info("Descriptions transfer for $feat:");
    update_descriptions($registry, $species, $feat, $old_data{$feat}, $opt{update});
  }
}

# Transfer xrefs
if ($opt{xrefs}) {
  for my $feat (@features) {
    $logger->info("Xrefs transfer for $feat:");
    update_xrefs($registry, $species, $feat, $old_data{$feat}, $opt{update});
  }
}

# Transfer versions
if ($opt{versions}) {
  for my $feat (@features) {
    $logger->info("Versions transfer for $feat:");
    update_versions($registry, $species, $feat, $old_data{$feat}, $opt{update});
  }
}

###############################################################################
#Retrieve all the features and store the descriptions and xrefs, the translation descriptions are not stored
sub get_feat_data {
  my ($registry, $species, $feature) = @_;
  
  my %feats;
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version // 1;
    my $description;
    if ($feature eq 'gene' or $feature eq 'transcript') {
      $description = $feat->description;
    }
    my $xrefs = $feat->get_all_DBEntries();
    $feats{$id} = { version => $version, description => $description, xrefs => $xrefs };
  }
  
  return \%feats;
}

sub update_versions {
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  return if $feature eq "translation";
  
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
  $logger->info("(Use --update to make the changes to the database))") if $update_count + $new_count > 0 and not $update;
}

#Update descriptions, excluding the feature type 'translation'
sub update_descriptions {
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  return if $feature eq "translation";
  
  my $update_count = 0;
  my $empty_count = 0;
  my $new_count = 0;
  
  my $adaptor = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$adaptor->fetch_all}) {
    my $id = $feat->stable_id;
    my $description = $feat->description;

    if (not defined $description) {
      my $old_feat = $old_feats->{$id};
      if (not $old_feat) {
        $new_count++;
        next;
      }
      my $old_description = $old_feat->{description};

      if (defined $old_description) {
        my $new_description = $old_description;
        $logger->debug("Transfer $feature $id description: $new_description");
        $update_count++;

        if ($update) {
          $feat->description($new_description);
          $adaptor->update($feat);
        }
      } else {
        $empty_count++;
      }
    }
  }
  
  $logger->info("$update_count $feature descriptions transferred");
  $logger->info("$empty_count $feature without description remain");
  $logger->info("$new_count new $feature, without description");
  $logger->info("(Use --update to update the descriptions in the database)") if $update_count > 0 and not $update;
}

#Xrefs are updated by mapping the gene stable_id
sub update_xrefs {
  # Transfer the xrefs from old entries to new entries, if those do not exist
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  
  my $update_count = 0;
  my $new_count = 0;
  my %no_transfer = ();
  my %yes_transfer = ();
  my $total_transfer = 0;

  my $aa = $registry->get_adaptor($species, "core", 'analysis');
  
  my $adaptor = $registry->get_adaptor($species, "core", $feature);
  my $xa = $registry->get_adaptor($species, "core", 'DBentry');
  for my $feat (@{$adaptor->fetch_all}) {
    my $id = $feat->stable_id;

    my $old_feat = $old_feats->{$id};
    if (not $old_feat) {
      $new_count++;
      next;
    }
    my $xrefs = $feat->get_all_DBEntries();
    my %xref_dict = map { $_->dbname => $_ } @$xrefs;

    my $old_xrefs = $old_feat->{xrefs};

    #A few xrefs in 'skip_analysis' are not transferred 
    for my $xref (@$old_xrefs) {
      my $dbname = $xref->dbname;
      $dbname = $alias_xrefs{$dbname} // $dbname;
      my $analysis = $xref->analysis;
      my $analysis_name = $analysis ? "$dbname (".$analysis->logic_name.")" : $dbname;
      if ($analysis and exists $skip_analysis{$analysis->logic_name}) {
        $logger->debug("NO TRANSFER for $feature $id xref:\t$analysis_name\twith ID " . $xref->primary_id);
        $no_transfer{$analysis_name}++;
        next;
      }

      $xref->dbname($dbname);
      if (not exists $xref_dict{$dbname}) {
        $logger->debug("Transfer $feature $id xref: $analysis_name with ID " . $xref->primary_id);
        $yes_transfer{$analysis_name}++;
        $update_count++;
        if ($update) {
          # Ensure we use an up to date analysis
          my $analysis = $xref->analysis;

          # No analysis attached? Give it a default one...
          if (not $analysis) {
            $logger->debug("No analysis for $analysis_name");
          } else {
            # Reload the analysis from this db
            my $analysis_name = $analysis->logic_name;
            $analysis = $aa->fetch_by_logic_name($analysis_name);
            $xref->analysis($analysis);
          }
          $xa->store($xref, $feat->dbID, $feature, 1); # ignore release
          $total_transfer++;
        }
      }
    }
  }
  
  $logger->info("$update_count $feature xrefs transferred");
  $logger->info("$new_count new $feature, without xref to transfer");
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
      $logger->info("(Use --update to update the xrefs in the database)");
    }
  } else {
    $logger->info("(No xrefs to update in the database)");
  }
}

#Retrieve other metadata 
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

    Transfer:
    --descriptions : Transfer the gene descriptions
    --versions     : Transfer and increment the gene versions, and init the others
    --xrefs        : Transfer the xrefs associated with genes
    
    Use this to make actual changes:
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
    "descriptions",
    "versions",
    "xrefs",
    "update",
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

