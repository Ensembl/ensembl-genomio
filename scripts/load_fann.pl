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

  load_fann.pl

=head1 SYNOPSIS

  Load JSON-stored functional annotation into coredb.

=head1 DESCRIPTION

  Add functional annotation from JSON dumps, conforming src/python/ensembl/io/genomio/data/schemas/functional_annotation.json

=head1 ARGUMENTS

  perl load_fann.pl
         -host
         -port
         -user
         -pass
         -dbname
         -json
         -display_db_default
         -feature_version_default
         -external_db_map
         -analysis_name
         -skip_ensembl_xrefs
         -help

=head1 EXAMPLE

  perl ./load_fann.pl \
    -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> \
    -dbname <core_db> \
    -analysis_name <logic_name> \
    -json species_functional_annotation.json

=cut

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use JSON;

my ($host, $port, $user, $pass, $dbname);
my ($filename, $display_db_default, $feature_version_default);
my ($external_db_map, $analysis_name);
my $skip_unknown_xref_source = 0;
my $skip_ensembl_xrefs = 1;
my $help = 0;
my %ensembl_xrefs = map { $_ => 1 } (
  "GO",
  "EntrezGene",
  "protein_id",
  "KEGG_Enzyme",
  "EMBL",
  "MEROPS",
  "PDB",
  "RefSeq_dna",
  "RefSeq_peptide",
  "STRING",
  "UniParc",
  "UniProtKB_all",
  "Uniprot/SWISSPROT",
  "Uniprot/SPTREMBL"
);

&GetOptions(
  'host=s'                     => \$host,
  'port=s'                     => \$port,
  'user=s'                     => \$user,
  'pass=s'                     => \$pass,
  'dbname=s'                   => \$dbname,
  'json=s'                     => \$filename,
  'display_db_default=s'       => \$display_db_default,
  'feature_version_default=i'  => \$feature_version_default,
  'external_db_map=s'          => \$external_db_map,
  'analysis_name=s'            => \$analysis_name,
  'skip_unknown_xref_source:i' => \$skip_unknown_xref_source,
  'skip_ensembl_xrefs:i'        => \$skip_ensembl_xrefs,
  'help|?'                     => \$help,
) or pod2usage(-message => "use -help", -verbose => 1);
pod2usage(-verbose => 2) if $help;

my $fh;
if (defined $filename) {
  (-e $filename) or die("-json $filename option used, but no such file. alas...\n");
  open($fh, "<", "$filename") or die("can't open json file $filename\n");
} else {
  $fh = *STDIN;
  warn "no -json filename specified. using STDIN as input.\n";
}

# Default db name for the display_xref and synonyms
$display_db_default //= 'BRC4_Community_Annotation';
$analysis_name //= 'brc4_import';

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname,
);

my $aa       = $dba->get_adaptor('Analysis');
my $analysis = $aa->fetch_by_logic_name($analysis_name);
if (! defined $analysis) {
  die "Analysis '$analysis_name' does not exist in the database.\n";
}

my $dbea = $dba->get_DBEntryAdaptor;

# Load and cache adaptors
my $_adaptors = {};
sub get_adaptor {
  my ($dba, $name) = @_;
  if (!exists $_adaptors->{$name}) {
    # warn "getting adaptor for \"$name\"\n";
    $_adaptors->{$name} = $dba->get_adaptor($name);
  }
  return $_adaptors->{$name};
}

# Slurp the json string
my $json_string;
{
  local $/; #Enable 'slurp' mode
  $json_string = <$fh>;
  close $fh;
  die("No json data") if not $json_string;
}

# Decode the json string in an array
my $data = array_ref(decode_json($json_string));

# Get external_db_map
my $extdb_map = load_external_db_map($external_db_map);

# Import each item in the array
for my $it (@$data) {
  my $do_update = 0;
  my ($id, $type) = map {$it->{$_}} qw/ id object_type /;
  if ($type eq 'transposable_element') {
    $type = 'gene';
  }
  my $lc_type = lc($type);
  my $adaptor = get_adaptor($dba, $type);
  if (not defined $adaptor) {
    warn qq/can't get adaptor for "$type" (id: "$id"). skipping...\n/;
    next;
  }

  # Sometimes we have stable IDs that have only different capitalisation.
  # This leads to adding xrefs to the [first] object returned by the `fetch_by_stable_id`.
  # Do deal with this a more complicated approach using `fetch_all_versions_by_stable_id` was used.
  # Not forgetting to failback to `fetch_by_stable_id` if there's no `fetch_all_versions_by_stable_id` option
  my $obj_list;
  if ($adaptor->can('fetch_all_versions_by_stable_id')) {
    # getting all the objects having the similar stable_id with any capitalisation and versions
    $obj_list = $adaptor->fetch_all_versions_by_stable_id($id);
  } elsif ($adaptor->can('fetch_by_stable_id')) {
    # failing back to simple single object call (for `Translations` mostly)
    $obj_list = [ $adaptor->fetch_by_stable_id($id) ];
  } else {
    warn qq/no way to fetch object for "$id" (type: "$type"). skipping...\n/;
    next;
  }
  if (not defined $obj_list) {
    warn qq/can't get objects list for "$id" (type: "$type"). skipping...\n/;
    next;
  }
  # grepping for the first object with the proper stable_id (versions are ignored here)
  my ($obj) = grep { $_ && $id eq ($_->stable_id()  // '') } @$obj_list;
  if (not defined $obj) {
    warn qq/can't get object for "$id" (type: "$type"). skipping...\n/;
    next;
  }
  # gene and transcript description
  if ($lc_type eq "gene" or $lc_type eq 'transcript') {
    $obj->description($it->{description}) if (exists $it->{description} && $it->{description} !~ m/^\s*$/);
    $do_update = 1;
  }
  # gene and transript versions
  if ($lc_type eq "gene" or $lc_type eq "transcript") {
    my $version = exists $it->{version}
        ? $it->{version}
        : $feature_version_default;
    if (defined $version) {
      $obj->version($version);
      $do_update = 1;
      # remove if fixed in core API
      update_version($dba, $lc_type, $id, $obj, $version);
    }
  }

  # xrefs
  my $ont = array_ref($it->{ontology_terms});
  my $xrefs_raw = array_ref($it->{xrefs});
  my $synonyms = array_ref($it->{synonyms});

  my ($display_xref, $syns) = get_syns($synonyms, $id, $type);
  my @xrefs = ( @$xrefs_raw, map { { id => $_, dbname => substr($_, 0, 2) } } @$ont );

  # Add display_xref to xrefs if it is not there
  my @display_xref_list = ();
  if ($display_xref) {
    my @with_disp_x = grep { $_->{id} eq $display_xref || ($_->{display_id} && $_->{display_id} eq $display_xref) } @xrefs;
    push @display_xref_list, @with_disp_x;

    if (!@display_xref_list) {
      my $dxref = {
        id => $display_xref,
        dbname => $display_db_default,
        info_type => 'DIRECT',
      };
      push @xrefs, $dxref;
    }
  }

  # Remove duplications
  @xrefs = unique_xrefs(@xrefs);
  # prepend missing GO and SO prefixes
  prepend_GO_SO(\@xrefs);
  # Do not load xrefs generated by Ensembl
  @xrefs = remove_ensembl_xrefs(\@xrefs, $extdb_map, $lc_type) if $skip_ensembl_xrefs;
  my $already_used = 0;
  my $stored_xref = undef;
  for my $xref (@xrefs) {
    # "attach" synonyms to the xref with the display_xref_name
    #  or to the first seen xref
    my $attach_syns = 0;
    if (defined $display_xref) {
      $attach_syns = $display_xref eq $xref->{id} || ($xref->{display_id} && $xref->{display_id} eq $display_xref);
    } else {
      $attach_syns = !$already_used;
      $already_used = 1;
    }

    # remove 'self-synonyms'
    my $add_list = $attach_syns
      ? [ grep {$_ ne $xref->{id} } @{$syns || []} ]
      : undef;
      
    # Used mapped external db name if it exists
    my $dbname = db_name_for_feature($extdb_map, $lc_type, $xref, $skip_unknown_xref_source);
    next if (!$dbname);

    my $xref_db_entry = store_xref(
      $dbea,
      $lc_type,
      $obj->dbID,
      $dbname,
      $xref->{id},
      $xref->{display_id} // $xref->{id},
      $add_list,
      $xref->{description},
      $xref->{info_type},
      $xref->{info_text}
    );
    # update 'display_xref' only for the first time or for the $set_display_xref_4
    if ( defined $display_xref && ( $display_xref eq $xref->{id} || $xref->{display_id} && $xref->{display_id} eq $display_xref ) ) {
      if ($lc_type eq "gene" || $lc_type eq "transcript") {
        $obj->display_xref($xref_db_entry);
        $do_update = 1;
        $stored_xref = $xref->{dbname}.':'.$xref->{id};
      } else {
        warn "not updating display_id for $type (id: \"$id\")\n";
      }
    }
  }

  # do update
  if ($do_update) {
    eval { $adaptor->update($obj) };
    if ($@) {
      my $xref_msg = defined $stored_xref? " display_xref_id for $stored_xref ": "";
      warn "failed to update object (id: \"$id\", type \"$type\")$xref_msg: $@\n";
    }
  }
}

$dba->dbc && $dba->dbc->disconnect_if_idle();
close($fh);

# utils
sub unique_xrefs {
  my (@xrefs) = @_;

  my %done;
  my @unique_xrefs;
  for my $xref (@xrefs) {
    my $key = $xref->{id} .":". $xref->{dbname};
    push @unique_xrefs, $xref if not $done{$key}++;
  }

  return @unique_xrefs;
}

sub prepend_GO_SO {
  my ($xrefs) = @_;
  for my $x (@$xrefs) {
    my $db = $x->{dbname};
    next if !($db eq "GO" || $db eq "SO");
    $x->{id} =~ s/^(?:$db:)?/$db:/;
    $x->{display_id} =~ s/^(?:$db:)?/$db:/ if (exists $x->{display_id} && $x->{display_id});
  }
}

sub remove_ensembl_xrefs {
  # Do not load xrefs generated by Ensembl
  my ($xrefs, $extdb_map, $lc_type) = @_;
  
  my @new_xrefs;
  for my $xref (@$xrefs) {
    my $dbname = db_name_for_feature($extdb_map, $lc_type, $xref, $skip_unknown_xref_source);
    if ($dbname && !exists $ensembl_xrefs{ $dbname }) {
      push @new_xrefs, $xref;
    }
    
  }
  
  return @new_xrefs;
}

sub norm_str {
  my ($raw) = @_;
  return $raw if (!defined $raw);

  my $out = "$raw";
  $out =~ s/^\s*//;
  $out =~ s/\s*$//;
  $out = undef if (length($out) == 0);
  return $out;
}

sub add_to_db_map {
  my ($db_map, $feature, $from, $val, $pat) = @_;
  $db_map->{$feature} //= {};
  $db_map->{$feature}->{$from} //= [];
  # We are adding same dbname but different pat to the database as an array
  my %uniq_dict = map {"val:".$_->{val}."with:".($_->{pat} // "_undef_") => 1 } @{$db_map->{$feature}->{$from}};
  my $k = "val:".$val."with:".($pat // "_undef_");
  if (!exists $uniq_dict{$k} ) {
    push @{$db_map->{$feature}->{$from} }, {val => $val, pat => $pat };
  }
}

sub load_external_db_map {
  my($path) = @_;
  
  return {} if not $path;
  
  my $db_map = { VALID => { _ANY_ => {}, _OTHER_ => {} }, IGNORE => { _ANY_ => {}, _OTHER_ => {} } };

  open my $map_fh, "<", $path or die "$!: $path";
  while (my $line = readline $map_fh) {
    chomp $line;
    $line =~ s/\s*#.*//;
    next if $line =~ /^\s*$/;

    my ($from_name, $to_name, $feature, $pat) = map { norm_str($_) } split(/\t/, $line);

    $feature = uc($feature // '_ANY_');
    next if ($feature eq "SEQ_REGION");

    $from_name = uc($from_name);

    if ($to_name eq '_IGNORE_') {
      add_to_db_map($db_map->{IGNORE}, $feature, $from_name, 1, $pat);
      add_to_db_map($db_map->{IGNORE}, '_OTHER_', $from_name, 1, undef) if ($feature ne "_ANY_");
    } else {
      add_to_db_map($db_map->{VALID}, $feature, $from_name, $to_name, $pat);
      add_to_db_map($db_map->{VALID}, '_OTHER_', $from_name, 1, undef) if ($feature ne "_ANY_");
    }
}
  return $db_map;
}

sub db_name_for_feature {
  my ($db_map, $lc_type, $xref, $skip_unknown_xref_source) = @_;
  return if (!$xref || !$xref->{dbname});
  
  my $raw_dbname = $xref->{dbname};
  my $xref_id = $xref->{id};
  my $res = $skip_unknown_xref_source ?  undef : $raw_dbname;

  return $res if !$db_map;
  my $feature = uc($lc_type);
  my $from_name = uc($raw_dbname);
  my $ignore_map = $db_map->{IGNORE};
  my $ignore_feature = exists $ignore_map->{$feature} && $ignore_map->{$feature}->{$from_name} || undef;
  my $ignore_any = exists $ignore_map->{_ANY_} && $ignore_map->{_ANY_}->{$from_name} || undef;
  my $ignore_other = exists $ignore_map->{_OTHER_} && $ignore_map->{_OTHER_}->{$from_name} || undef;

  my $valid_map = $db_map->{VALID};
  my $valid_feature = exists $valid_map->{$feature} && $valid_map->{$feature}->{$from_name} || undef;
  my $valid_any = exists $valid_map->{_ANY_} && $valid_map->{_ANY_}->{$from_name} || undef;
  my $valid_other = exists $valid_map->{_OTHER_} && $valid_map->{_OTHER_}->{$from_name} || undef;

  # check if there's a specific ignore rule
  if ($ignore_feature && @$ignore_feature) {
    for my $case (@$ignore_feature) {
      my $pat = $case->{pat};
      return if (!defined $pat);
      return if ($xref_id =~ m/$pat/);
    }
  }
  # check if there's a specific valid rule
  if ($valid_feature && @$valid_feature) {
    for my $case (@$valid_feature) {
      my $pat = $case->{pat};
      return $case->{val} if (!defined $pat);
      return $case->{val} if ($xref_id =~ m/$pat/);
    } 
  }
  # check if mentioned anywhere else and no global validness; no pattern checked
  return if (($valid_other && @$valid_other) && !($valid_any && @$valid_any));
  # check global ignore
  return if ($ignore_any && @$ignore_any);
  # then check global validness
  if ($valid_any && @$valid_any) {
    return $valid_any->[0]->{val};
  }
  # return raw name or undef based on $skip_unknown_xref_source flag
  return $res;
}

sub get_syns {
  my ($synonyms, $id, $type) = @_;
  my @default_ones = map { $_->{synonym} } grep { ref($_) eq "HASH" && ($_->{default}) } @$synonyms;
  if (scalar(@default_ones) > 1) {
    warn "not a single default synonym for \"$type\" \"$id\"\n";
  }
  my @syns = map { ref($_) ne "HASH" ? $_ : $_->{synonym} } @$synonyms;
  return ($default_ones[0], \@syns)
}

sub array_ref {
  my ($o) = @_;

  return [] if (!defined($o));
  return $o if (ref($o) eq "ARRAY");
  return [ $o ];
}

sub store_xref {
  my ($dbea, $object_type, $id, $external_db_name, $external_id, $external_display, $synonyms, $description, $info_type, $info_text) = @_;
  # make an xref
  my $entry = new Bio::EnsEMBL::DBEntry(
    -adaptor     => $dbea,
    -primary_id  => $external_id,
    -display_id  => $external_display,
    -dbname      => $external_db_name,
    -description => $description,
    -info_type   => $info_type,
    -info_text   => $info_text,
    -analysis    => $analysis,
  );
  # add synonyms
  $entry->{synonyms} = $synonyms if ($synonyms and ref $synonyms eq "ARRAY");
  # store xref
  my $ignore_release = 1;
  return $entry if $dbea->store( $entry, $id, $object_type, $ignore_release);
  return;
}

sub update_version {
  my ($dba, $type, $id, $obj, $version) = @_;

  eval {
    my $sth = $dba->dbc->prepare(" UPDATE $type "
                           . " SET version = ? "
                           . " WHERE ${type}_id = ? ");
    $sth->execute($version, $obj->dbID);
    $sth->finish();
  };
  warn "failed to update object's version (id: \"$id\", type \"$type\", version \"$version\"): $@\n" if ($@);
}
