#!/usr/bin/env perl

=pod

=head1 NAME

  load_fann.pl

=head1 SYNOPSIS

  Load JSON-stored functional annottation into coredb.

=head1 DESCRIPTION

  Add functional annotation from JSON dumps, conforming functional_annotation_schema.json

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
         -help

=head1 EXAMPLE

  perl ./load_fann.pl \
    -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> \
    -dbname <core_db> \
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
my $help = 0;

&GetOptions(
  'host=s'                     => \$host,
  'port=s'                     => \$port,
  'user=s'                     => \$user,
  'pass=s'                     => \$pass,
  'dbname=s'                   => \$dbname,
  'json=s'                     => \$filename,
  'display_db_default=s'       => \$display_db_default,
  'feature_version_default=i'  => \$feature_version_default,
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
$display_db_default //= 'VB_Community_Annotation';

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname,
);

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

while (<$fh>) {
  chomp;
  my $data = array_ref(decode_json($_));
  for my $it (@$data) {
    my $do_update = 0;
    my ($id, $type) = map {$it->{$_}} qw/ id object_type /;
    my $lc_type = lc($type);

    my $adaptor = get_adaptor($dba, $type);
    if (not defined $adaptor) {
      warn qq/can't get adaptor for "$type" (id: "$id"). skipping...\n/;
      next;
    }
    my $obj = $adaptor->fetch_by_stable_id($id);
    if (not defined $obj) {
      warn qq/can't get object for "$id" (type: "$type"). skipping...\n/;
      next;
    }

    # gene description
    if ($lc_type eq "gene") {
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
    if ($display_xref and not grep { $_->{id} eq $display_xref } @xrefs) {
      my $dxref = {
        id => $display_xref,
        dbname => $display_db_default,
        info_type => 'DIRECT',
      };
      push @xrefs, $dxref;
    }

    my $already_used = 0;
    my $stored_xref = undef;
    for my $xref (@xrefs) {
      # "attach" synonyms to the xref with the display_xref_name
      #  or to the first seen xref
      my $attach_syns = 0;
      if (defined $display_xref) {
        $attach_syns = $display_xref eq $xref->{id};
      } else {
        $attach_syns = !$already_used;
        $already_used = 1;
      }

      # remove 'self-synonyms'
      my $add_list = $attach_syns
        ? [ grep {$_ ne $xref->{id} } @{$syns || []} ]
        : undef;

      my $xref_db_entry = store_xref(
        $dbea,
        $lc_type,
        $obj->dbID,
        $xref->{dbname},
        $xref->{id},
        $xref->{id},
        $add_list,
        $xref->{description},
        $xref->{info_type},
        $xref->{info_text}
      );

      # update 'display_xref' only for the first time or for the $set_display_xref_4
      if (defined $display_xref && $display_xref eq $xref->{id}) {
        $obj->display_xref($xref_db_entry);
        $do_update = 1;
        $stored_xref = $xref->{dbname}.':'.$xref->{id};
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
}

$dba->dbc && $dba->dbc->disconnect_if_idle();
close($fh);


# utils
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
