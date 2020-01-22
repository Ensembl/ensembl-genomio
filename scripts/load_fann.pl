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
my ($filename);
my $help = 0;

&GetOptions(
  'host=s'      => \$host,
  'port=s'      => \$port,
  'user=s'      => \$user,
  'pass=s'      => \$pass,
  'dbname=s'    => \$dbname,
  'json=s'    => \$filename,
  'help|?'      => \$help,
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

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname,
);

my $dbea = $dba->get_DBEntryAdaptor;
my $_a = {};
sub a {
  my ($dba, $name) = @_;  
  if (!exists $_a->{$name}) {
    # warn "getting adaptor for \"$name\"\n";
    $_a->{$name} = $dba->get_adaptor($name);
  }
  return $_a->{$name};
}


while (<$fh>) {
  chomp;
  my $data = array_ref(decode_json($_));
  for my $it (@$data) {
    my ($id, $type) = map {$it->{$_}} qw/ id object_type /;

    my $a = a($dba, $type);
    if (not defined $a) {
      warn qq/can't get adaptor for "$type" (id: "$id"). skipping...\n/;
      next;
    }
    my $obj = $a->fetch_by_stable_id($id);
    if (not defined $obj) {
      warn qq/can't get object for "$id" (type: "$type"). skipping...\n/;
      next;
    }

    if (lc($type) eq "gene") {
      $obj->description($it->{description}) if (exists $it->{description} && $it->{description} !~ m/^\s*$/); 
      # update 
      eval { $a->update($obj) };
      if ($@) {
        warn "failed to update object (id: \"$id\", type \"$type\"): $@\n";
      } 
    }
    
    # xrefs
    my $ont = array_ref($it->{ontology_terms});
    my $xrefs_raw = array_ref($it->{xrefs}); 
    my $synonyms = array_ref($it->{synonyms});

    my ($set_display_xref_4, $external_syns) = get_syns($synonyms, $id, $type);

    my @xrefs = ( @$xrefs_raw, map { { id => $_, dbname => substr($_, 0, 2) } } @$ont );
    for my $xref (@xrefs) {
      # use synonyms only first time or for the specific reference (the first one) (only once, if ids match)
      my $syns = (!defined $set_display_xref_4 || $set_display_xref_4 eq $xref->{id})? $external_syns : undef; 
      $set_display_xref_4 = undef if (defined $set_display_xref_4 && $set_display_xref_4 eq $xref->{id});

      my $xref_db_entry = store_xref(
        $dbea,
        lc("$type"), $obj->dbID,
        $xref->{dbname}, $xref->{id}, $xref->{id},
        $syns,
        $xref->{description}, $xref->{info_type}, $xref->{info_text}
      );
      
      # update 'display_xref' only for the first time or for the $set_display_xref_4
      if (defined $xref_db_entry && defined $syns && scalar(@$syns) > 0) {
        eval { $obj->display_xref($xref_db_entry) };
        if ($@) {
          warn "Failed to update display_xref_id for $type $id (", $xref->{dbname}.":".$xref->{id} ,")\n";
        } else {
          $a->update($obj);
        }
      }
      # add synonyms only for the xref specified by $set_display_xref_4 or the first one
      $external_syns = undef if (! defined $set_display_xref_4);
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

