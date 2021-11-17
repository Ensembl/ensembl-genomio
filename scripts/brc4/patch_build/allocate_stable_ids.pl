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
use Bio::EnsEMBL::DBEntry;
use Try::Tiny;
use File::Path qw(make_path);
use LWP::UserAgent;
use HTTP::Request;
use JSON;

{
  package OSID_service;
  
  use Moose;
  
  has 'url' => (
    is => 'rw',
    isa => 'Str',
  );
  has 'user' => (
    is => 'rw',
    isa => 'Str',
  );
  has 'pass' => (
    is => 'rw',
    isa => 'Str',
  );
  has 'species' => (
    is => 'rw',
    isa => 'Str',
  );
  has 'species_id' => (
    is => 'rw',
    isa => 'Int',
  );

  sub request {
    my ($self, $type, $page, $data, $value) = @_;
    
    my $json = JSON->new->allow_nonref;
    my $ua = LWP::UserAgent->new();
    my $url = $self->url() . '/' . $page;
    $url .= '/' . $value if $value;
    
    my $request;
    if ($type eq 'get') {
      $url .= "?" . join("&", map { "$_=$data->{$_}" } (sort keys %$data));
      $request = HTTP::Request->new(GET => $url);
    } elsif ($type eq 'post') {
      $request = HTTP::Request->new(POST => $url);
      $request->header('content-type' => 'application/json');
      $request->content($json->encode($data));
    } elsif ($type eq 'patch') {
      $request = HTTP::Request->new(PATCH => $url);
      $request->header('content-type' => 'application/json');
      $request->content($json->encode($data));
    }
    $request->authorization_basic($self->user(), $self->pass());
    my $res = $ua->request($request);
    
    if ($res->is_success) {
      return $json->decode($res->content());
    } else {
      die $res->code . " : '". $res->message() . "'";
    }
  }
  
  sub request_get {
    my ($self, $page, $data) = @_;
    return $self->request('get', $page, $data);
  }
  sub request_post {
    my ($self, $page, $data) = @_;
    return $self->request('post', $page, $data);
  }
  sub request_patch {
    my ($self, $page, $data, $value) = @_;
    return $self->request('patch', $page, $data, $value);
  }
  
  # Given a species name, load its id from the OSID service
  sub connect {
    my ($self, $species) = @_;
    
    my $data = $self->request_get("organisms", { organismName => $species });
    
    if (@$data == 1) {
      my $species_id = $data->[0]->{organismId};
      $logger->debug("Species $species has id $species_id");
      $self->species_id($species_id);
    } else {
      die "No data found for species $species. Make sure it is in the OSID server.";
    }
  }
  
  sub get_gene_ids {
    my ($self, $number) = @_;
    
    $logger->info("Requesting $number gene ids...");

    my $data = $self->request_post("idSets", {
        organismId => $self->species_id,
        generateGenes => $number
      });

    if ($data) {
      my $set_id = $data->{idSetId};
      my $ids_list = $data->{generatedIds};
      die "No gene ids generated (requested $number)" if @$ids_list == 0;
      my @ids = map { $_->{geneId} } @$ids_list;
      return ($set_id, \@ids);
    } else {
      die "ERROR: no data, when I expected generated ids and id set";
    }
  }

  sub get_transcripts_ids {
    my ($self, $set_id, $map) = @_;
    
    $logger->info("Requesting transcripts ids...");

    my $data = $self->request_patch("idSets", $map, $set_id);

    if ($data) {
      my $ids_list = $data->{generatedIds};
      die "No transcripts ids generated" if @$ids_list == 0;
      
      my %tr_ids;
      for my $gene (@$ids_list) {
        my $gene_id = $gene->{geneId};
        my $trs = $gene->{transcripts};
        my $prots = $gene->{proteins};
        
        $tr_ids{$gene_id} = {transcripts => $trs, proteins => $prots};
      }
      return \%tr_ids;
      
    } else {
      die "ERROR: no data, when I expected generated ids and id set";
    }
  }
}

{
  package OSID_service_dev;
  use Moose;
  extends 'OSID_service';
  
  sub connect {
    my ($self) = @_;
    
    say "Connection to FAKE OSID service";
  }

  sub get_gene_ids {
    my ($self, $count) = @_;
    
    my $prefix = 'FAKEID_';
    my @ids = ();
    for my $i (1..($count+1)) {
      push @ids, $prefix . $i;
    }
    
    return 1, \@ids;
  }

  sub get_transcripts_ids {
    my ($self, $set_id, $map) = @_;
    
    my %tr_ids;
    my @tr_ids = ();
    my @prot_ids = ();
    for my $gene_data(@$map) {
      my $gene_id = $gene_data->{geneId};
      my $count =  $gene_data->{transcripts};
      for my $i (1..($count+1)) {
        push @tr_ids, $gene_id . '.R' . $i;
        push @prot_ids, $gene_id . '.P' . $i;
      }
      $tr_ids{$gene_id} = {transcripts => \@tr_ids, proteins => \@prot_ids};
    }
    return \%tr_ids;
  }
}

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $osid;
if ($opt{update}) {
  $osid = OSID_service->new(
    url  => $opt{osid_url},
    user => $opt{osid_user},
    pass => $opt{osid_pass},
  );
} else {
  $osid = OSID_service_dev->new();
}
$osid->connect($opt{organism});
allocate_genes($osid, $registry, $opt{species}, $opt{update}, $opt{xref_source}, $opt{prefix});

###############################################################################

sub allocate_genes {
  my ($osid, $registry, $species, $update, $xref_source, $prefix) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  my $tra = $registry->get_adaptor($species, "core", "transcript");
  my $prota = $registry->get_adaptor($species, "core", "translation");
  my $dbenta = $registry->get_adaptor($species, "core", "dbentry");
  
  my $genes_count = 0;
  my $transcripts_count = 0;
  my $translations_count = 0;
  
  my @genes = @{ $ga->fetch_all() };
  
  if ($prefix) {
    my $nold = scalar(@genes);
    @genes = grep { $_->stable_id =~ /^$prefix/ } @genes;
    my $nnew = scalar(@genes);
    $logger->info("Reduce list using prefix $prefix: from $nold genes to $nnew");
  }
  $logger->info(scalar(@genes) . " genes that will have a new stable_id allocated");
  
  # How to get all the ids

  # Part 1: get gene ids, store in a map for later
  my ($set_id, $gene_ids) = $osid->get_gene_ids(scalar(@genes));
  
  my %gene_map;
  my @tran_count;
  for my $gene (@genes) {
    # Allocate gene id
    my $gene_id = shift @$gene_ids;
    $gene_map{$gene->stable_id} = $gene_id;
    
    # Count transcript for that gene, for later transcripts allocation
    my @transcripts = @{ $gene->get_all_Transcripts };

    # Make content for OSID patch
    push @tran_count, {
      'geneId' => $gene_id,
      'transcripts' => scalar(@transcripts)
    };
  }
  
  # Part 2: get transcripts ids
  my $tr_ids = $osid->get_transcripts_ids($set_id, \@tran_count);
  
  # Part 3: Actually replace the ids
  for my $gene (@genes) {
    $genes_count++;

    my $gene_id = $gene_map{$gene->stable_id};
    
    my @transcripts = @{ $gene->get_all_Transcripts };

    my $new_tran_ids = $tr_ids->{$gene_id}->{transcripts};
    my $new_prot_ids = $tr_ids->{$gene_id}->{proteins};
    
    $logger->debug("Use gene id $gene_id to replace " . $gene->stable_id);

    # Store changes
    if ($update) {
      my $old_gene_id = $gene->stable_id;
      $gene->stable_id($gene_id);
      $ga->update($gene);
      
      # Xref
      if ($xref_source) {
        my $dbentry = Bio::EnsEMBL::DBEntry->new(
          -adaptor => $dbenta,
          -dbname => $xref_source,
          -primary_id => $old_gene_id
        );
        $dbenta->store($dbentry, $gene->dbID, 'Gene');
      }
    }
    
    for my $tr (@transcripts) {
      $transcripts_count++;

      my $tran_id = shift @$new_tran_ids;
      my $prot_id = shift @$new_prot_ids;

      $logger->debug("$gene_id: Use tran id $tran_id to replace " . $tr->stable_id);

      # Store changes
      if ($update) {
        my $old_tr_id = $tr->stable_id;
        $tr->stable_id($tran_id);
        $tra->update($tr) if $update;

        # Xref
        if ($xref_source) {
          my $dbentry = Bio::EnsEMBL::DBEntry->new(
            -adaptor => $dbenta,
            -dbname => $xref_source,
            -primary_id => $old_tr_id
          );
          $dbenta->store($dbentry, $tr->dbID, 'Transcript');
        }
      }

      my $prot = $tr->translation();
      if ($prot) {
        $translations_count++;

        $logger->debug("$gene_id: Use prot id $prot_id to replace " . $prot->stable_id);

        # Store changes
        if ($update) {
          my $old_prot_id = $prot->stable_id;
          
          # There is no update for translations in the Ensembl API, so do it via SQL
          my $sth = $prota->dbc->prepare( sprintf("UPDATE translation SET stable_id ='%s' WHERE translation_id=%s", $prot_id, $prot->dbID) );
          $sth->execute();
          
          # Xref
          if ($xref_source) {
            my $dbentry = Bio::EnsEMBL::DBEntry->new(
              -adaptor => $dbenta,
              -dbname => $xref_source,
              -primary_id => $old_prot_id
            );
            $dbenta->store($dbentry, $prot->dbID, 'Translation');
          }
        }
      }
    }
  }
  
  $logger->info("$genes_count allocated ids for genes");
  $logger->info("$transcripts_count allocated ids for transcripts");
  $logger->info("$translations_count allocated ids for translations");
  $logger->info("(Use --update to make the changes to the database)") if not $update;
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
    Allocate stable ids for all genes in a given core db, from an OSID system.

    --registry <path> : Ensembl registry
    --species <str>   : production_name of one species
    --osid_url <str>
    --osid_user <str>
    --osid_pass <str> : OSID connection details

    --organism <str>  : species name in OSID, if it is not the production_name
    
    --update          : Do the actual changes (default is no OSID call and no db changes)
    --prefix <str>    : Only replace ids for genes with this prefix [optional]
    --xref_source <str>: Keep the old ids under this xref name, e.g. "RefSeq" [optional]
    
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
    "organism=s",
    "osid_url=s",
    "osid_user=s",
    "osid_pass=s",
    "update",
    "prefix:s",
    "xref_source:s",
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  $opt{organism} //= $opt{species};
  usage("OSID details needed") if not $opt{osid_url} and not $opt{osid_user} and not $opt{osid_pass};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

