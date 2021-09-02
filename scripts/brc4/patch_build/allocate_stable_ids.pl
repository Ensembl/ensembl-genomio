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
  
  # Given a species name, load its id from the OSID service
  sub connect {
    my ($self, $species) = @_;
    
    my $ua = LWP::UserAgent->new();
    my $url = $self->url() . '/organisms';
    $url .= "?organismName=$species";
    my $request = HTTP::Request->new(GET => $url);
    $request->authorization_basic($self->user(), $self->pass());
    my $res = $ua->request($request);
    
    if ($res->is_success) {
      my $json = JSON->new->allow_nonref;
      my $data = $json->decode($res->content());
      
      if (@$data == 1) {
        $self->species_id($data->[0]->{organismId});
      } else {
          die "No data found for species $species. Make sure it is in the OSID server.";
      }
      die;
    } else {
      die $res->content();
    }
    
    my $species_id = 0;
    
    say "Species $species has id $species_id";
    $self->species_id($species_id);
  }
  sub get_gene_ids {
  }
  sub get_transcripts_ids {}
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
    
    return \@ids;
  }

  sub get_transcripts_ids {
    my ($self, $gene_id, $count) = @_;
    
    my @ids = ();
    for my $i (1..($count+1)) {
      push @ids, $gene_id . '_t' . $i;
    }
    
    return \@ids;
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
$osid->connect($opt{species});
allocate_genes($osid, $registry, $opt{species}, $opt{update});

###############################################################################
sub osid_connection {}

sub allocate_genes {
  my ($osid, $registry, $species, $update) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  
  my $genes_count = 0;
  my $transcripts_count = 0;
  my $translations_count = 0;
  
  my @genes = @{ $ga->fetch_all() };
  
  my $gene_ids = $osid->get_gene_ids(scalar(@genes));
  
  for my $gene (@genes) {
    $genes_count++;

    my $gene_id = shift @$gene_ids;
    say "Use gene id $gene_id to replace " . $gene->stable_id;
    
    my @transcripts = @{ $gene->get_all_Transcripts };
    
    my $tr_ids = $osid->get_transcripts_ids($gene_id, scalar(@transcripts));
      for my $tr (@transcripts) {
        $transcripts_count++;

        my $tr_id = shift @$tr_ids;
        say "Use transcript id $tr_id to replace " . $tr->stable_id;
        
        my $trl = $tr->translation();
        if ($trl) {
          my $trl_id = $tr_id;
          $trl_id =~ s/_t(\d+)$/_p$1/;
          say "Use translation id $trl_id to replace " . $trl->stable_id;
          $translations_count++;
        }
    }
    
    if ($update) {
      #
    }
    last;
  }
  
  say STDERR "$genes_count allocated ids for genes";
  say STDERR "$transcripts_count allocated ids for transcripts";
  say STDERR "$translations_count allocated ids for translations";
  say STDERR "(Use --update to make the changes to the database)" if not $update;
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
    
    --update          : Do the actual changes (default is no OSID call and no db changes)
    
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
    "osid_url=s",
    "osid_user=s",
    "osid_pass=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage("OSID details needed") if not $opt{osid_url} and not $opt{osid_user} and not $opt{osid_pass};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

