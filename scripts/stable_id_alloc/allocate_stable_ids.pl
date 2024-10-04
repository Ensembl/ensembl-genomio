#!/usr/env perl
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
use Bio::EnsEMBL::DBEntry;
use Try::Tiny;
use File::Path qw(make_path);
use LWP::UserAgent;
use Time::Piece;
use HTTP::Request;
use JSON;

my $default_analysis_name = 'brc4_import';


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
      die "No gene ids generated (requested $number)" if $number > 0 and @$ids_list == 0;
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

    return 1, [] if $count == 0;
    
    my $prefix = 'FAKEID_';
    my @ids = ();
    for my $i (1..($count+1)) {
      push @ids, $prefix . $i;
    }
    
    return 1, \@ids;
  }

  sub get_transcripts_ids {
    my ($self, $set_id, $map) = @_;
    
    my %tr_ids_out;
    for my $gene_data (@$map) {
      my @tr_ids = ();
      my @prot_ids = ();
      my $gene_id = $gene_data->{geneId};
      my $count = $gene_data->{transcripts};
      for my $i (1..($count+1)) {
        push @tr_ids, $gene_id . '.R' . $i;
        push @prot_ids, $gene_id . '.P' . $i;
      }
      $tr_ids_out{$gene_id} = {transcripts => \@tr_ids, proteins => \@prot_ids};
    }
    return \%tr_ids_out;
  }
}

###############################################################################

sub load_gene_list {
  my ($gene_file) = @_;

  return if not $gene_file;

  my @gene_ids;
  open my $file_fh, "<", $gene_file;
  while (my $line = readline $file_fh) {
    chomp($line);
    push @gene_ids, $line;
  }
  return \@gene_ids;
}

sub load_transcript_list {
  # Get a dict of gene ids and transcript ids, where the key is gene_id and the value an array of transcript_ids
  # for that gene_id
  my ($transcript_file) = @_;

  return if not $transcript_file;

  my %ids;
  open my $file_fh, "<", $transcript_file;
  while (my $line = readline $file_fh) {
    chomp($line);
    my ($gene_id, $tr_id) = split("\t", $line);
    if ($ids{$gene_id}) {
      push @{$ids{$gene_id}}, $tr_id;
    } else {
      $ids{$gene_id} = [$tr_id];
    }
  }
  return \%ids;
}

sub allocate_genes {
  my ($osid, $registry, $species, $update, $xref_source, $analysis_name, $prefix, $after_date, $output_map, $genes_list, $trs_list) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  my $tra = $registry->get_adaptor($species, "core", "transcript");
  my $prota = $registry->get_adaptor($species, "core", "translation");
  my $dbenta = $registry->get_adaptor($species, "core", "dbentry");
  my $analysa = $registry->get_adaptor($species, "core", "analysis");

  my $import_an = $analysa->fetch_by_logic_name($analysis_name);
  die("No analysis $analysis_name found") if not $import_an;
  
  my $genes_count = 0;
  my $transcripts_count = 0;
  my $translations_count = 0;
  
  my @genes;
  if ($genes_list) {
    @genes = map { $ga->fetch_by_stable_id($_) } @$genes_list;
  } elsif ($trs_list) {
    @genes = map { $ga->fetch_by_stable_id($_) } sort keys %$trs_list;
  } else {
    @genes = @{ $ga->fetch_all() };
  }
  
  if ($prefix) {
    my $nold = scalar(@genes);
    @genes = grep { $_->stable_id =~ /^$prefix/ } @genes;
    my $nnew = scalar(@genes);
    $logger->info("Reduce list using prefix $prefix: from $nold genes to $nnew");
  }
  if ($after_date) {
    my $nold = scalar(@genes);
    @genes = grep { localtime($_->created_date)->datetime gt $after_date } @genes;
    my $nnew = scalar(@genes);
    $logger->info("Reduce list using date $after_date: from $nold genes to $nnew");
  }
  $logger->info(scalar(@genes) . " genes will have a new stable_id allocated");

  # Check that all genes have been found
  my $not_defined = scalar(grep {not defined $_} @genes);
  if ($not_defined) {
    die("$not_defined genes to replace are not defined\n");
  }
  
  # How to get all the ids

  # Transcript list = don't need to get gene ids, we already have them
  my @tran_count;
  my $set_id;
  my %gene_map;
  if ($trs_list) {
    ($set_id) = $osid->get_gene_ids(0);  # Get a set id without asking for gene ids
    for my $gene (@genes) {
      my $gene_id = $gene->stable_id;
      my $trs = $trs_list->{$gene_id};
      # Make content for OSID patch
      push @tran_count, {
        'geneId' => $gene_id,
        'transcripts' => scalar @$trs
      };
    }
  } else {
    # Part 1: get gene ids, store in a map for later
    ($set_id, my $gene_ids) = $osid->get_gene_ids(scalar(@genes));
    
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
    save_gene_map(\%gene_map, $output_map);
  }
  
  # Part 2: get transcripts ids
  my $tr_ids = $osid->get_transcripts_ids($set_id, \@tran_count);
  
  # Part 3: Actually replace the ids
  for my $gene (@genes) {

    my $gene_id = $trs_list ? $gene->stable_id : $gene_map{$gene->stable_id};

    my $new_tran_ids = $tr_ids->{$gene_id}->{transcripts};
    my $new_prot_ids = $tr_ids->{$gene_id}->{proteins};
    
    if (not $trs_list) {
      $genes_count++;
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
            -primary_id => $old_gene_id,
            -display_id => $old_gene_id
          );
          
          # Add analysis ENA for this
          $dbentry->analysis($import_an);

          $dbenta->store($dbentry, $gene->dbID, 'Gene');
        }
      }
    }
    
    # Now update the transcripts
    my @transcripts = @{ $gene->get_all_Transcripts };

    # Restrict to the transcripts in the list
    if ($trs_list) {
      my %tr_ids = map { $_ => 1 } @{$trs_list->{$gene_id}};
      @transcripts = grep { $tr_ids{$_->stable_id} } @transcripts;
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
        $logger->debug("Reset ID from $old_tr_id to $tran_id: " . $tr->stable_id);

        if ($update) {
          $tra->update($tr);

          # Xref
          if ($xref_source) {
            my $tr_xref_source = $xref_source;
            if ($xref_source eq 'RefSeq') {
              if ($old_tr_id =~ /^[NX]M_/) {
                $tr_xref_source = 'RefSeq_mRNA';
              } else {
                $tr_xref_source = '';
              }
            } elsif ($xref_source eq 'GenBank') {
              $tr_xref_source = 'GenBank_transcript';
            }

            if ($tr_xref_source) {
              my $dbentry = Bio::EnsEMBL::DBEntry->new(
                -adaptor => $dbenta,
                -dbname => $tr_xref_source,
                -primary_id => $old_tr_id,
                -display_id => $old_tr_id
              );
              $dbenta->store($dbentry, $tr->dbID, 'Transcript');
            }
          }
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
            my $tl_xref_source = $xref_source;
            if ($xref_source eq 'RefSeq') {
              $tl_xref_source = 'RefSeq_peptide';
            } elsif ($xref_source eq 'GenBank') {
              $tl_xref_source = 'GenBank_translation';
            }
            my $dbentry = Bio::EnsEMBL::DBEntry->new(
              -adaptor => $dbenta,
              -dbname => $tl_xref_source,
              -primary_id => $old_prot_id,
              -display_id => $old_prot_id
            );
            my $ignore_release = 1;
            $dbenta->store($dbentry, $prot->dbID, 'Translation', $ignore_release);
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

sub save_gene_map {
  my ($map, $map_path) = @_;
  return if not $map_path;

  open my $map_fh, ">", $map_path;
  for my $old_id (sort keys %$map) {
    my $new_id = $map->{$old_id};
    my @line = ($old_id, $new_id, "", "release_name", "date");
    print $map_fh join("\t", @line) . "\n";
  }
  close $map_fh;
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }


  $help .= <<"EOF";
    Allocate stable ids for all genes in a given core db, from an OSID system.

    --registry <path> : Ensembl registry
    --osid_url <str>
    --osid_user <str>
    --osid_pass <str> : OSID connection details
    --species <str>   : production_name of one species

    --organism <str>  : species name in OSID, if it is not the production_name
    
    --update          : Do the actual changes (default is no OSID call and no db changes)

    Selection of features to allocate:
    --prefix <str>    : Only replace ids for genes with this prefix [optional]
    --after_date <str>: Only replace ids for genes created after this date (ISO format) [optional]
    
    Selection from a list:
    --gene_list <path>: List of gene_ids to replace from a file.
    --transcript_list <path>: List of transcript ids to replace from a file (need gene_id as first column).

    Metadata:
    --xref_source <str>: Keep the old ids under this xref name, e.g. "RefSeq" [optional]
    --analysis_name <str>: Name of the analysis to use for the renamed xrefs (default: $default_analysis_name)
    
    --mock_osid       : Get IDs from the Fake OSID server (for dev testing with update)

    --output_map <path>: Save the feature mapping in a tab file (old_id, new_id, "", "release_name", "date"). You will need to replace the last 2.

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
    "mock_osid",
    "prefix:s",
    "after_date:s",
    "xref_source:s",
    "analysis_name:s",
    "output_map:s",
    "gene_list:s",
    "transcript_list:s",
    "help",
    "verbose",
    "debug",
  );
  $opt{analysis_name} = $default_analysis_name if not $opt{analysis_name};

  usage("Registry needed") if not $opt{registry};
  usage("Species needed") if not $opt{species};
  usage("Output mapping file needed") if not $opt{output_map};
  $opt{organism} //= $opt{species};
  if (not $opt{mock_osid}) {
    usage("OSID details needed") if not ($opt{osid_url} and $opt{osid_user} and $opt{osid_pass});
  }
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{registry};
$registry->load_all($reg_path);

my $osid;
if ($opt{update} and not $opt{mock_osid}) {
  $osid = OSID_service->new(
    url  => $opt{osid_url},
    user => $opt{osid_user},
    pass => $opt{osid_pass},
  );
} else {
  $osid = OSID_service_dev->new();
}
$osid->connect($opt{organism});

my $transcript_list = load_transcript_list($opt{transcript_list});
my $genes_list = load_gene_list($opt{gene_list});
allocate_genes($osid, $registry, $opt{species}, $opt{update}, $opt{xref_source}, $opt{analysis_name}, $opt{prefix}, $opt{after_date}, $opt{output_map}, $genes_list, $transcript_list);

__END__

