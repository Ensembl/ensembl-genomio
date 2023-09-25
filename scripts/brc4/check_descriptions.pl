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

use Data::Dumper;

###############################################################################
# MAIN
my $sus_threshold = 100;
main();

sub main {
  # Get command line args
  my %opt = %{ opt_check() };

  $logger->info("Load registry");
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($opt{registry}, 1);

  my @header = qw(
  build
  component
  orgabbrev
  species
  genes
  genes_desc
  genes_xref
  tr_transferable
  tr_putative
  tr_untransferable  
  );
  if ($opt{translations}) {
    my @tl_header = qw(
    tl_transferable
    tl_putative
    tl_untransferable  
    );
    @header = (@header, @tl_header);
  }
  say("#" . join("\t", @header));

  my @all_species = ($opt{species}) || @{$registry->get_all_species()};
  for my $species (sort @all_species) {
    my $ma = $registry->get_adaptor($species, "core", "MetaContainer");
    my $component = get_meta_value($ma, 'BRC4.component');
    if ($opt{component} and $opt{component} ne $component) {
      $ma->dbc->disconnect_if_idle();
      next;
    }

    my $count = check_genes($registry, $species);
    my $build = get_build($ma, $species);
    my $org = get_meta_value($ma, 'BRC4.organism_abbrev');

    $ma->dbc->disconnect_if_idle();

    my @line = (
      $build,
      $component,
      $org,
      $species,
      $count->{gene_total},
      $count->{gene_full},
      $count->{gene_xref},

      $count->{empty_tr_transferable},
      $count->{empty_tr_putative},
      $count->{empty_tr_untransferable},
    );

    if ($opt{translations}) {
      my @tl_line = (
        $count->{empty_tl_transferable},
        $count->{empty_tl_putative},
        $count->{empty_tl_untransferable},
      );
      @line = (@line, @tl_line);
    }
    say join("\t", @line);    
  }
}

sub get_build {
  my ($ma, $key) = @_;

  my $dbname = $ma->dbc->dbname;
  if ($dbname =~ /_(\d+)_\d+_\d+$/) {
    return $1;
  }
}

sub get_meta_value {
  my ($ma, $key) = @_;

  my ($value) = @{ $ma->list_value_by_key($key) };

  return $value;
}

sub check_genes {
  my ($registry, $species, $check_translations) = @_;
  
  my $ga = $registry->get_adaptor($species, "core", "gene");
  my $dbname = $ga->dbc->dbname;
  $logger->info("Database:\t$dbname");
  $logger->info("Species:\t$species");
  
  my %count = (
    empty_tr_putative => 0,
    empty_tr_transferable => 0,
    empty_tr_untransferable => 0,
    gene_and_transcript_empty => 0,
    all_empty => 0,
    gene_full => 0,
    gene_xref => 0,
    gene_empty => 0,
    gene_total => 0,
  );

  if ($check_translations) {
    my %tl_count = (
    empty_tl_putative => 0,
    empty_tl_transferable => 0,
    empty_tl_untransferable => 0,
    );
    %count = (%count, %tl_count);
  }
  
  # To check some names are not highly repeated (e.g. hypothetical protein)
  my %tname;
  my %tlname;
  
  GENE: for my $gene (@{$ga->fetch_all()}) {
    my $stable_id = $gene->stable_id;
    
    $count{gene_total}++;
    my $g_status = check_description_status($gene->description);
    if ($g_status->{full}) {
      $count{gene_full}++;
      next GENE;
    } elsif ($g_status->{xref}) {
      $count{gene_xref}++;
    } else {
      $count{gene_empty}++;
    }
    
    my @transcripts = @{$gene->get_all_Transcripts()};

    # Check transcript and translation description
    my $tdesc = "";
    my $tldesc = "";
    for my $transc (@transcripts) {
      my $cur_tdesc = $transc->description;
      
      # Some genomes have the same description for different transcripts in the same gene, but with ", variant" added
      $cur_tdesc =~ s/,( transcript)? variant( \d+| [A-Z]\d*)?$// if $cur_tdesc;
      my $cur_check = check_description_status($cur_tdesc);
      if ($cur_check->{putative}) {
        $cur_tdesc = "hypothetical protein";
      }

      if (not $tdesc) {
        $tdesc = $cur_tdesc;
      } elsif ($cur_tdesc and lc($tdesc) ne lc($cur_tdesc)) {
        $logger->info("Gene $stable_id has several transcript descriptions:\n\t'$tdesc'\n\t'".$cur_tdesc ."'");
        $count{empty_tr_untransferable}++;
        next GENE;
      }
      
      if ($check_translations) {
        # Translation product from Uniprot
        my $cur_tldesc = get_translation_product($transc);
        
        if (not $tldesc) {
          $tldesc = $cur_tldesc;
        } elsif ($cur_tldesc and lc($tldesc) ne lc($cur_tldesc)) {
          $logger->info("Gene $stable_id has several translation products:\n\t'$tldesc'\n\t'".$cur_tldesc ."'");
          $count{empty_tl_untransferable}++;
          next GENE;
        }
      }
    }
    
    if ($tdesc) {
      $tname{$tdesc}++;
      my $t_status = check_description_status($tdesc);
      
      if ($t_status->{putative}) {
        $count{empty_tr_putative}++;
      } else {
        $count{empty_tr_transferable}++;
      }
    } elsif ($check_translations and $tldesc) {
      $tlname{$tldesc}++;
      my $tl_status = check_description_status($tldesc);
      
      if ($tl_status->{putative}) {
        $count{empty_tl_putative}++;
      } else {
        $count{empty_tl_transferable}++;
      }
    } else {
      $count{all_empty}++;
    }
  }
  
  check_repeated_names(\%tname);
  
  for my $stat (sort keys %count) {
    $logger->info(sprintf("\t%7d %s", $count{$stat}, $stat));
  }
  
  return \%count;
}

sub get_translation_product {
  my ($transcript) = @_;
  
  my $translation = $transcript->translation;
  return if not $translation;
  
  my @dblinks = @{ $translation->get_all_DBEntries('Uniprot%') };
  my %uniq_descs = map { $_->description => 1 } grep { $_->description } @dblinks;
  my @descs = grep { $_ } sort keys %uniq_descs;
  
  if (@descs == 1) {
    return shift @descs;
  } elsif (@descs > 1) {
    my $tr_id = $transcript->stable_id;
    my $descs_str = join " | ", @descs;
    $logger->debug("Warning: several Uniprot products for this transcript '$tr_id': $descs_str\n");
    return shift @descs;
  }
}

sub check_repeated_names {
  my ($name) = @_;
  
  my @sus = grep { $name->{$_} > $sus_threshold } keys %$name;
  for my $desc (sort { $name->{$b} <=> $name->{$a} } @sus) {
    my $status = check_description_status($desc);
    if (not $status->{putative}) {
      $logger->info("The transcript description '$desc' is repeated $name->{$desc} times");
    }
  }
}

sub check_description_status {
  my ($desc) = @_;

  my %status = (
    full => 0,
    empty => 0,
    xref => 0,
    putative => 0,
  );
  
  if ($desc) {
    if ($desc =~ /\[Source:/) {
      $status{xref} = 1;
    } elsif ($desc =~ /(conserved)? *(hypothetical|uncharacterized|predicted|putative)( *conserved)?[ _]+protein(, conserved)?|protein of unknown function/i) {
      $status{putative} = 1;
    } else {
      $status{full} = 1;
    }
  } else {
    $status{empty} = 1;
  }
  
  return \%status;
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
    Check genes and transcript descriptions in a core

    --registry <path> : Ensembl registry for the core database

    Optional:
    --species <str>   : production_name of one species
    --component <str> : restrict check to species of this component
    --translations    : check how many genes without descriptions from transcript could use a translation Uniprot xref
    
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
    "component=s",
    "translations",
    "help",
    "verbose",
    "debug",
  );

  usage("Registry needed") if not $opt{registry};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

