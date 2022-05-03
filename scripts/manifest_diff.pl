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
use JSON;
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 
use File::Basename qw(dirname);
use File::Spec::Functions qw(catfile);
use Test::More;
use Test::Deep;
use Test::Differences;

###############################################################################
my $seqr_key_default = "name";

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my ($file1, $file2) = ($opt{in1}, $opt{in2});

if ($file1 =~ /manifest.json$/ and $file2 =~ /manifest.json$/) {
  diag "Compare manifests";
  compare_manifests($file1, $file2, \%opt);
} else {
  say "Unrecognized files compared";
}

sub compare_manifests {
  my ($man1, $man2, $opt) = @_;
  
  my $dir1 = dirname($man1);
  my $dir2 = dirname($man2);
  my $manifest1 = get_manifest($man1);
  my $manifest2 = get_manifest($man2);

  # Get seq_region name mapping
  my $map1 = get_map($manifest1, $dir1, $opt{seqr_key});
  my $map2 = get_map($manifest2, $dir2, $opt{seqr_key});

  for my $name (sort keys %$manifest1) {
    my $file1 = $manifest1->{$name};
    my $file2 = $manifest2->{$name};

    if ($file1 and $file2) {
      my $path1 = catfile($dir1, $file1);
      my $path2 = catfile($dir2, $file2);

      if ($name =~ /^fasta/ and $opt->{do_fasta}) {
        diag "Compare files $file1 and $file2";
        my $is_dna = ($name =~ /dna/);
        compare_fasta($path1, $path2, $is_dna, $map1, $map2);
      }
      if ($name eq 'gff3' and $opt->{do_gff3}) {
        diag "Compare files $file1 and $file2";
        compare_gff3($path1, $path2, $map1, $map2);
      }
      if ($opt->{do_json}) {
        if ($name eq 'genome') {
          diag "Compare files $file1 and $file2";
          my $data1 = get_sorted_json($path1);
          my $data2 = get_sorted_json($path2);

          # Remove quite a lot of keys
          $data1 = clean_genome_meta($data1);
          $data2 = clean_genome_meta($data2);

          eq_or_diff($data1, $data2, "Json content identical");
        }
        if ($name eq 'seq_region') {
          diag "Compare files $file1 and $file2";
          my $data1 = get_sorted_json($path1, $opt{seqr_key});
          my $data2 = get_sorted_json($path2, $opt{seqr_key});
          $data1 = remove_keys($data1, ["EBI_seq_region_name"]);
          $data2 = remove_keys($data2, ["EBI_seq_region_name"]);
          compare_entries($data1, $data2, $opt{seqr_key});
        }
        if ($name eq 'functional_annotation') {
          diag "Compare files $file1 and $file2";
          my $data1 = get_sorted_json($path1);
          my $data2 = get_sorted_json($path2);
          $data1 = remove_keys($data1, ["xrefs", "version"]);
          $data2 = remove_keys($data2, ["xrefs", "version"]);
          $data1 = unique_array_keys($data1, "xrefs", "id", "dbname");
          $data2 = unique_array_keys($data2, "xrefs", "id", "dbname");
          compare_entries($data1, $data2, "id");
        }
      }
    }
  }
  done_testing();
}

sub clean_genome_meta {
  my ($data) = @_;

  delete $data->{annotation};
  delete $data->{species}->{division};
  delete $data->{species}->{production_name};
  delete $data->{assembly}->{name};

  delete $data->{BRC4};
  delete $data->{species}->{BRC4_organism_abbrev};

  $data = move_provider($data, "name");
  $data = move_provider($data, "url");

  return $data;
}

sub move_provider {
  my ($data, $key) = @_;

  if ($data->{assembly}->{"provider_".$key} and not $data->{provider}->{$key}) {
    $data->{provider}{$key} = $data->{assembly}->{"provider_".$key};
    delete $data->{assembly}->{"provider_".$key};
  }

  return $data;
}

sub get_manifest {
  my ($man_path) = @_;

  my $json = get_json($man_path);
  my %data;
  for my $file (keys %$json) {
    $data{$file} = $json->{$file}->{file};
  }
  return \%data;
}

sub get_map {
  my ($manifest, $dir, $key) = @_;
  
  my %map;
  my $seqr_file = $manifest->{seq_region};
  if ($seqr_file) {
    my $seqr_path = catfile($dir, $seqr_file);
    if (-e $seqr_path) {
      my $seqr = get_json($seqr_path);
      %map = map { $_->{name} => $_->{$key} } grep { exists $_->{$key} } @$seqr;
    }
  }

  return \%map;
}

sub get_seqs {
  my ($fasta) = @_;

  my %seqs;
  my ($seq, $id);
  open my $fh, "<", $fasta;
  while (my $line = readline($fh)) {
    chomp $line;
    if ($line =~ /^>\s*(.+)\s*$/) {
      if ($id) {
        push @{ $seqs{$seq} }, $id;
      }
      $id = $1;
      $seq = '';
    } elsif (not $line =~ /^\s*$/) {
      $seq .= $line;
    }
  }
  if ($id) {
    push @{ $seqs{$seq} }, $id;
  }
  close $fh;
  
  return \%seqs;
}

sub compare_fasta {
  my ($f1, $f2, $is_dna, $map1, $map2) = @_;

  if (not $is_dna) {
    $map1 = {};
    $map2 = {};
  }
  my $seqs1 = get_fasta($f1, $map1);
  my $seqs2 = get_fasta($f2, $map2);

  my $n1 = scalar(keys %$seqs1);
  my $n2 = scalar(keys %$seqs2);
  cmp_ok($n1, '>', 0, "Fasta 1 has sequences ($n1)");
  cmp_ok($n2, '>', 0, "Fasta 2 has sequences ($n2)");

  my $count_ambiguous = 0;
  my $count_different_seq = 0;
  my $count_different_length = 0;
  my $count_identical = 0;

  for my $id (sort keys %$seqs1) {
    if ($seqs2->{$id}) {
      my $seq1 = $seqs1->{$id};
      my $seq2 = $seqs2->{$id};

      # Diff lengths
      if (length($seq1) != length($seq2)) {
          $count_different_length++;
          diag("$id different lengths: " . length($seq1) . " vs " . length($seq2));
      } elsif ($seq1 ne $seq2) {
        if ($is_dna) {
          # Check if the difference is only for ambiguous chars
          my $n_ambiguous = diff_is_ambiguous($seq1, $seq2);
          if ($n_ambiguous) {
            $count_ambiguous++;
          } else {
            $count_different_seq++;
            diag("$id different sequences");
          }
        }
      } else {
        $count_identical++;
      }

      delete $seqs1->{$id};
      delete $seqs2->{$id};
    } else {
      diag("$id is only in fasta 1");
    }
  }
  for my $id (sort keys %$seqs2) {
      diag("$id is only in fasta 2");
  }

  my $nonly1 = scalar(keys %$seqs1);
  my $nonly2 = scalar(keys %$seqs2);
  cmp_ok($nonly1, '==', 0, "Fasta 1 has no sequence specific to it");
  cmp_ok($nonly2, '==', 0, "Fasta 2 has no sequence specific to it");
  cmp_ok($count_different_length, '==', 0, "All sequences are the same length");
  cmp_ok($count_different_seq, '==', 0, "All sequences are the same sequence, expect for ambiguous letters");

  #say("$count_identical identical sequences") if $count_identical;
  #say("$count_ambiguous identical except for ambiguous nucleotides") if $count_ambiguous;
}

sub diff_is_ambiguous {
  my ($seq1, $seq2) = @_;

  # Check length
  return if length($seq1) != length($seq2);
  my @seq1 = split(//, $seq1);
  my @seq2 = split(//, $seq2);

  # Compare each amino acid
  my $n_ambiguous = 0;
  for (my $i = 0; $i < length($seq1); $i++) {
    my $aa1 = $seq1[$i];
    my $aa2 = $seq2[$i];

    if ($aa1 ne $aa2 and is_ambiguous($aa1) and is_ambiguous($aa2)) {
      $n_ambiguous++;
    }
  }

  return $n_ambiguous;
}

sub is_ambiguous {
  my ($aa) = @_;

  return not $aa =~ /[CGTA]/i;
}

sub get_fasta {
  my ($fasta, $map) = @_;

  my %seqs;
  my ($seq, $id);
  open my $fh, "<", $fasta;
  while (my $line = readline($fh)) {
    chomp $line;
    if ($line =~ /^>\s*(.+)\s*$/) {
      if ($id) {
        $id = $map->{$id} if $map->{$id};
        $seqs{$id} = $seq;
      }
      $id = $1;
      $seq = '';
    } elsif (not $line =~ /^\s*$/) {
      $seq .= $line;
    }
  }
  if ($id) {
    $id = $map->{$id} if $map->{$id};
    $seqs{$id} = $seq;
  }
  close $fh;
  
  return \%seqs;
}

sub get_sorted_json {
  my ($json_path, $delete_version) = @_;

  my $json = get_json($json_path);
  
  # Sort arrays
  $json = sort_json($json, $delete_version);

  return $json;
}

sub get_json {
  my ($json_path) = @_;

  my $json_string;
  {
    local $/; #Enable 'slurp' mode
    open(my $fh, "<", "$json_path");
    $json_string = <$fh>;
    close $fh;
  }
  die("No json data") if not $json_string;

  my %data;
  my $json = decode_json($json_string);
}

sub sort_json {
  my ($json, $key) = @_;

  $key //= "name";

  if (ref($json) eq '') {
    return $json;
  } elsif (ref($json) eq 'HASH') {
    for my $key (keys %$json) {
      $json->{$key} = sort_json($json->{$key});
    }
    return $json;
  } elsif (ref($json) eq 'ARRAY') {
    my @new_array = ();

    # Sort array
    if (@$json) {
      if (ref($json) eq '') {
        $json = [sort @$json];
      } elsif (ref($json) eq 'HASH') {
        # Sort for functional_annotation
        if ($json->[0]->{object_type} and $json->[0]->{id}) {
          $json = [ sort { $a->{object_type} cmp $b->{object_type} and $a->{id} cmp $b->{id} } @$json ];
        # Sort for others
        } elsif ($json->[0]->{$key}) {
          $json = [ sort { $a->{$key} cmp $b->{$key} } @$json ];
        }
      }

      # Sort each element in the array itself
      for (my $i = 0; $i < @$json; $i++) {
        $json->[$i] = sort_json($json->[$i]);
      }
    }
    return $json;
  } else {
    return $json;
  }
}

sub remove_keys {
  my ($data, $keys) = @_;
  
  for my $entry (@$data) {
    for my $key (@$keys) {
      delete $entry->{$key} if $entry->{$key};
    }
  }
  return $data;
}

sub unique_array_keys {
  my ($data, $key, @tags) = @_;

  for my $entry (@$data) {
    my $array = $entry->{$key};
    if ($array and ref($array) eq 'ARRAY') {
      my %done;
      my @new_array;
      for my $item (@$array) {
        my $key = join(':', map { $item->{$_} } @tags);
        push @new_array, $item if not $done{$key}++;
      }
      $entry->{$key} = \@new_array;
      my $diff = scalar(@$array) - scalar(keys %done);
    }
  }
  return $data;
}

sub compare_entries {
  my ($data1, $data2, $key) = @_;
  
  # Extract entries that are in common
  ($data1, $data2) = list_diff($data1, $data2, $key);

  # Compare those entries that are in common
  my @sorted1 = sort { $a->{$key} cmp $b->{$key} } @$data1;
  my @sorted2 = sort { $a->{$key} cmp $b->{$key} } @$data2;

  for (my $i = 0; $i < @sorted1; $i++) {
    my $en1 = $sorted1[$i];
    my $en2 = $sorted2[$i];

    my $value1 = $en1->{$key} ? "$en1->{$key}" : $en1->{name};
    my $value2 = $en2->{$key} ? "$en2->{$key}" : $en2->{name};
    my $value = "$value1/$value2";

    # Special
    if ($en1->{description}) {
      $en1->{description} =~ s/ \[Source:VB Community Annotation\]//;
    }
    if ($en2->{description}) {
      $en2->{description} =~ s/ \[Source:VB Community Annotation\]//;
    }

    eq_or_diff($en1, $en2, "Json content identical for $value");
  }
}

sub list_diff {
  my ($data1, $data2, $key) = @_;

  die "No $key in data1" if not $data1->[0]->{$key};
  die "No $key in data2" if not $data2->[0]->{$key};
  
  my %ids1 = map { $_->{$key} => 1 } @$data1;
  my %ids2 = map { $_->{$key} => 1 } @$data2;
  
  my %common;
  for my $id (keys %ids1) {
    if ($ids2{$id}) {
      $common{$id} = 1;
      delete $ids1{$id};
      delete $ids2{$id};
    }
  }
  
  my $nfile1 = scalar(keys %ids1);
  my $nfile2 = scalar(keys %ids2);

  cmp_ok($nfile1, "==", 0, "No entry found only in file 1");
  for my $id (sort keys %ids1) {
    diag "ONLY in file 1: $id";
  }

  cmp_ok($nfile2, "==", 0, "No entry found only in file 2");
  for my $id (sort keys %ids2) {
    diag "ONLY in file 2: $id";
  }

  $data1 = [grep { $common{ $_->{$key} } } @$data1];
  $data2 = [grep { $common{ $_->{$key} } } @$data2];
  return ($data1, $data2);
}

######################################################################
sub merge_types {
  my  ($data, $main_type, $merged_type) = @_;

  if ($data->{$main_type}) {
    if ($data->{$merged_type}) {
      $data->{$main_type} = { %{$data->{$main_type}}, %{$data->{$merged_type}} };
      delete $data->{$merged_type};
    }
  } else {
    $data->{$main_type} = $data->{$merged_type};
    delete $data->{$merged_type};
  }
  return $data;
}

sub compare_gff3 {
  my ($gff1, $gff2, $map1, $map2) = @_;

  my $data1 = get_gff3($gff1, $map1);
  my $data2 = get_gff3($gff2, $map2);
  
  # Changes made during import
  $data1 = merge_types($data1, "gene", "ncRNA_gene");
  $data2 = merge_types($data2, "gene", "ncRNA_gene");
  $data1 = merge_types($data1, "exon", "pseudogenic_exon");
  $data2 = merge_types($data2, "exon", "pseudogenic_exon");
  

  for my $type (sort keys %$data1) {
    my $entries1 = $data1->{$type};
    my $count1 = scalar(keys %$entries1);

    if (exists $data2->{$type}) {
      my $entries2 = $data2->{$type};
      my $count2 = scalar(keys %$entries2);

      diag "$type: $count1 vs $count2";

      delete $data2->{$type};

      # Check features coordinates
      my (@coords_diff, %only1, %only2);
      for my $id (keys %$entries1) {
        my $feat1 = $entries1->{$id};
        my $feat2 = $entries2->{$id};
        if (not $feat2) {
          $only1{$id} = $feat1;
          delete $entries1->{$id};
          next;
        }
        delete $entries1->{$id};
        delete $entries2->{$id};

        if (
               $feat1->{chr} ne $feat2->{chr}
            or $feat1->{start} != $feat2->{start}
            or $feat1->{end}   != $feat2->{end}
        ) {
          push @coords_diff, sprintf("%s (%s:%d-%d\tvs\t%s:%d-%d)",
            ($id,
              $feat1->{chr},
              $feat1->{start},
              $feat1->{end},
              $feat2->{chr},
              $feat2->{start}, 
              $feat2->{end}
            ));
        }
      }
      for my $id (keys %$entries2) {
        $only2{$id} = $entries2->{$id};
        delete $entries2->{$id};
      }
      my ($only1_coord, $only2_coord, $redundant) = diff_coord_entries(\%only1, \%only2);

      @coords_diff = sort @coords_diff;
      cmp_ok(scalar(@coords_diff), "==", 0, "$type: same id = same coordinates");
      for my $str (@coords_diff) { diag("$type: different coordinates: $str") };

      # Compare by id, only for genes, transcripts and translations
      if ($type eq 'gene' or $type eq 'mRNA' or $type eq 'CDS') {
        cmp_ok(scalar(keys %only1), "==", 0, "$type: none found only in gff 1 by id");
        for my $str (sort keys %only1) { diag("$type: Only in gff1: $str") };

        cmp_ok(scalar(keys %only2), "==", 0, "$type: none found only in gff 2 by id");
        for my $str (sort keys %only2) { diag("$type: Only in gff2: $str") };

        # Compare by coordinates
        cmp_ok($redundant, "==", 0, "$type: No redundant features");
      } else {
        diag("$type: $redundant redundant features");
      }

      cmp_ok(scalar(@$only1_coord), "==", 0, "$type: none found only in gff 1");
      for my $str (@$only1_coord) { diag("$type: Only in gff1 by coordinates: $str") };

      cmp_ok(scalar(@$only2_coord), "==", 0, "$type: none found only in gff 2");
      for my $str (@$only2_coord) { diag("$type: Only in gff2 by coordinates: $str") };
    }
    else {
      fail("$type only found in gff1");
      for my $str (sort keys %$entries1) { diag("$type: Only in gff1: $str") };
    }
    delete $data1->{$type};
  }
  
  for my $type (keys %$data2) {
    my $entries2 = $data2->{$type};
    fail("$type only found in gff2");
    for my $str (sort keys %$entries2) { diag("$type: Only in gff2: $str") };
  }

  return;
}

sub diff_coord_entries {
  my($e1, $e2) = @_;
  
  # Compare by coords, not id
  my $redundant = 0;
  my (%c1, %c2);
  for my $id (keys %$e1) {
    my $e = $e1->{$id};
    my $coord = sprintf("%s:%d-%d", map { $e->{$_} } qw(chr start end));
    if ($c1{$coord}) {
      #say STDERR "Already a feature at $coord: $c1{$coord} (new: $id)";
      $redundant++;
    }
    $c1{$coord} = $id;
  }

  for my $id (keys %$e2) {
    my $e = $e2->{$id};
    my $coord = sprintf("%s:%d-%d", map { $e->{$_} } qw(chr start end));
    $c2{$coord} = $id;
  }
  
  for my $coord (keys %c1) {
    if (exists $c2{$coord}) {
      delete $c2{$coord};
      delete $c1{$coord};
    }
  }

  my @only1 = map { "$c1{$_} ($_)" } sort keys %c1;
  my @only2 = map { "$c2{$_} ($_)" } sort keys %c2;
  
  return \@only1, \@only2, $redundant;
}

sub get_gff3 {
  my ($gff_path, $map) = @_;

  my $fh;
  if ($gff_path =~ /\.gz$/) {
    open($fh, '-|', "gzip -dc $gff_path");
  } else {
    open($fh, "<", $gff_path);
  }

  my $data;
  while (my $line = readline $fh) {
    next if $line =~ '^#';
    chomp($line);

    my ($chr, $source, $type, $start, $end, $a, $b, $c, $attr) = split(/\t/, $line);
    my %attr = map { split(/=/, $_) } split(/;/, $attr);
    my $id = $attr{ID} // "$chr-$start-$end";
    my $biotype = $attr{biotype};

    $chr = $map->{$chr} if $map and $map->{$chr};

    # Actually just ignore UTRs: those are not loaded anyway
    if (($type eq 'five_prime_UTR' or $type eq 'three_prime_UTR')) {
      next;
    }

    my $entry = {
      chr => $chr,
      start => $start,
      end => $end,
    };
    
    $data->{$type}->{$id} = $entry;
  }
  return $data;
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
    Compare two sets of files given their manifest, assuming they follow BRC4 specifications.
    
    --in1 <path>      : path the the first manifest file
    --in2 <path>      : path the the second manifest file

    OPTIONS:
    --seqr_key <str>  : Key to use to map the seq_regions (default: $seqr_key_default)
    
    SPECIFIC TESTS:
    --do_fasta
    --do_gff3
    --do_json

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
    "in1=s",
    "in2=s",
    "do_fasta",
    "do_gff3",
    "do_json",
    "seqr_key=s",
    "help",
    "verbose",
    "debug",
  ) or usage();

  if (not ($opt{do_fasta} or $opt{do_gff3} or $opt{do_json})) {
    $opt{do_fasta} = 1;
    $opt{do_gff3} = 1;
    $opt{do_json} = 1;
  }
  $opt{seqr_key} //= $seqr_key_default;
  usage()                if $opt{help};
  usage("Manifest 1 needed")  if not $opt{in1};
  usage("Manifest 2 needed")  if not $opt{in2};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

