#!/usr/env perl
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
# MAIN
# Get command line args
my %opt = %{ opt_check() };

my ($file1, $file2) = ($opt{in1}, $opt{in2});

if ($file1 =~ /\.fa(sta)?$/ and $file2 =~ /\.fa(sta)?$/) {
  diag "Compare fasta files";
  compare_fasta($file1, $file2);
} elsif ($file1 =~ /manifest.json$/ and $file2 =~ /manifest.json$/) {
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

  my $dna_map;
  if ($opt->{map}) {
    $dna_map = get_map($manifest1, $manifest2, $dir1, $dir2);
    #say "Mapped DNA sequences: " . scalar(keys %$dna_map);
  }

  for my $name (sort keys %$manifest1) {
    my $file1 = $manifest1->{$name};
    my $file2 = $manifest2->{$name};

    if ($file1 and $file2) {
      my $path1 = catfile($dir1, $file1);
      my $path2 = catfile($dir2, $file2);

      if ($name =~ /^fasta/ and $opt->{do_fasta}) {
        diag "Compare files $file1 and $file2";
        my $is_dna = ($name =~ /dna/);
        compare_fasta($path1, $path2, $is_dna);
      }
      if ($name eq 'gff3' and $opt->{do_gff3}) {
        diag "Compare files $file1 and $file2";
        compare_gff3($path1, $path2, $dna_map);
      }
      if ($opt->{do_json}) {
        if ($name eq 'genome') {
          diag "Compare files $file1 and $file2";
          my $data1 = get_sorted_json($path1);
          my $data2 = get_sorted_json($path2);
          #delete $data2->{assembly}->{name};
          #delete $data2->{species}->{division};
          #delete $data2->{species}->{production_name};

          my $deep = 1;
          compare_json($data1, $data2, $deep);
        }
        if ($name eq 'seq_region') {
          diag "Compare files $file1 and $file2";
          my $deep = 0;
          my $data1 = get_sorted_json($path1);
          my $data2 = get_sorted_json($path2);
          $data1 = remove_keys($data1, "BRC4_seq_region_name", "EBI_seq_region_name");
          $data2 = remove_keys($data2, "BRC4_seq_region_name", "EBI_seq_region_name");
          compare_entries($data1, $data2, "name");
        }
        if ($name eq 'functional_annotation') {
          diag "Compare files $file1 and $file2";
          my $deep = 0;
          my $data1 = get_sorted_json($path1);
          my $data2 = get_sorted_json($path2);
          #$data1 = remove_keys($data1, "xrefs", "version", "is_pseudogene");
          #$data2 = remove_keys($data2, "xrefs", "version", "is_pseudogene");
          $data1 = remove_keys($data1, "version");
          $data2 = remove_keys($data2, "version");
          compare_entries($data1, $data2, "id");
        }
      }
    }
  }
  done_testing();
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
  my ($m1, $m2, $dir1, $dir2) = @_;

  # First get the fasta files
  my ($fasta1) = grep { $_ eq 'fasta_dna' } keys %$m1;
  my ($fasta2) = grep { $_ eq 'fasta_dna' } keys %$m2;

  die "No dna fasta in " . join(", ", sort keys %$m1) if not $fasta1;
  die "No dna fasta in " . join(", ", sort keys %$m2) if not $fasta2;

  my $seqs1 = get_seqs(catfile($dir1, $m1->{$fasta1}));
  my $seqs2 = get_seqs(catfile($dir2, $m1->{$fasta2}));

  # Match the ids
  my %map;
  for my $seqh (keys %$seqs1) {
    if (exists $seqs2->{$seqh}) {
      my @ids1 = sort @{ $seqs1->{$seqh} };
      my @ids2 = sort @{ $seqs2->{$seqh} };
      if (@ids1 == 1 and @ids2 == 1) {
        $map{$ids2[0]} = $ids1[0];
      } else {
        for (my $i = 0; $i < @ids1; $i++) {
          $map{$ids2[$i]} = $ids1[$i];
        }
      }
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
  my ($fasta1, $fasta2, $is_dna) = @_;

  fasta_diff($fasta1, $fasta2, $is_dna);
}

sub fasta_diff {
  my ($f1, $f2, $is_dna) = @_;

  my $seqs1 = get_fasta($f1);
  my $seqs2 = get_fasta($f2);

  my $n1 = scalar(keys %$seqs1);
  my $n2 = scalar(keys %$seqs2);
  cmp_ok($n1, '>', 0, "Fasta 1 has sequences ($n1)");
  cmp_ok($n2, '>', 0, "Fasta 1 has sequences ($n2)");
  cmp_ok($n1, '==', $n2, "Fasta 1 and 2 have the same number of sequences");

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
  cmp_ok($nonly1, '==', 0, "Fasta 1 has no sequence specific to it ($nonly1)");
  cmp_ok($nonly2, '==', 0, "Fasta 2 has no sequence specific to it ($nonly2)");
  cmp_ok($count_different_length, '==', 0, "All sequences are the same length ($count_different_length different)");
  cmp_ok($count_different_seq, '==', 0, "All sequences are the same sequence, expect for ambiguous letters ($count_different_seq different)");

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
  my ($fasta) = @_;

  my %seqs;
  my ($seq, $id);
  open my $fh, "<", $fasta;
  while (my $line = readline($fh)) {
    chomp $line;
    if ($line =~ /^>\s*(.+)\s*$/) {
      if ($id) {
        $seqs{$id} = $seq;
      }
      $id = $1;
      $seq = '';
    } elsif (not $line =~ /^\s*$/) {
      $seq .= $line;
    }
  }
  if ($id) {
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
  my ($json, $delete_version) = @_;

  if (ref($json) eq '') {
    return $json;
  } elsif (ref($json) eq 'HASH') {
    for my $key (keys %$json) {
      if ($key eq 'version' and $delete_version) {
        delete $json->{$key};
      } else {
        $json->{$key} = sort_json($json->{$key}, $delete_version);
      }
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
        } elsif ($json->[0]->{name}) {
          $json = [ sort { $a->{name} cmp $b->{name} } @$json ];
        }
      }

      # Sort each element in the array itself
      for (my $i = 0; $i < @$json; $i++) {
        $json->[$i] = sort_json($json->[$i], $delete_version);
      }
    }
    return $json;
  } else {
    return $json;
  }
}

sub compare_json {
  my ($data1, $data2, $deep) = @_;

  if ($deep) {
    eq_or_diff($data1, $data2, "Json content identical");
  } else {
    cmp_deeply($data1, $data2, "Json content identical");
  }
}

sub remove_keys {
  my ($data, @keys) = @_;
  
  for my $entry (@$data) {
    for my $key (@keys) {
      delete $entry->{$key} if $entry->{$key};
    }
  }
  return $data;
}

sub compare_entries {
  my ($data1, $data2, $key) = @_;
  
  my $ndata1 = scalar @$data1;
  my $ndata2 = scalar @$data2;
  cmp_ok($ndata1, '==', $ndata2, "Same number of entries ($ndata1 vs $ndata2)");

  # Extract entries that are in common
  ($data1, $data2) = list_diff($data1, $data2, $key);

  # Compare those entries that are in common
  my @sorted1 = sort { $a->{$key} cmp $b->{$key} } @$data1;
  my @sorted2 = sort { $a->{$key} cmp $b->{$key} } @$data2;

  for (my $i = 0; $i < @sorted1; $i++) {
    my $en1 = $sorted1[$i];
    my $en2 = $sorted2[$i];
    my $value = "$en1->{$key}";
    $value .= "/$en2->{$key}" if $en1->{$key} ne $en2->{$key};
    cmp_deeply($en1, $en2, "Json content identical for $value");
  }
}

sub list_diff {
  my ($data1, $data2, $key) = @_;
  
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

  cmp_ok($nfile1, "==", 0, "No entry found only in file 1 ($nfile1)");
  for my $id (sort keys %ids1) {
    diag "ONLY in file 1: $id";
  }

  cmp_ok($nfile2, "==", 0, "No entry found only in file 2 ($nfile2)");
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
  my ($gff1, $gff2, $dna_map) = @_;

  my $data1 = get_gff3($gff1);
  my $data2 = get_gff3($gff2, $dna_map);

  # Compare biotypes counts
  say "Differences:";
  my %stats = (
    same => [],
    diff => [],
    only1 => [],
    only2 => [],
  );

  my (@all_only1, @all_only2);
  
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

      delete $data2->{$type};

      # Check features coordinates
      my $coords_diff = 0;
      for my $id (keys %$entries1) {
        my $feat1 = $entries1->{$id};
        my $feat2 = $entries2->{$id};
        next if not $feat2;
        if (
               $feat1->{chr} ne $feat2->{chr}
            or $feat1->{start} != $feat2->{start}
            or $feat1->{end}   != $feat2->{end}
        ) {
          $coords_diff++;
          say(sprintf("Different (only first one shown): %s:%d-%d vs %s:%d-%d", ($feat1->{chr}, $feat1->{start}, $feat1->{end}, $feat2->{chr}, $feat2->{start}, $feat2->{end}))) if $coords_diff == 1;
        }
      }

      # Check length count
      if ($count1 == $count2) {
        my $line = "\tSAME\t$type ($count1)";
        $line .= " (diff coord: $coords_diff)" if $coords_diff;
        push @{$stats{same}}, $line;
      } else {
        # List differences
        my ($only1, $only2, $redundant) = diff_entries($entries1, $entries2);

        if (@$only1 or @$only2) {
          push @all_only1, map { "$type\t$_" } @$only1;
          push @all_only2, map { "$type\t$_" } @$only2;

          my $line = "\tDIFF\t$type ($count1 != $count2)";
          $line .= " (diff coord: $coords_diff)" if $coords_diff;
          push @{$stats{diff}}, $line;
        } else {
          my $line = "\tSAME\t$type ($count1)";
          $line .= " (diff coord: $coords_diff)" if $coords_diff;
          $line .= " ($redundant redundant)" if $redundant;
          push @{$stats{same}}, $line;
        }
      }
    }
    else {
      push @{$stats{only1}}, "\tONLY1\t$type ($count1)";
      push @all_only1, map { "$type\t$_" } sort keys %$entries1;
    }
    delete $data1->{$type};
  }
  
  for my $type (keys %$data2) {
    my $entries2 = $data2->{$type};
    my $count2 = scalar(keys %$entries2);
    push @{$stats{only2}}, "\tONLY2\t$type ($count2)";
    push @all_only2, map { "$type\t$_" } sort keys %$entries2;
  }

  for my $name (qw(same diff only1 only2)) {
    my @lines = @{ $stats{$name} };

    for my $line (@lines) {
      say $line;
    }
  }

  @all_only1 = sort @all_only1;
  @all_only2 = sort @all_only2;

  say "Details, only in 1:" if @all_only1;
  for my $str (@all_only1) { say "\t$str" }

  say "Details, only in 2:" if @all_only2;
  for my $str (@all_only2) { say "\t$str" }

  return;
}

sub diff_entries {
  my($e1, $e2) = @_;
  
  # Compare by coords, not id
  my $redundant;
  my (%c1, %c2);
  for my $id (keys %$e1) {
    my $e = $e1->{$id};
    my $coord = join(" ", map { $e->{$_} } qw(chr start end));
    if ($c1{$coord}) {
      #say STDERR "Already a feature at $coord: $c1{$coord} (new: $id)";
      $redundant++;
    }
    $c1{$coord} = $id;
  }
  for my $id (keys %$e2) {
    my $e = $e2->{$id};
    my $coord = join(" ", map { $e->{$_} } qw(chr start end));
    $c2{$coord} = $id;
  }
  
  for my $coord (keys %c1) {
    if (exists $c2{$coord}) {
      delete $c2{$coord};
      delete $c1{$coord};
    }
  }

  my @only1 = map { "$_\t$c1{$_}" } sort keys %c1;
  my @only2 = map { "$_\t$c2{$_}" } sort keys %c2;
  
  my $c = @only1 + @only2;
  
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

    # Remove small UTRs like UPenn
    if (($type eq 'five_prime_UTR' or $type eq 'three_prime_UTR') and $end - $start < 3) {
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
  $help .= <<'EOF';
    SCRIPTDESCRIPTION
    
    --in1 <path>
    --in2 <path>
    
    --do_fasta
    --do_gff3
    --do_json

    --map
    
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
    "map",
    "help",
    "verbose",
    "debug",
  ) or usage();

  if (not ($opt{do_fasta} or $opt{do_gff3} or $opt{do_json})) {
    $opt{do_fasta} = 1;
    $opt{do_gff3} = 1;
    $opt{do_json} = 1;
  }
  usage()                if $opt{help};
  usage("Manifest 1 needed")  if not $opt{in1};
  usage("Manifest 2 needed")  if not $opt{in2};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}
__END__


