=head1 LICENSE

Copyright [2009-2015] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Pipeline::Runnable::BRC4::DumpAGP;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');
use File::Spec::Functions qw(catdir catfile);

sub run {
  my ($self) = @_;
  my $species      = $self->param_required('species');
  
  my $gap_type = 'scaffold';
  my $linkage  = 'yes';
  my $evidence = 'paired-ends';
  
  
  my $dba = $self->core_dba();
  my $ama = $dba->get_adaptor('AssemblyMapper');
  my $csa = $dba->get_adaptor('CoordSystem');
  my $sa  = $dba->get_adaptor('Slice');

  my @coord_maps = $self->get_coord_maps($dba);

  my $sub_dir = $self->create_dir('agp');
  my %agp_files;
  foreach my $pair (@coord_maps) {
    my $first_level_cs = $csa->fetch_by_dbID($pair->[0]);
    my $second_level_cs = $csa->fetch_by_dbID($pair->[1]);
  
    my $mapper = $ama->fetch_by_CoordSystems($first_level_cs, $second_level_cs);
    my $slices = $sa->fetch_all($first_level_cs->name, $first_level_cs->version);

    my $map_name = $first_level_cs->name() . "-" . $second_level_cs->name();
    my $agp_file = catfile($sub_dir, $self->production_name() . '_' . $map_name . '.agp');
    $agp_files{$map_name} = $agp_file;
    open(my $out_fh, '>', $agp_file) or $self->throw("Cannot open file $agp_file: $!");
    foreach my $slice (sort {$a->seq_region_name cmp $b->seq_region_name} @$slices) {
      my @seq_level_coords =
        $mapper->map(
          $slice->seq_region_name,
          $slice->start,
          $slice->end,
          $slice->strand,
          $first_level_cs
        );
      
      my $asm_start = 1;
      my $cmp_count = 0;
      
      foreach my $seq_level_coord (@seq_level_coords) {
        my $length = $seq_level_coord->end - $seq_level_coord->start + 1;
        my @line = ($slice->seq_region_name);
        
        if ($seq_level_coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          my $seq_level_name = $ama->seq_ids_to_regions([$seq_level_coord->id]);
          my $orientation = ($seq_level_coord->strand eq -1) ? '-' : '+';
          
          push @line, (
            $asm_start,
            $asm_start + $length - 1,
            ++$cmp_count,
            'W',
            $$seq_level_name[0],
            $seq_level_coord->start,
            $seq_level_coord->end,
            $orientation
          );
          
        } elsif ($seq_level_coord->isa('Bio::EnsEMBL::Mapper::Gap')) {          
          push @line, (
            $seq_level_coord->start,
            $seq_level_coord->end,
            ++$cmp_count,
            'N',
            $length,
            $gap_type,
            $linkage,
            $evidence
          );
          
        }
        print $out_fh join("\t", @line)."\n";
        
        $asm_start += $length;
      }
      close($out_fh);
    }
  }
  
  $self->param("agp_files", \%agp_files);
}

sub write_output {
  my ($self) = @_;
  $self->dataflow_output_id(
    {
      'agp_files' => $self->param('agp_files'),
    }, 2);
}

sub get_coord_maps {
  my ($self, $dba) = @_;

  my $pairs_sql = '
    SELECT sa.coord_system_id, sc.coord_system_id
    FROM assembly a
    LEFT JOIN seq_region sa ON a.asm_seq_region_id = sa.seq_region_id
    LEFT JOIN seq_region sc ON a.cmp_seq_region_id = sc.seq_region_id
    LEFT JOIN coord_system ca ON sa.coord_system_id = ca.coord_system_id
    LEFT JOIN coord_system cc ON sc.coord_system_id = cc.coord_system_id
    WHERE ca.attrib LIKE "%default_version%"
      AND cc.attrib LIKE "%default_version%"
    GROUP BY sa.coord_system_id, sc.coord_system_id;
  ';
  
  my $dbh = $dba->dbc->db_handle();

  my $sth = $dbh->prepare($pairs_sql);
  $sth->execute();

  my @pairs;
  while (my @pair = $sth->fetchrow_array()) {
    push @pairs, \@pair;
  }

  return @pairs;
}

1;
