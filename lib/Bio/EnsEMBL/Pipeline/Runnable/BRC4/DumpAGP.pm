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

sub param_defaults {
  my ($self) = @_;
  
  return {
    %{$self->SUPER::param_defaults},
    'seq_level' => 'contig',
  };
}

sub fetch_input {
  my ($self) = @_;
  
  my $sub_dir = $self->create_dir('agp');
  my $agp_file = catfile($sub_dir, $self->production_name() . '.agp');
  $self->param('agp_file', $agp_file);
}

sub run {
  my ($self) = @_;
  my $species      = $self->param_required('species');
  my $agp_file     = $self->param_required('agp_file');
  my $seq_level    = $self->param_required('seq_level');
  
  my $gap_type = 'scaffold';
  my $linkage  = 'yes';
  my $evidence = 'paired-ends';
  
  open(my $out_fh, '>', $agp_file) or $self->throw("Cannot open file $agp_file: $!");
  
  my $dba = $self->core_dba();
  my $ama = $dba->get_adaptor('AssemblyMapper');
  my $csa = $dba->get_adaptor('CoordSystem');
  my $sa  = $dba->get_adaptor('Slice');
  
  my $seq_level_cs = $csa->fetch_by_name($seq_level);  
  my @top_level_names = $self->top_level_names($dba);
  
  foreach my $top_level_name (sort {$a cmp $b} @top_level_names) {
    my $top_level_cs = $csa->fetch_by_name($top_level_name);
    my $mapper = $ama->fetch_by_CoordSystems($top_level_cs, $seq_level_cs);
    my $slices = $sa->fetch_all($top_level_name);
    
    foreach my $slice (sort {$a->seq_region_name cmp $b->seq_region_name} @$slices) {
      my @seq_level_coords =
        $mapper->map(
          $slice->seq_region_name,
          $slice->start,
          $slice->end,
          $slice->strand,
          $top_level_cs
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
    }
  }
  
  close($out_fh);
}

sub write_output {
  my ($self) = @_;
  $self->dataflow_output_id({ 'agp_file' => $self->param('agp_file')}, 2);
}

sub top_level_names {
  my ($self, $dba) = @_;
  
  my $top_level_sql = '
    SELECT cs.name FROM 
      coord_system cs inner join 
      seq_region sr using (coord_system_id) 
      inner join seq_region_attrib sra using (seq_region_id) inner join 
      attrib_type at using (attrib_type_id) 
    where at.code="toplevel" 
    group by cs.name;
  ';
  
  my $dbh = $dba->dbc->db_handle();
  my $top_level_names = $dbh->selectcol_arrayref($top_level_sql);
  
  return @$top_level_names;
}

1;
