package Bio::EnsEMBL::Pipeline::Runnable::DumpSeqRegionJson;

use strict;
use warnings;

use JSON;
use File::Path qw(make_path);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Production::Pipeline::Common::Base');

sub fetch_input {
  my ($self) = @_;

  $self->param( 'dba', $self->core_dba() );

  return;
}

sub param_defaults {
  return {};
}

sub run {
  my ($self) = @_;

  print $self->param('species');
  $self->write_json();
  return;
}

sub write_json {
  my ($self) = @_;

  my $sub_dir = $self->create_dir('json');
  $self->info(
          "Processing " . $self->production_name() . " into $sub_dir" );
  my $dba = $self->core_dba();

  # Get seq_regions
  my $sa = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "slice" );
  my $syna = Bio::EnsEMBL::Registry->get_adaptor(
    $self->production_name, "core", "seqregionsynonym" );

  # Get all top level seq regions
  my $slices = $sa->fetch_all('toplevel');

  my @seq_regions;
  foreach my $slice (@$slices) {
    my $syns = $syna->get_synonyms( $slice->get_seq_region_id() );
    my $seq_region = {
      name => $slice->seq_region_name(),
      coord_system_level => $slice->coord_system_name(),
      length => $slice->length(),
    };

    # Additional metadata
    $seq_region->{synonyms} = [ map { $_->name } @$syns ] if @$syns;
    $seq_region->{circular} //= $slice->is_circular;

    push @seq_regions, $seq_region;
  }

  $dba->dbc()->disconnect_if_idle();

  # Write data to json
  $self->write_json_file($sub_dir, \@seq_regions);
} ## end sub write_json

sub write_json_file {
  my ($self, $sub_dir, $data) = @_;
  my $json_file_path =
    $sub_dir . '/' . $self->production_name() . '.json';
  $self->info("Writing to $json_file_path");
  open my $json_file, '>', $json_file_path or
    die "Could not open $json_file_path for writing";
  print $json_file encode_json($data);
  close $json_file;
  $self->info("Write complete");
  return;
}

1;
