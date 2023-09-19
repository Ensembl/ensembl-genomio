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


=pod

=head1 NAME

Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AnalysisSetup

=head1 DESCRIPTION

Add a new analysis to a database.

If the analysis exists already, the default behaviour is to rename the
existing analysis, and insert a fresh analysis. This is useful if you want
to easily compare new and old results associated with that analysis, but means
that you are responsible for subsequently tidying up the database yourself,
since you probably don't want both new and old analyses in a production db.

Alternatively, if the analysis exists you can delete it, and any rows
connected to it (delete_existing=1). You are responsible for providing
an arrayref of the tables which have the relevant analysis_id foreign
keys (linked_tables).

The script will look up descriptions from the production database, by default;
this relies on the production server (usually pan-prod) being in the registry.
Lookup parameters will be over-ridden if description parameters are explicitly
given, or can be set as undefined (production_lookup=0).

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::AnalysisSetup;

use strict;
use warnings;

#use base ('Bio::EnsEMBL::Hive::Process');
use base ('Bio::EnsEMBL::Pipeline::Runnable::EG::LoadGFF3::Base');
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Hive::Utils::URL qw/ parse /;

sub param_defaults {
  return {
    # basic
    logic_name         => undef,
    module             => undef,
    production_lookup  => 1,
    # aux
    delete_existing    => 0,
    logic_rename       => undef,
    keep_logic_name   => 0,
    linked_tables      => [],
  };
}

sub fetch_input {
  my $self = shift @_;
  
  my $logic_name = $self->param_required('logic_name');
  $self->param('program', $logic_name) unless $self->param_is_defined('program');
  
  if (not $self->param('delete_existing') and not $self->param('keep_logic_name')) {
    $self->param('logic_rename', "$logic_name\_bkp") unless $self->param_is_defined('logic_rename');
  }
}

sub run {
  my $self = shift @_;
  my $logic_name = $self->param_required('logic_name');
  
  my $dba = $self->get_type_dba();
  my $dbh = $dba->dbc->db_handle;
  my $aa = $dba->get_adaptor('Analysis');
  my $analysis = $aa->fetch_by_logic_name($logic_name);
  
  if (defined $analysis) {
    if ($self->param('delete_existing')) {
      my $analysis_id = $analysis->dbID;
      foreach my $table (@{$self->param('linked_tables')}) {
        my $sql = "DELETE FROM $table WHERE analysis_id = $analysis_id";
        my $sth = $dbh->prepare($sql) or throw("Failed to delete rows using '$sql': ".$dbh->errstr);
        $sth->execute or throw("Failed to delete rows using '$sql': ".$sth->errstr);
      }
      $aa->remove($analysis);
    } elsif (not $self->param('keep_logic_name')) {
      my $logic_rename = $self->param_required('logic_rename');
      my $renamed_analysis = $aa->fetch_by_logic_name($logic_rename);
      if (defined $renamed_analysis) {
        $self->throw(
          "Cannot rename '$logic_name' to '$logic_rename' because '$logic_rename' already exists.\n".
          "Either provide a different 'logic_rename' parameter, or delete the '$logic_rename' analysis.\n"
        );
      } else {
        $analysis->logic_name($logic_rename);
        $aa->update($analysis);
      }
    }
  }
  
  if ($self->param('production_lookup')) {
    $self->production_updates;
  }
  
  my $new_analysis = $self->create_analysis;
  $aa->store($new_analysis);
 
  $dba->dbc->disconnect_if_idle();
}

sub create_analysis {
  my ($self) = @_;
  
  my $analysis = Bio::EnsEMBL::Analysis->new(
    -logic_name      => $self->param('logic_name'),
    -db              => $self->param('db'),
    -db_version      => $self->param('db_version'),
    -db_file         => $self->param('db_file'),
    -program         => $self->param('program'),
    -program_version => $self->param('program_version'),
    -program_file    => $self->param('program_file'),
    -parameters      => $self->param('parameters'),
    -module          => $self->param('module'),
    -module_version  => $self->param('module_version'),
    -gff_source      => $self->param('gff_source'),
    -gff_feature     => $self->param('gff_feature'),
    -description     => $self->param('description'),
    -display_label   => $self->param('display_label'),
    -displayable     => $self->param('displayable'),
    -web_data        => $self->param('web_data'),
  );
  
  return $analysis;
}

sub production_updates {
  my ($self) = @_;
  my $logic_name = $self->param('logic_name');
  my $dba        = $self->production_dba();
  my $dbc        = $dba->dbc;
  my $dbh        = $dbc->db_handle();
  my %properties;
  
  # Load generic, non-species-specific, analyses
  my $sth = $dbh->prepare(
    'SELECT ad.description, ad.display_label, 1, wd.data '.
    'FROM analysis_description ad '.
     'LEFT OUTER JOIN web_data wd ON ad.web_data_id = wd.web_data_id '.
    'WHERE ad.is_current = 1 '.
    'AND ad.logic_name = ? '
  );
  $sth->execute($logic_name);
  
  $sth->bind_columns(\(
    $properties{'description'},
    $properties{'display_label'},
    $properties{'displayable'},
    $properties{'web_data'},
  ));
  $sth->fetch();
  
  $properties{'web_data'} = eval { $properties{'web_data'} } if defined $properties{'web_data'};
  
  # Explicitly passed parameters do not get overwritten.
  foreach my $property (keys %properties) {
    if (! $self->param_is_defined($property)) {
      $self->param($property, $properties{$property});
    }      
  }

  $dbc->disconnect_if_idle();
}

1;
