.. See the NOTICE file distributed with this work for additional information
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

.. ensembl-genomio documentation master file, created by
   sphinx-quickstart on Tue May  3 16:24:49 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===========================================
Ensembl-genomio
===========================================

A repository dedicated to pipelines used to turn basic genomic data into formatted 
Ensembl core databases. Also allow users to dump core databases into various formats.

File formats handled : FastA, gff3, JSON (*following BRC4 specifications*).

Check out the :doc:`usage` section for further information of requirements to
run ensembl-genomio pipelines.

Ehive pipelines
-------------------------------------------
#. 'Genome loader': Creates an Ensembl core database from a set of flat files.
#. 'Genome dumper': Dumps flat files from an Ensembl core database.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Contents
--------
Check out :ref:`installation <install>` section for further information on how 
to install the project.

.. toctree::
   :maxdepth: 1
   :caption: Index

   usage
   install
   modules
   license
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
