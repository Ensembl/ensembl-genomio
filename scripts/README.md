* `get_dna_fasta_for.pl` -- Get toplevel dna sequences for the list of ids in stdin
* `load_fann.pl` -- Xref linked load of functional annotation.
* `gff_stats.py` -- check GFF validity based on config/gff_metaparser/valid_structures.conf, log, try to fix
  ls data/*/*gff*gz | xargs -n 1 python3 gff_stats.py 2>&1 | tee gff3_seen_cases.log

