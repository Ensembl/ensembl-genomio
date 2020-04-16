import gzip
import re

from collections import defaultdict
from Bio import SeqIO


class IdTrimmer:
  def __init__(self, re_rules = dict()):
    self._rules = None
    self.compile_rules(re_rules)

  def compile_rules(self, re_rules):
    if re_rules:
      self._rules = { k.lower():re.compile(v) for k,v in re_rules.items() }

  def normalize(self, id_str, type = None):
    id_str = str(id_str)
    _type = type.lower()
    if not self._rules:
      return id_str
    if type is None:
      _type = "any"
    if _type in self._rules:
      return self._rules[_type].sub("", id_str)
    if "any" in self._rules:
      return self._rules["any"].sub("", id_str)
    return id_str

  def __call__(self, *args, **kwargs):
    return self.normalize(*args, **kwargs)


class PfxTrimmer(IdTrimmer):
  def __init__(self, trim_str):
    self._rules = None
    if not trim_str:
      return
    rules = defaultdict(list)
    for trim_pair in filter(None,trim_str.split(",")):
      tag, *pat = trim_pair.split(":",1)
      if not pat:
        tag, pat = "ANY", [tag]
      rules[tag].append(re.escape(pat[0]))
    rules = { k: r"^(?:%s)" % ("|".join(v)) for k, v in rules.items() }
    self.compile_rules(rules)


class ExtMapper:
  def __init__(self, tag,  map_file = None, map_str = None):
    pass

  def map(self, val):
    pass

  def __call__(self, *args, **kwargs):
    return self.map(*args, **kwargs)

class SeqLenDict:
  def __init__(self, fna_file = None):
    self._len = None
    self.load_from_file(fna_file)

  def load_from_file(self, fasta_file):
    if not fasta_file:
      return
    _open = fasta_file.endswith(".gz") and gzip.open or open
    with _open(fasta_file, 'rt') as fasta:
      fasta_parser = SeqIO.parse(fasta, "fasta")
      self._len = dict()
      for rec in fasta_parser:
        self._len[rec.name] = len(rec)

  def get_len(self, srid):
    if self._len and srid in self._len:
      return self._len[srid]
    return None

  def __call__(self, srid):
    return self.get_len(srid)


class UpdatingLen:
  def __init__(self, val):
    self._update = False
    self._val = val
    if not self._val:
      self._update = True

  def is_updating(self):
    return self._update

  def update(self, val, stop_on_success = False):
    if not self._update:
      return
    self._val = val
    if self._val and stop_on_success:
      self._update = False

  def __call__(self, *args, **kwargs):
    return self.update(*args, **kwargs)

  def __int__(self):
    return self._val or int(0)

  def __len__(self):
    return self.__int__()


