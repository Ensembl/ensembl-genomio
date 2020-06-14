import sys

class FixAction:
  def __init__(self, raw, rule):
     self._raw = raw
     self._rule_name = rule.NAME
     self._rule_lineno = rule._lineno
     self._action = self.parse(self._raw)

  def __bool__(self):
    return self._action is not None

  def __repr__(self):
    if self._action is None:
      return str(None)
    return "action: %s of %s at %s" % (
      self._raw, self._rule_name, self._rule_lineno
    )

  def parse(self, raw):
    return raw.split("/")
    #if not self._action: self._action = None
    #if not self._action[0]: self._action = None

    #gene/@MRNA/@CDS	SUB	gene/transcript.biotype=@MRNA/@CDS
    #
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/pre_miRNA/-/exon
    #gene/primary_transcript/mirna/exon	SUB	ncRNA_gene/-/miRNA/exon
    #mirna	SUB	+ncRNA_gene/miRNA.biotype=miRNA/+exon
    #
    # balanced number of all - added == number of old


