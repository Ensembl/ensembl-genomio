# finding most commomn prefix
class PfxTr:
  def __init__(self):
    self.letters = dict()
    self.count = int(0)

  def add(self, s):
    if s is None or s == "":
      return
    root = self
    for c in [""] + list(s):
      if c not in root.letters:
        root.letters[c] = PfxTr()
      root = root.letters[c]
      root.count += 1

  def get_max(self):
     (pfx, max_cnt) = ("", 0)
     root = self
     while root.letters:
       c = max(root.letters.keys(), key = lambda x: root.letters[x].count)
       if root.letters[c].count < max_cnt:
         return(pfx, max_cnt)
       max_cnt = root.letters[c].count
       pfx += c
       root = root.letters[c] 
     return (pfx, max_cnt)

