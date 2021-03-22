import sys
import os.path
import logging
import gzip

from mouse.gffEntry import GffEntry

################################################################################
#

################################################################################
#
# class GFF -- holds the contents of the GFF file (that we're interested in,
#              anyway).

class GffFile:

  TRANSCRIPTS = set(["mRNA","miRNA","ncRNA","pre_miRNA","pseudogene",
                     "rRNA","snRNA","snoRNA","tRNA","C_gene_segment",
                     "J_gene_segment","D_gene_segment","V_gene_segment",
                     "gene_segment"])

  TRANSCRIPT_PARTS = set(["exon","five_prime_UTR","three_prime_UTR","CDS",
                          "intron"])

  TOP_LEVEL = set(["gene","transposable_element","miRNA","miRNA_5p",
                   "miRNA_3p","miRNA_precursor"])

  def __init__(self):
    self.entries = {}
    self.parents = {}
    self.offspring = {}

  @classmethod
  def open(self,fn):
    split = os.path.splitext(fn)
    if split[1] == ".gz":
      fd = gzip.open(fn,"rt")
    else:
      fd = open(fn)
    return fd

  def load(self,fn,fd):
    for line in fd:
      if line == "##FASTA\n":
        break
      elif line[0] == "#":
        continue
      else:
        entry = GffEntry.parse(line)
        if entry == None:
          continue
        if entry.oid in self.entries and entry.category != "transposable_element":
          logging.getLogger().warning("SKIPPING duplicate id: %s" % (entry.oid,))
        else:
          self.entries[entry.oid] = entry
          self.parents[entry.oid] = entry.parent
          if entry.parent != None and len(entry.parent) > 0:
            for par in entry.parent:
              if par not in self.offspring:
                self.offspring[par] = set()
              self.offspring[par].add(entry.oid)
          elif entry.category == 'mRNA' and entry.oid not in self.parents:
            sys.stderr.write("mRNA missing parents: %s" % (entry.oid))
          elif entry.category == 'mRNA' and len(self.parents[oid]) == 0:
            sys.stderr.write("mRNA zero parents: %s" % (entry.oid))
    fd.close()

#  def parse(self,line):
#    flds = line.strip().split('\t')
#    if flds[2] in ('gene','exon','mRNA'):
#      self.log.debug("Entry: '%s' '%s'" % (flds[2],line.strip()))
#    if flds[2] in self.exclude:
#      return None
#    self.log.debug("zork! %s" % (flds[2],))
#    aList = flds[8].split(';')
#    aDict = {}
#    for a in aList:
#      aflds = a.split('=')
#      if aflds[0] == "Parent":
#        aDict[aflds[0]] = set(aflds[1].split(','))
#      else:
#        aDict[aflds[0]] = aflds[1]
#    x = Entry(flds[0],flds[2],flds[3],flds[4],flds[6],aDict)
#    return Entry(flds[0],flds[2],flds[3],flds[4],flds[6],aDict)
#    self.log.debug("class x: %s" % (type(x),))
#    return x

  def fixName(self,name):
    x = name.find('{')
    if x >= 0:
      name = name[0:x]
    if name[0:6] == "snRNA:":
      name = name[6:]
    elif name[0:7] == "snoRNA:":
      name = name[7:]
    elif name[0:5] == "tRNA:":
      name = name[5:]
    return name

  def makeTrackSet(self,base,ucsc=False):
    fds = {}
    for tag in trackSet:
      fds[tag] = open("%s_%s.bed" % (base,tag),'w')
    for entry in self.entries.values():
      if entry.type in ("miRNA","transposable_element"):
        print(entry.bedWithName(self.fixName(entry.name())),file=fds[entry.type])
      elif entry.type == "gene":
        kids = self.offspring[entry.fid()]
        for k in kids:
          thing = BedStruct(k,self.fixName(k.name()),self.offspring.get(k.fid(),[]))
          print(thing.bed(),file=fds[k.type])
    for fd in fds.values():
      fd.close()

  def makePipelineSet(self,base,ucsc=False):
    fds = {}
    for tag in trackSet:
      fds[tag] = open("%s_%s.bed" % (base,tag),'w')
    for entry in self.entries.values():
      if entry.type in ("miRNA","transposable_element"):
        print(entry.bed(ucsc=ucsc),file=fds[entry.type])
      elif entry.type == "gene":
        if entry.fid() not in self.offspring:
          continue
        kids = self.offspring[entry.fid()]
        gtype = kids[0].type
        print(entry.bed(ucsc=ucsc),file=fds[gtype])
    for fd in fds.values():
      fd.close()

  def topLevelTypes(self,fd):
    noParent = set()
    for e in self.entries.values():
      if e.parent == None or len(e.parent) == 0:
        noParent.add(e.category)
    for fid in noParent:
      print(fid,file=fd)

  def parentTypes(self,fd):
    pt = {}
    for e in self.entries.values():
      if e.parent != None and len(e.parent) > 0:
        etype = e.category
        for p in e.parent:
          if etype not in pt:
            pt[etype] = set()
          if p not in self.entries:
            pt[etype].add("orphan")
            sys.stderr.write("Orphan: %s %s %s\n" % (e.category,e.oid,p))
          else:
            pt[etype].add(self.entries[p].category)
    for e,p in iter(pt.items()):
      fd.write("%s : %s\n" % (e,",".join(p)))

  def hasName(self,fd):
    cats = set()
    for e in self.entries.values():
      if e.name() != None:
        cats.add(e.type)
    catsL = list(cats)
    catsL.sort()
    for c in catsL:
      print(c,file=fd)

  def hasFullName(self,fd):
    cats = set()
    for e in self.entries.values():
      if e.fullname() != None:
        cats.add(e.type)
    catsL = list(cats)
    catsL.sort()
    for c in catsL:
      print(c,file=fd)

  def stats(self,fd):
    tset = {}
    for e in self.entries.values():
      t = e.type
      if t not in tset:
        tset[t] = 0
      tset[t] += 1
    for t,v in iter(tset.items()):
      fd.write("%09d\t%s\n" % (v,t))

  def checkConsistent(self,fd):
    # first check whether each gene has all kids the same type
    for parent,kids in iter(self.offspring.items()):
      p_entry = self.entries[parent]
      if p_entry.type == "gene":
        ktypes = set()
        for k in kids:
          ktypes.add(k.type)
        if len(ktypes) > 1:
          fd.write("%s: kid types: %s\n",parent,",".join(ktypes))
        elif len(ktypes) == 0:
          fd.write("%s: gene with no offspring\n",parent)

  def maxKids(self,fd):
    # report the maximum number of children of each type a parent type has
    types = {}
    for parent,kids in iter(self.offspring.items()):
      p_entry = self.entries[parent]
      kmap = {}
      for k in kids:
        kmap[k.type] = kmap.get(k.type,0) + 1
      if p_entry.type not in types:
        types[p_entry.type] = {}
      gset = types[p_entry.type]
      for k,v in iter(kmap.items()):
        cmax = gset.get(k,0) 
        if v > cmax:
          gset[k] = v
        if v == 82:
          print(p_entry.bed())
    for t,tmap in iter(types.items()):
      print(t,file=fd)
      for k,v in tmap.items():
        print("    %s : %d" % (k,v),file=fd)

  def makeAnno(self,category,ucsc,outFD):
    sys.stderr.write("cat=%s ucsc=%s\n" % (category,ucsc))
    if category in GffFile.TOP_LEVEL:
      for (key,e) in self.entries.items():
        if e.category == category:
          attr = "gene_id \"%s\"; transcript_id \"\"; gene_name \"%s\";" % (key,e.name())
          outFD.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (e.chrom,e.source,e.category,e.left,e.right,e.score,e.strand,e.frame,attr))
    elif category in GffFile.TRANSCRIPTS:
      written = set() # for exons we've already written
      for (key,e) in self.entries.items():
        if e.category == category:
          if e.source == "RNAcentral" and e.parent == None:
            continue
          if len(e.parent) > 1:
            sys.stderr.write("PARENT: transcript %s has parents '%s'\n" % (key,",".join(e.parent)))
          gene_id = next(iter(e.parent)) # should only be one...
          gene_name = self.entries[gene_id].name()
          transcript_id = e.name()
          if key in self.offspring:
            for exonId in self.offspring[key]:
              if exonId not in written:
                exon = self.entries[exonId]
                if exon.category == "exon":
                  attr = "gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";" % (gene_id,transcript_id,gene_name)
                  outFD.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (exon.chrom,exon.source,exon.category,exon.left,exon.right,exon.score,exon.strand,exon.frame,attr))
                  written.add(exonId)
          else:
            attr = "gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";" % (gene_id,transcript_id,gene_name)
            outFD.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (e.chrom,e.source,e.category,e.left,e.right,e.score,e.strand,e.frame,attr))
    elif category in GffFile.TRANSCRIPT_PARTS:
      for (key,e) in self.entries.items():
        if e.category == category:
          grandparents = set()
          for par in e.parent:
            obj = self.entries[par]
            if obj.parent != self.parents[par]:
              sys.stderr.write("Parental mismatch: %s: obj.parent %s != parents[] %s" % (par,obj.parent,self.parents[par]))
            if obj.parent != None and len(obj.parent) > 0:
              grandparent = next(iter(self.entries[par].parent)) # we know there's only one parent per transcript
              grandparents.add(grandparent)
            elif obj.parent == None:
              sys.stderr.write("Missing grandparent: %s %s %s parents='%s'\n" % (obj.category,par,obj.name(),",".join(obj.parent)))
            elif len(obj.parent) == 0:
              sys.stderr.write("Zero grandparent: %s %s %s parents='%s'\n" % (obj.category,par,obj.name(),",".join(obj.parent)))
#              grandparents.add(par)
          if len(grandparents) > 1:
            sys.stderr.write("GRANDPARENTS: %s %s has multiple grandparents '%s', choosing one arbitrarily\n" % (category,key,",".join(grandparents)))
            for par in e.parent:
              parobj = self.entries[par]
              sys.stderr.write("    Parent: %s %s %s\n" % (parobj.category,par,parobj.name()))
          if len(grandparents) > 0:
            gene_id = next(iter(grandparents))
            grandparent = self.entries[gene_id]
            transcript_id = next(iter(e.parent))
            parent = self.entries[transcript_id]
            attr = "gene_id \"%s\"; transcript_id \"%s\"; gene_name \"%s\";" % (gene_id,transcript_id,grandparent.name())
            outFD.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (e.chrom,e.source,e.category,e.left,e.right,e.score,e.strand,e.frame,attr))
