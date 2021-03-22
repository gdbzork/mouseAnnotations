import sys
import os.path
import argparse
import logging
import re

################################################################################

class GffEntry:

  EXCLUDE = set(["BAC_cloned_genomic_insert","DNA_motif","RNAi_reagent",
                 "TF_binding_site","breakpoint","chromosome","chromosome_band",
                 "complex_substitution","deletion","enhancer","exon_junction",
                 "golden_path_region","insertion_site",
                 "insulator","match","match_part","mature_peptide",
                 "modified_RNA_base_feature","oligonucleotide",
                 "origin_of_replication","orthologous_region","orthologous_to",
                 "pcr_product","point_mutation","polyA_site",
                 "protein","protein_binding_site","region","regulatory_region",
                 "repeat_region","rescue_fragment","sequence_variant",
                 "silencer","syntenic_region","tandem_repeat",
                 "transcription_start_site",
                 "transposable_element_insertion_site",
                 "uncharacterized_change_in_nucleotide_sequence",
                 "biological_region","supercontig","start_codon",
                 "stop_codon"])

  # we have INCLUDE as well as EXCLUDE so we can detect categories that
  # we didn't know were there at all
  INCLUDE = set(["CDS","exon","five_prime_UTR","gene","mRNA","miRNA","ncRNA",
                 "pre_miRNA","pseudogene","rRNA","snRNA","snoRNA",
                 "three_prime_UTR","tRNA","transposable_element","lnc_RNA",
                 "transcript","ncRNA_gene","pseudogenic_transcript","scRNA",
                 "intron","C_gene_segment","J_gene_segment","D_gene_segment",
                 "V_gene_segment","gene_segment","miRNA_5p","miRNA_3p",
                 "miRNA_precursor","noncoding_exon"])

  ATTR_KEY_PAT = re.compile("^(\w+)$")

  # includes human and mouse conversions
  CHROM2UCSC = {"2L": "chr2L",
                "2R": "chr2R",
                "3L": "chr3L",
                "3R": "chr3R",
                "1": "chr1",
                "2": "chr2",
                "3": "chr3",
                "4": "chr4",
                "5": "chr5",
                "6": "chr6",
                "7": "chr7",
                "8": "chr8",
                "9": "chr9",
                "10": "chr10",
                "11": "chr11",
                "12": "chr12",
                "13": "chr13",
                "14": "chr14",
                "15": "chr15",
                "16": "chr16",
                "17": "chr17",
                "18": "chr18",
                "19": "chr19",
                "X": "chrX",
                "Y": "chrY",
                "MT": "chrM",
                "mitochondrion_genome": "chrM"}

################################################################################

  @staticmethod
  def parseRNAcentral(attrs):
    attrList = attrs.split(";")
    attrMap = {}
    for attr in attrList:
      (k,v) = attr.split(" ")
      v = v.replace("\"","")
      if k == "Parent":
        v = set([v])
      attrMap[k] = v
    return attrMap
    
  @staticmethod
  def parseAttributes(attrs):
    attrList = attrs.split(";")
    attrMap = {}
    for attr in attrList:
      (key,val) = attr.strip().split("=")
      valList = val.split(",")
#      valList = list(map(lambda x:x.split(":")[-1],valList))
      if len(valList) > 1 or key == "Parent":
        attrMap[key] = set(valList)
      else:
        attrMap[key] = valList[0]
    return attrMap

  @staticmethod
  def parse(line,ucsc=True):
    flds = line.strip().split("\t")
    entryType = flds[2]
    if entryType in GffEntry.EXCLUDE:
      return None
    elif entryType not in GffEntry.INCLUDE:
      logging.getLogger().warning("UNKNOWN entryType '%s'" % (entryType,))
      return None
    if flds[1] == "RNAcentral":
      attrMap = GffEntry.parseRNAcentral(flds[8])
    else:
      attrMap = GffEntry.parseAttributes(flds[8])

    if flds[2] == 'mRNA' and attrMap['Parent'] == None:
      sys.stderr.write("Line '%s': parent missing from mRNA\n" % (line.strip(),))
    if flds[2] == 'mRNA' and (attrMap['Parent'] != None and len(attrMap['Parent']) == 0):
      sys.stderr.write("Line '%s': parent set empty\n" % (line.strip(),))
    chrom = flds[0]
    if ucsc:
      chrom = GffEntry.CHROM2UCSC[chrom] if chrom in GffEntry.CHROM2UCSC else chrom
    if flds[1] == "RNAcentral":
      entry = GffEntry(chrom,flds[1],attrMap['type'],flds[3],flds[4],flds[5],flds[6],flds[7],attrMap)
    else:
      entry = GffEntry(chrom,flds[1],flds[2],flds[3],flds[4],flds[5],flds[6],flds[7],attrMap)
    return entry

################################################################################

  def __init__(self,chrom,source,fType,left,right,score,strand,frame,attrDict):
    self.chrom = chrom
    self.source = source
    self.category = fType
    self.left = int(left)
    self.right = int(right)
    self.score = score
    self.strand = strand
    self.frame = frame
    self.attrs = attrDict
    self.oid = self.attrs["ID"] if "ID" in self.attrs else id(self)
    if self.oid == id(self) and fType not in ("exon","CDS","three_prime_UTR","five_prime_UTR","intron"):
      logging.getLogger().warning("ENTRY missing ID: %s:%d-%d %s %s" % (self.chrom,self.left,self.right,self.source,self.category))
    self.parent = self.attrs["Parent"] if "Parent" in self.attrs else None

  def name(self):
    return self.attrs.get("Name",None)

  def fullname(self):
    return self.attrs.get("fullname",None)

  def attribute(self,attr):
    return self.attrs.get(attr,"None")

  def attrNames(self):
    return self.attrs.keys()
