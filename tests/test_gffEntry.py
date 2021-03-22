import sys
import unittest
import io
import logging

sys.path.insert(0,"..")

from flybase.gffEntry import GffEntry

class TestGffEntry(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    # log to string if possible
    #TestGffEntry.logBuffer = io.StringIO()
    #TestGffEntry.hdlr = logging.StreamHandler(stream=TestGffEntry.logBuffer)
    #logging.getLogger().addHandler(TestGffEntry.hdlr)
    pass

  def setUp(self):
    # log to string if possible
    self.logBuffer = io.StringIO()
    self.hdlr = logging.StreamHandler(stream=self.logBuffer)
    logging.getLogger().addHandler(self.hdlr)

  @classmethod
  def tearDownClass(cls):
    #logging.getLogger().removeHandler(TestGffEntry.hdlr)
    pass

  def tearDown(self):
    logging.getLogger().removeHandler(self.hdlr)

  def test_sanity(self):
    self.assertEqual(True,True)

  def test_parseAttrSimple(self):
    attrs = "ID=FBtr0300689;Name=CG11023-RB;Parent=FBgn0031208"
    attrMap = GffEntry.parseAttributes(attrs)
    self.assertEqual("FBtr0300689",attrMap["ID"])
    self.assertEqual("CG11023-RB",attrMap["Name"])
    self.assertEqual(set(["FBgn0031208"]),attrMap["Parent"])

  def test_parseAttrList(self):
    attrs = "ID=FBgn0031208_d7236e17720;Name=CG11023;to_species=Danio rerio;to_name=Drer\si:dkey-21c19.3;Dbxref=ZFIN:ZDB-GENE-030131-6462,EntrezGene:100007636,Ensembl_Danio_rerio:ENSDARG00000078574;diopt_source=eggNOG,Inparanoid,OrthoDB,Phylome,RoundUp,TreeFam"
    attrMap = GffEntry.parseAttributes(attrs)
    self.assertTrue("RoundUp" in attrMap["diopt_source"])
    self.assertTrue("eggNOG" in attrMap["diopt_source"])
    self.assertTrue("Inparanoid" in attrMap["diopt_source"])
    self.assertTrue("OrthoDB" in attrMap["diopt_source"])
    self.assertTrue("Phylome" in attrMap["diopt_source"])
    self.assertTrue("RoundUp" in attrMap["diopt_source"])
    self.assertTrue("TreeFam" in attrMap["diopt_source"])

  def test_parseSimple(self):
    line = "2L	FlyBase	CDS	7680	8116	.	+	0	Parent=FBtr0300689,FBtr0300690\n"
    entry = GffEntry.parse(line)
    self.assertTrue("FBtr0300689" in entry.parent)
    self.assertTrue("FBtr0300690" in entry.parent)
    self.assertEqual(entry.chrom,"chr2L")
    self.assertEqual(entry.source,"FlyBase")
    self.assertEqual(entry.left,7680)
    self.assertEqual(entry.oid,id(entry))

  def test_parseUnknown(self):
    line = "2L	FlyBase	Zork	7680	8116	.	+	0	Parent=FBtr0300689,FBtr0300690\n"
    entry = GffEntry.parse(line)
    self.assertEqual(entry,None)
    self.assertEqual("UNKNOWN entryType 'Zork'\n",self.logBuffer.getvalue())
