import unittest
import sys
import io
import logging

sys.path.insert(0,"..")

class TestUtils(unittest.TestCase):

  def logSetup(self):
    self.logstream = io.StringIO("")
    self.loghdlr = logging.StreamHandler(self.logstream)
    self.loghdlr.setLevel(logging.DEBUG)
    rootlog = logging.getLogger()
    try:
      self.stdhdlr = rootlog.handlers[0]
    except:
      self.stdhdlr = None
    rootlog.addHandler(self.loghdlr)
    self.oldLogLev = rootlog.getEffectiveLevel()
    rootlog.setLevel(logging.DEBUG)
    if self.stdhdlr != None:
      rootlog.removeHandler(self.stdhdlr)

#  def getLogs(self):
#    return self.logstream

  def logReset(self):
    rootlog = logging.getLogger()
    if self.stdhdlr != None:
      rootlog.addHandler(self.stdhdlr)
    rootlog.removeHandler(self.loghdlr)
    rootlog.setLevel(self.oldLogLev)

