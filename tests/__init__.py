import unittest
import numpy as np
from os import path
import bioitools
from bioitools import smoothers
from bioitools import correlations

class TestBioitools(unittest.TestCase):
	def setUp(self):
		self.faiFile = path.join(path.dirname(__file__),'test.fa.fai')
		self.goodBed = path.join(path.dirname(__file__),'test.bedgraph')
		self.badBed = path.join(path.dirname(__file__),'bad.bedgraph')
	def test_parseFai(self):
		self.assertEqual(bioitools.ParseFai(self.faiFile), {'Chr1':30427671, 'Chr2':19698289, 'Chr3':23459830})
	def test_bounds(self):
		self.assertEqual(bioitools._ChromBounds([1,1,1,2,2,3,4]), {1:(0,3), 2:(3,5), 3:(5,6), 4:(6,7)})
	def test_contiguous(self):
		self.assertEqual(bioitools.smoothers.isContiguous(np.array([0,10,20]), np.array([10,20,30]), {1:(0,3)}), True)
		self.assertEqual(smoothers.isContiguous(np.array([0,10,20]), np.array([10,30,40]), {1:(0,3)}), False)
	def test_reader(self):
		c,s,e,v = bioitools.ParseBedgraph(self.goodBed)
		np.testing.assert_array_equal(v,np.array([1,2,3,1,2,3,4,0]))
		np.testing.assert_allclose(smoothers.haar(v),np.array([2,2,2,2,2,2,4,0]), atol=0.00001)
		cDict = bioitools._ChromBounds(c)
		self.assertEqual(cDict,{'1':(0,8)})
		self.assertEqual(smoothers.isContiguous(s,e,cDict),True)
	def test_badFile(self):
		c,s,e,v = bioitools.ParseBedgraph(self.badBed)
		chromDict = bioitools._ChromBounds(c)
		self.assertEqual(smoothers.isContiguous(s,e,chromDict),False)
	def test_phiCorr(self):
		A = np.array([0,0,1,1,0,0],dtype=np.bool)
		self.assertEqual(correlations.phiCorr(A,A), 1)
		self.assertEqual(correlations.phiCorr(A,np.logical_not(A)), -1)

if __name__ == '__main__':
    unittest.main()
