import unittest
import numpy as np
from os import path
import bioitools
from bioitools import smoothers
from bioitools import correlations
from bioitools import fileParsers
from bioitools import normalize
from bioitools import mergers
from bioitools import operations
from bioitools import filters
from StringIO import StringIO
import sys

class TestBioitools(unittest.TestCase):
	def setUp(self):
		self.faiFile = path.join(path.dirname(__file__),'test.fa.fai')
		self.goodBed = path.join(path.dirname(__file__),'test.bedgraph')
		self.good2Bed = path.join(path.dirname(__file__),'test2.bedgraph')
		self.badBed = path.join(path.dirname(__file__),'bad.bedgraph')
	def test_parseFai(self):
		self.assertEqual(fileParsers.fai(self.faiFile), {'Chr1':30427671, 'Chr2':19698289, 'Chr3':23459830})
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
	def test_average(self):
		c,s,e,v = bioitools.ParseBedgraph(self.goodBed)
		np.testing.assert_array_equal(v,np.array([1,2,3,1,2,3,4,0]))
		np.testing.assert_array_equal(smoothers.reflectArray(v,3)[0],np.array([2,1,2,3,1,2,3,4,0,4]))
		ref4 = np.array([2,1,2,3,1,2,3,4,0,4,3])
		np.testing.assert_array_equal(smoothers.reflectArray(v,4)[0], ref4)
		avg4 = np.mean([ ref4[i:i+4] for i in range(len(ref4)-3)], 1)
		np.testing.assert_allclose(smoothers.average(v,3),np.array([5/3.0, 2, 2, 2, 2, 3, 7/3.0, 8/3.0]), atol=0.00001)
		np.testing.assert_allclose(smoothers.average(v,4), avg4, atol=0.00001)
	def test_badFile(self):
		c,s,e,v = bioitools.ParseBedgraph(self.badBed)
		chromDict = bioitools._ChromBounds(c)
		self.assertEqual(smoothers.isContiguous(s,e,chromDict),False)
	def test_phiCorr(self):
		A = np.array([0,0,1,1,0,0],dtype=np.bool)
		self.assertEqual(correlations.phiCorr(A,A), 1)
		self.assertEqual(correlations.phiCorr(A,np.logical_not(A)), -1)
	def test_rpgcNormalize(self):
		sys.stdout = StringIO()
		normalize.rpgcNormalize(self.goodBed)
		retVal = sys.stdout.getvalue()
		goodVal = '''1	0	10	5
1	10	20	10
1	20	30	15
1	30	40	5
1	40	50	10
1	50	60	15
1	60	70	20
1	70	80	0
'''
		self.assertMultiLineEqual(retVal, goodVal)
	def test_rpccNormalize(self):
		# Need to improve this test with another chromosome
		sys.stdout = StringIO()
		normalize.rpccNormalize(self.goodBed)
		retVal = sys.stdout.getvalue()
		goodVal = '''1	0	10	5
1	10	20	10
1	20	30	15
1	30	40	5
1	40	50	10
1	50	60	15
1	60	70	20
1	70	80	0
'''
		self.assertMultiLineEqual(retVal, goodVal)
	def test_merge_mean(self):
		good2 = [self.goodBed, self.goodBed]
		with open(self.goodBed,'r') as F:
			self.assertEqual(list(mergers.mergeBG(good2, method='mean')), map(lambda x: x.rstrip('\n'), F.readlines()))
		goodVal	=	'''1	0	10	2
1	10	20	2
1	20	30	3
1	30	40	2
1	40	50	2
1	50	60	3
1	60	70	4
1	70	80	0'''
		self.assertMultiLineEqual('\n'.join(list(mergers.mergeBG([self.goodBed, self.good2Bed], method='mean'))), goodVal)
	def test_merge_max(self):
		good2 = [self.goodBed, self.goodBed]
		with open(self.goodBed,'r') as F:
			self.assertEqual(list(mergers.mergeBG(good2, method='max')), map(lambda x: x.rstrip('\n'), F.readlines()))
	def test_operation_divide(self):
		good2 = [self.goodBed, self.goodBed]
		goodVal	=	'''1	0	10	1
1	10	20	1
1	20	30	1
1	30	40	1
1	40	50	1
1	50	60	1
1	60	70	1'''
		self.assertMultiLineEqual('\n'.join(list(operations.divide(*good2))), goodVal)
	def test_operation_add(self):
		good2 = [self.goodBed, self.goodBed]
		goodVal	=	'''1	0	10	2
1	10	20	4
1	20	30	6
1	30	40	2
1	40	50	4
1	50	60	6
1	60	70	8
1	70	80	0'''
		self.assertMultiLineEqual('\n'.join(list(operations.add(*good2))), goodVal)
	def test_operation_subtract(self):
		good2 = [self.goodBed, self.goodBed]
		goodVal	=	'''1	0	10	0
1	10	20	0
1	20	30	0
1	30	40	0
1	40	50	0
1	50	60	0
1	60	70	0
1	70	80	0'''
		self.assertMultiLineEqual('\n'.join(list(operations.subtract(*good2))), goodVal)
	def test_operation_multiply(self):
		good2 = [self.goodBed, self.goodBed]
		goodVal	=	'''1	0	10	1
1	10	20	4
1	20	30	9
1	30	40	1
1	40	50	4
1	50	60	9
1	60	70	16
1	70	80	0'''
		self.assertMultiLineEqual('\n'.join(list(operations.multiply(*good2))), goodVal)
	def test_filters_outlying(self):
		goodVal = '''1	60	70	4
1	70	80	0'''
		gen = filters.outlying(self.goodBed, std_mul=1, min_run=2)
		self.assertMultiLineEqual('\n'.join(list(gen)), goodVal)
		gen = filters.outlying(self.goodBed, std_mul=1, min_run=5)
		self.assertMultiLineEqual('\n'.join(list(gen)), '')
		


if __name__ == '__main__':
    unittest.main()
