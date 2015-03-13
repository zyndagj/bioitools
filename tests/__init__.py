import unittest
import bioitools
from os import path

class TestBioitools(unittest.TestCase):
	def setUp(self):
		self.faiFile = path.join(path.dirname(__file__),'test.fa.fai')
	def test_parseFai(self):
		self.assertEqual(bioitools.ParseFai(self.faiFile), {'Chr1':30427671, 'Chr2':19698289, 'Chr3':23459830})

if __name__ == '__main__':
    unittest.main()
