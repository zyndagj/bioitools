def parseFai(inFile):
	'''
	Parses a fa.fai into a python dictionary

	Paramteters
	================================
	inFile FILE	fai file
	'''
	return dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), open(inFile,'r').readlines())))
