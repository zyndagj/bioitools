def fastq(inFile):
	IF = open(inFile,'r')
	name1 = IF.readline()
	while name1:
		seq = IF.readline().rstrip('\n')
		name2 = IF.readline()
		qual = IF.readline().rstrip('\n')
		yield((name1,seq,qual))
		name1 = IF.readline()
	IF.close()

def fasta(infile):
	header = ''
	seqBuild = ''
	for line in open(infile,'r'):
		tmp = line.rstrip('\n')
		if tmp[0] == '>':
			if header and seqBuild:
				yield((header, seqBuild))
			seqBuild = ''
			header = tmp[1:]
		else:
			seqBuild += tmp
	yield((header, seqBuild))

def fai(inFile):
	'''
	Parses a fa.fai into a python dictionary

	Paramteters
	================================
	inFile	FILE	fai file
	
	Output
	================================
	{chr:size, chr:size... }
	'''
	return dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), open(inFile,'r').readlines())))

def bedgraph(inFile):
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		chrom = tmp[0]
		start = int(tmp[1])
		end = int(tmp[2])
		rest = tmp[3:]
		yield((chrom, start, end, rest))
