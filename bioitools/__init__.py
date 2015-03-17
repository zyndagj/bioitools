#!/usr/bin/env python

import os
import numpy as np

def _ChromBounds(chromList):
	'''
	Calculates the locations of each chromosome in a bed file.
	'''
	chromDict = {}
	current = chromList[0]
	first = 0
	for i in xrange(len(chromList)):
		if chromList[i] != current:
			chromDict[current] = (first,i)
			current = chromList[i]
			first = i
	chromDict[current] = (first,i+1)
	return chromDict

def ParseBedgraph(inFile):
	'''
	Parses a bedgraph (.bedgraph, .bg) into a tuple, where each record is a tuple corresponding to (chr, start, end, value). This will only parse the first four columns.

	Parameters
	================================
	inFile	FILE	bedgraph
	'''
	_checkFile(inFile,[".bedgraph",".bg"])
	chrList = []
	startList = []
	endList = []
	valList = []
	for line in open(inFile,'r'):
		tmp = line.rstrip('\n').split('\t')
		chrList.append(tmp[0])
		startList.append(int(tmp[1]))
		endList.append(int(tmp[2]))
		valList.append(float(tmp[3]))
	return chrList, np.array(startList,dtype=np.uint32), np.array(endList,dtype=np.uint32), np.array(valList,dtype=np.float32)

def PrintBedgraph(c, s, e, v):
	lC = len(c)
	if len(s) != lC or len(e) != lC or len(v) != lC:
		raise("Something wrong with input to print")
	for i in xrange(lC):
		print "%s\t%u\t%u\t%.4f" %(c[i],s[i],e[i],v[i])

def ParseFai(inFile):
	'''
	Parses a fa.fai into a python dictionary

	Paramteters
	================================
	inFile	FILE	fai file
	'''
	_checkFile(inFile,[".fai"])
	return dict(map(lambda y: (y[0], int(y[1])), map(lambda y: y.split('\t'), open(inFile,'r').readlines())))

def _checkFile(inFile, extList):
	if not os.path.exists(inFile):
		raise ValueError("File doesn't exist")
	if os.path.splitext(inFile)[1] not in extList:
		raise ValueError("Expects a %s file."%(extension))
