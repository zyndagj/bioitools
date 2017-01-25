#!/usr/bin/env python

import os
import numpy as np
import fileParsers

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

def TranslateSamFlag(flag):
	Message = ("paired in sequencing","a proper pair","unmapped","mate is mapped","forward","mate is reversed","the first read","the second read","not primary","fails qc","a duplicate","a supplrementary alignment")
	outMessage = "Got flag %i\nThe read is "%(flag)
	tArray = map(int,bin(flag)[::-1].split('b')[0])
	blocks = []
	for i in xrange(len(tArray)):
		if tArray[i]:
			blocks.append(Message[i])
	print outMessage+', '.join(blocks)+'.'

def ParseBedgraph(inFile, vOnly=False):
	'''
	Parses a bedgraph (.bedgraph, .bg) into a tuple, where each record is a tuple corresponding to (chr, start, end, value). This will only parse the first four columns.

	Parameters
	================================
	inFile	FILE	bedgraph
	'''
	chrList = []
	startList = []
	endList = []
	valList = []
	if vOnly:
		for chrom, start, end, rest in fileParsers.bedgraph(inFile):
			valList.append(float(rest[0]))
		return np.array(valList, dtype=np.float32)
	else:
		for chrom, start, end, rest in fileParsers.bedgraph(inFile):
			chrList.append(chrom)
			startList.append(start)
			endList.append(end)
			valList.append(float(rest[0]))
		return chrList, np.array(startList,dtype=np.uint32), np.array(endList,dtype=np.uint32), np.array(valList,dtype=np.float32)

def PrintBedgraph(c, s, e, v):
	lC = len(c)
	if len(s) != lC or len(e) != lC or len(v) != lC:
		raise("Something wrong with input to print")
	for i in xrange(lC):
		print "%s\t%u\t%u\t%.4f" %(c[i],s[i],e[i],v[i])
