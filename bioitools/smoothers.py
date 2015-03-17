#!/usr/bin/env python

import pywt
import numpy as np

def haar(vals, p=80):
	hC = pywt.wavedec(vals,'haar')
	cutVal = np.percentile(np.abs(np.concatenate(hC)), p)
	for A in hC:
		A[np.abs(A) < cutVal] = 0
	tVals = pywt.waverec(hC,'haar')
	return tVals[:len(vals)]

def hann(vals, windowSize=20):
	# Smooths using a hanning window
	w = np.hanning(windowSize)
	return np.convolve(w/w.sum(), vals, mode='same')

def isContiguous(starts, ends, chromDict):
	'''
	Tests to make sure that the locations in the bedgraph are contiguous and non-overlapping.
	
	Parameters
	==============================
	starts		array of start positions
	ends		array of end positions
	chromDict	dictionary of chromosome bounds
	'''
	for chrom in chromDict:
		cS, cE = chromDict[chrom]
		if not (starts[cS+1:cE]==ends[cS:cE-1]).all():
			return False
	return True

def runSmoother(method, p, b, vals, chromDict):
	if method == 'haar':
		sMethod = haar
		arg = p
	elif method == 'hann':
		sMethod = hann
		arg = b
	smoothV = np.zeros(len(vals))
	for chrom in chromDict:
		cS, cE = chromDict[chrom]
		smoothV[cS:cE] = sMethod(vals[cS:cE], arg)
	return smoothV
