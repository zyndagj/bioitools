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

def reflectArray(vals, windowSize=20):
	# if w is odd
	if windowSize % 2:
		# extend int(w/2) each side
		wB, wE = (int(windowSize/2), int(windowSize/2))
	else:
		# extend int(w/2)-1 and int(w/2)
		wB, wE = (int(windowSize/2)-1, int(windowSize/2))
	return np.r_[vals[wB:0:-1],vals,vals[-2:-wE-2:-1]], wB, wE

def hann(vals, windowSize=20):
	# Smooths using a hanning window
	if len(vals) < windowSize: return vals
	w = np.hanning(windowSize)
	reflected, wB, wE = reflectArray(vals, windowSize)
	smoothed = np.convolve(w/w.sum(), reflected, mode='valid')
	return smoothed

def average(vals, windowSize=20):
	# Smooths using a moving average window
	if len(vals) < windowSize: return vals
	w = np.ones(windowSize)
	reflected, wB, wE = reflectArray(vals, windowSize)
	smoothed = np.convolve(w/w.sum(), reflected, mode='valid')
	return smoothed

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
	D = {'haar':(haar, p),\
		'hann':(hann, b),\
		'average':(average, b)\
	}
	smoothV = np.zeros(len(vals))
	for chrom in sorted(chromDict.keys()):
		cS, cE = chromDict[chrom]
		cVals = vals[cS:cE]
		tmp = D[method][0](cVals, D[method][1])
		assert(len(tmp) == int(cE-cS))
		smoothV[cS:cE] = tmp
	return smoothV
