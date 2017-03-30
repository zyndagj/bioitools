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
	return np.r_[vals[windowSize-1:0:-1],vals,vals[-2:-windowSize-1:-1]]

def hann(vals, windowSize=20):
	# Smooths using a hanning window
	w = np.hanning(windowSize)
	reflected = reflectArray(vals, windowSize)
	return np.convolve(w/w.sum(), reflected, mode='valid')

def average(vals, windowSize=20):
	# Smooths using a hanning window
	w = np.ones(windowSize)
	reflected = reflectArray(vals, windowSize)
	return np.convolve(w/w.sum(), reflected, mode='valid')

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
	for chrom in chromDict:
		cS, cE = chromDict[chrom]
		cVals = vals[cS:cE]
		tmp = D[method][0](cVals, D[method][1])
		assert(len(tmp) == cE-cS)
		smoothV[cS:cE] = tmp
	return smoothV
