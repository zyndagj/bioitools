#!/usr/bin/env python

import numpy as np
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as h

def phiCorr(x,y):
	x1 = x == 1
	y1 = y == 1
	n11 = np.sum(np.logical_and(x1,y1))
	n10 = np.sum(np.logical_and(x1,np.logical_not(y1)))
	n01 = np.sum(np.logical_and(np.logical_not(x1),y1))
	n00 = np.sum(np.logical_not(np.logical_or(x1,y1)))
	N = np.array([[n11,n10],[n01,n00]])
	nDot = np.hstack((np.sum(N,axis=0),np.sum(N,axis=1)))
	return (n11*n00-n10*n01)/np.sqrt(np.multiply.reduce(nDot))

def repToBool(vals):
	out = np.zeros(len(vals),dtype=np.bool)
	out[vals > 0] = 1
	return out
