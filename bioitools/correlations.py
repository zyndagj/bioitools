#!/usr/bin/env python

import numpy as np
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as h
import os

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

def makeRepliCorr(inFiles, outFile, makePlot):
	labels = map(lambda x: os.path.split(x)[1], inFiles)
	D = []
	for f in inFiles:
		c,s,e,v = bioitools.ParseBedgraph(f)
		D.append(repToBool(v))
	nD = np.array(D)
	pd = distance.squareform(distance.pdist(nD, phiCorr))
	clusters = h.linkage(pd, method='complete')
	if args.p:
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gs
		fig = plt.figure(figsize=(10,3))
		hmGS = gs.GridSpec(1,2,wspace=0, hspace=0, width_ratios=[0.2,1])
		denAX = fig.add_subplot(hmGS[0,0])
		den = h.dendrogram(clusters, orientation="right")
		plt.axis('off')
		hmAX = fig.add_subplot(hmGS[0,1])
		axi = plt.imshow(nD[den['leaves']], aspect='auto', origin='lower', interpolation='nearest', cmap='RdBu')
		hmAX.set_yticks(range(nD.shape[0]))
		hmAX.yaxis.set_ticks_position('right')
		hmAX.set_xticklabels("")
		hmAX.set_yticklabels(labels[den['leaves']])
		for l in hmAX.get_xticklines()+hmAX.get_yticklines():
			l.set_markersize(0)
		plt.savefig(args.o)
	print clusters
