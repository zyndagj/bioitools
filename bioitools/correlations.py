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
	return (n11*n00-n10*n01)/np.sqrt(np.multiply.reduce(nDot, dtype=np.float32), dtype=np.float32)

def phiDist(x,y):
	'''
	Returns 1-(Phi Correlation Coefficient)
	'''
	return 1.0-phiCorr(x,y)

def repToBool(vals):
	out = np.zeros(len(vals),dtype=np.bool)
	out[vals > 0] = 1
	return out

def makeRepliCorr(figType, inFiles, outFile, savePlot, renderPlot, distMethod):
	import bioitools
	labels = np.array(map(lambda x: os.path.splitext(os.path.split(x)[1])[0], inFiles))
	D = []
	if savePlot or renderPlot:
		if not renderPlot:
			import matplotlib
			matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gs
		# Calculate distance
		from multiprocessing import Pool
		p = Pool(16)
		from functools import partial
		pPB = partial(bioitools.ParseBedgraph, vOnly=True)
		ret = p.map(pPB, inFiles)
		p.close()
		p.join()
		#for i in inFiles:
		for i in range(len(inFiles)):
			v = ret[i]
			#c,s,e,v = bioitools.ParseBedgraph(f)
			if distMethod == 'pearson':
				D.append(v)
			elif distMethod == 'phi':
				D.append(repToBool(v))
			else:
				sys.exit("Not a valid distance method")
		nD = np.array(D)
		if distMethod == 'pearson':
			d = distance.pdist(nD, 'correlation')
		elif distMethod == 'phi':
			d = distance.pdist(nD, phiDist)
		else:
			sys.exit("Not a valid distance method")
		pd = distance.squareform(d)
		clusters = h.linkage(pd, method='complete')
		if figType == 'genome':
			# Build figure
			fig = plt.figure(figsize=(15,3))
			hmGS = gs.GridSpec(1,2,wspace=0, hspace=0, left=0.01, right=0.85, width_ratios=[1,15])
			denAX = fig.add_subplot(hmGS[0,0])
			den = h.dendrogram(clusters, orientation="right")
			plt.axis('off')
			hmAX = fig.add_subplot(hmGS[0,1])
			axi = plt.imshow(nD[den['leaves']], aspect='auto', origin='lower', interpolation='nearest', cmap='YlGnBu')
			hmAX.set_yticks(range(nD.shape[0]))
			hmAX.yaxis.set_ticks_position('right')
			hmAX.set_xticklabels("")
			hmAX.set_yticklabels(labels[den['leaves']])
			for l in hmAX.get_xticklines()+hmAX.get_yticklines():
				l.set_markersize(0)
			cbGS = gs.GridSpec(1,1,left=0.97, right=0.98)
			cbAX = fig.add_subplot(cbGS[0,0])
			plt.colorbar(axi, cax=cbAX, use_gridspec=True)
		elif figType == 'matrix':
			# Build figure
			fH = max((int(pd.shape[0]/2),3))
			fig = plt.figure(figsize=(2+fH, fH))
			hmGS = gs.GridSpec(1,2,wspace=0, hspace=0, left=0.01, right=0.75, bottom=0.12, top=(1.0-0.3/fH), width_ratios=[1,fH+1])
			denAX = fig.add_subplot(hmGS[0,0])
			den = h.dendrogram(clusters, orientation="left")
			plt.axis('off')
			hmAX = fig.add_subplot(hmGS[0,1])
			M = pd[den['leaves'], :]
			M = M[:,den['leaves']]
			axi = plt.imshow(1-M, aspect='auto', origin='lower', interpolation='nearest', cmap='YlGnBu', vmin=0, vmax=1)
			#axi = plt.pcolor(1-pd[den['leaves']], cmap='YlGnBu', linewidth=2, edgecolors='w')
			#hmAX.set_yticks(np.arange(nD.shape[0])+0.5)
			hmAX.set_yticks(np.arange(nD.shape[0]))
			hmAX.yaxis.set_ticks_position('right')
			hmAX.set_xticks(np.arange(nD.shape[0]))
			hmAX.set_xticklabels(labels[den['leaves']], rotation=45)
			hmAX.set_yticklabels(labels[den['leaves']])
			plt.title("Sample Correlation")
			for l in hmAX.get_xticklines()+hmAX.get_yticklines():
				l.set_markersize(0)
			cbGS = gs.GridSpec(1,1,left=0.89, right=0.92)
			cbAX = fig.add_subplot(cbGS[0,0])
			plt.colorbar(axi, cax=cbAX, use_gridspec=True)
			# Print correlation matrix
			sLabels=list(labels[den['leaves']])
			print '\t'.join([' ']+sLabels)
			for i in range(len(M)):
				print '\t'.join([sLabels[i]]+map(str, 1-M[i].round(2)))
		else:
			raise ValueError("Incorrect figure type: %s"%(figType))
		if savePlot:
			plt.savefig(outFile)
		if renderPlot:
			plt.show()
