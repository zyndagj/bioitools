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

def makeRepliCorr(figType, inFiles, outFile, savePlot, renderPlot, distMethod, titleString):
	import bioitools
	labels = np.array(map(lambda x: os.path.splitext(os.path.split(x)[1])[0], inFiles))
	D = []
	if savePlot or renderPlot:
		import matplotlib
		if not renderPlot:
			matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		import matplotlib.gridspec as gs
		def get_font_width(nChars, font_size, dpi):
			f = plt.figure(1, figsize=(10,10), dpi=dpi)
			r = f.canvas.get_renderer()
			t = plt.text(0.5, 0.5, 'G'*nChars, fontsize=font_size)
			right_width = t.get_window_extent(renderer=r).width
			plt.close(f)
			return right_width
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
			from matplotlib import rcParams
			rcParams['pdf.fonttype'] = 42
			rcParams['svg.fonttype'] = 'none'
			old_settings = np.seterr(all='ignore')
			# Number of inputs
			num_rows = pd.shape[0]
			maxLabelChars = max(map(len, labels))
			dpi = 100 if renderPlot else 200
			padding = dpi/8
			bottom_pad = padding*3
			den_width = dpi
			main_px = dpi*2/3*num_rows if num_rows <= 4 else dpi/2*num_rows
			#font_size = 12 if num_rows <= 4 else 9
			font_size = max(min(int(12 - 0.25 * (num_rows-6)), 13), 8)
			right_width = get_font_width(maxLabelChars, font_size, dpi)
			top_height = np.ceil(float(right_width)/np.sqrt(2))
			title_space = dpi*1/2
			cbar_height = dpi/4
			total_width_px = padding + den_width + padding + main_px + padding + right_width + padding
			tW = float(total_width_px)
			total_height_px = bottom_pad + cbar_height + padding + main_px + padding + top_height + title_space 
			tH = float(total_height_px)
			# Initialize figure
			fig = plt.figure(1,figsize=(total_width_px/dpi, total_height_px/dpi), dpi=dpi)
			#fig = plt.figure(figsize=(11, 9.5)) if num_rows < 20 else plt.figure(figsize=(11.0*1.5, 9.5*1.5))
			# Write title
			if titleString:
				plt.suptitle(titleString)
			else:
				plt.suptitle("Sample Correlation")
			plt.axis('off')
			# Caluclate and set font size
			#font_size = max(min(int(13 - 0.25 * (num_rows-6)), 13), 5)
			rcParams.update({'font.size': font_size})
			# Make dendrogram
			axdendro = fig.add_axes([padding/tW, (bottom_pad+cbar_height+padding)/tH, den_width/tW, main_px/tH])
    			#axdendro.set_axis_off()
			den = h.dendrogram(clusters, orientation="left")
			plt.axis('off')
			#axdendro.set_xticks([])
			#axdendro.set_yticks([])
			# Sort matrix
			index = den['leaves']
			M = 1-pd[index, :]
			M = M[:, index]
			# Calculate vmin
			vmin = 0 if M.min() >= 0 else -1
			# Draw color matrix
			axmatrix = fig.add_axes([(padding*2+den_width)/tW, (bottom_pad+cbar_height+padding)/tH, main_px/tW, main_px/tH])
			#axmatrix = fig.add_axes([0.13, 0.1, 0.6, 0.7])
			cmatrix = axmatrix.pcolormesh(M, edgecolors='white', cmap='PiYG', vmax=1, vmin=vmin)
			plt.box(on=None)
    			axmatrix.set_xlim(0, num_rows)
    			axmatrix.set_ylim(0, num_rows)
			# Write yticks
    			axmatrix.yaxis.tick_right()
    			axmatrix.set_yticks(np.arange(num_rows) + 0.5)
    			axmatrix.set_yticklabels(labels[index])
			# Write xticks
			axmatrix.xaxis.set_tick_params(labeltop='on')
			axmatrix.xaxis.set_tick_params(labelbottom='off')
			axmatrix.set_xticks(np.arange(M.shape[0]) + 0.5)
			axmatrix.set_xticklabels(labels[index],rotation=45,ha='left')
			# Style ticks
			axmatrix.tick_params(axis='x', which='both', bottom='off', top='off')
			axmatrix.tick_params(axis='y', which='both', left='off', right='off')
			# Plot colorbar.
			axcolor = fig.add_axes([(padding*2+den_width)/tW, (bottom_pad)/tH, main_px/tW, cbar_height/tH])
			#axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
			cbar = plt.colorbar(cmatrix, cax=axcolor, orientation='horizontal')
			cbar.solids.set_edgecolor("face")
			for row in range(num_rows):
				for col in range(num_rows):
					textColor = 'white' if abs(M[row,col]) >= 0.5 else 'black'
					axmatrix.text(row + 0.5, col + 0.5, "{:.2f}".format(M[row, col]), ha='center', va='center', color=textColor)
			# Print correlation matrix
			sLabels=list(labels[den['leaves']])
			print '\t'.join([' ']+sLabels)
			for i in range(len(M)):
				print '\t'.join([sLabels[i]]+map(str, M[i].round(2)))
		else:
			raise ValueError("Incorrect figure type: %s"%(figType))
		if savePlot:
			plt.savefig(outFile)
		if renderPlot:
			plt.show()
