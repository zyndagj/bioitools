#!/usr/bin/env python
#
###############################################################################
# Author: Greg Zynda
# Last Modified: 04/23/2018
###############################################################################
# BSD 3-Clause License
# 
# Copyright (c) 2018, Texas Advanced Computing Center
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

def rpccNormalize(inFile, level=1):
	from collections import defaultdict
	from operator import itemgetter
	BGF = open(inFile, 'r')
	noComments = filter(lambda x: x[0] != '#', BGF)
	genome = defaultdict(list)
	# Read Data
	for line in noComments:
		tmp = line.rstrip('\n').split('\t')
		s, e, v = int(tmp[1]), int(tmp[2]), float(tmp[3])
		genome[tmp[0]].append((s,e,v))
	# Process each chrom
	sChroms = sorted(genome.keys())
	for chrom in sChroms:
		total = sum(map(itemgetter(2), genome[chrom]))
		nBases = sum(map(lambda x: x[1]-x[0], genome[chrom]))
		scale = nBases/float(total)
		for record in genome[chrom]:
			s, e, v = record
			oV = v*scale
			if int(oV) == oV:
				print '%s\t%i\t%i\t%i'%(chrom, s, e, oV)
			else:
				print '%s\t%i\t%i\t%.3f'%(chrom, s, e, oV)

def rpgcNormalize(inFile, level=1):
	from collections import defaultdict
	from operator import itemgetter
	BGF = open(inFile, 'r')
	noComments = filter(lambda x: x[0] != '#', BGF)
	genomeList = []
	# Read Data
	for line in noComments:
		tmp = line.rstrip('\n').split('\t')
		c, s, e, v = tmp[0], int(tmp[1]), int(tmp[2]), float(tmp[3])
		genome[tmp[0]].append((c,s,e,v))
	# Process each chrom
	total = sum(map(itemgetter(3), genome[chrom]))
	nBases = sum(map(lambda x: x[2]-x[1], genome[chrom]))
	scale = nBases/float(total)
	for record in genomeList:
			chrom, s, e, v = record
			oV = v*scale
			if int(oV) == oV:
				print '%s\t%i\t%i\t%i'%(chrom, s, e, oV)
			else:
				print '%s\t%i\t%i\t%.3f'%(chrom, s, e, oV)

if __name__ == "__main__":
	main()

