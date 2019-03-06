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

from bioitools.fileParsers import bedgraph
from bioitools import ParseBedgraph
import numpy as np

def outlying(inFile, method='mean', std_mul=2, min_run=5, ddof=1, sqrt=True):
	agg_dict = {'mean':np.mean, 'median':np.median}
	agg = agg_dict[method]
	vA = ParseBedgraph(inFile, vOnly=True)
	if sqrt: vA = np.sqrt(vA)
	midVal = agg(vA)
	std = np.std(vA, ddof=ddof)
	minThresh = midVal-std_mul*float(std)
	maxThresh = midVal+std_mul*float(std)
	tmpStack = []
	for c,s,e,v in bedgraph(inFile):
		res = float(v[0])
		tres = res
		if sqrt:
			tres = np.sqrt(res)
		if tres > maxThresh or (tres < minThresh and tres > 0):
			if int(res) == res:
				outStr = '%s\t%i\t%i\t%i'%(c,s,e,res)
			else:
				outStr = '%s\t%i\t%i\t%f'%(c,s,e,res)
			tmpStack.append(outStr)
		else:
			if len(tmpStack) >= min_run:
				for i in tmpStack: yield i
			tmpStack = []
	if len(tmpStack) >= min_run:
		for i in tmpStack: yield i
