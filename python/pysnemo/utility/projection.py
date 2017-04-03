# projection.py - Utilities for projecting events into 1D histograms
#
# Copyright (c) 2011 Andrew J. Bennieston <A.J.Bennieston@warwick.ac.uk>
# Copyright (c) 2011 The University of Warwick
#
# This file is part of pySNemo.
#
# pySNemo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pySNemo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with pySNemo.  If not, see <http://www.gnu.org/licenses/>.

"""
projection : Utilities for projecting events onto 1D lines
=======================

"""
import math

from ROOT import TH1I

__all__ = ['Projector']

class Projector(object):
	def __init__(self, hits=None):
		self.hits = hits

	def _gen_hist(self, hitlist, coord):
		minval = min(hitlist) - 1
		maxval = max(hitlist) + 1
		nbins = int(math.ceil((maxval - minval) / 3.0))
		name = 'proj_' + coord
		title = 'proj_' + coord + ';' + coord + ';nhits'
		hist = TH1I(name, title, nbins, minval, maxval)
		for hit in hitlist:
			hist.Fill(hit)
		return hist

	def projectX(self):
		if len(self.hits):
			x_hits = [h[0] for h in self.hits]
			return self._gen_hist(x_hits, 'x')

	def projectY(self):
		if len(self.hits):
			y_hits = [h[1] for h in self.hits]
			return self._gen_hist(y_hits, 'y')

	def projectZ(self):
		if len(self.hits):
			z_hits = [h[2] for h in self.hits]
			return self._gen_hist(z_hits, 'z')

	def getProjectedMaxima(self):
		hx = self.projectX()
		hy = self.projectY()
		hz = self.projectZ()
		xmax = hx.GetMaximum()
		ymax = hy.GetMaximum()
		zmax = hz.GetMaximum()
		return (xmax, ymax, zmax)

	def getCoordOfLeastMaximum(self):
		maxima = self.getProjectedMaxima()
		sorted_list = sorted(enumerate(maxima), key=lambda x: x[1])
		min_coord = sorted_list[0][0]
		return min_coord

