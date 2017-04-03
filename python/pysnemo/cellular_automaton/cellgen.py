# cellgen - A cell generator for the 3D cellular automaton
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
cellgen : Cell Generator
========================

This module implements a cell generator for the 3D cellular
automaton.
"""

__all__ = ['CellGenerator3D']

import math
import operator

from scipy.spatial import KDTree
from pysnemo.cellular_automaton.configuration import config
from pysnemo.cellular_automaton import utility
from pysnemo.cellular_automaton.cell import Cell
from pysnemo.cellular_automaton.point import Hit

def kahn_slope(points):
	"""Calculate the least square solution to a 3D line fit
	on the direction cosine, from P.C. Kahn, Comp. Chem. 13 (1989) 191
	Input: List of two triplets for 2 points on the line
	"""
	(x,y,z) = [points[0][i] - points[1][i] for i in range(3)]
	distance = math.sqrt(x**2 + y**2 + z**2)
	denom = math.sqrt((x*distance)**2 + (y*distance)**2 + (z*distance)**2)
	alphax = x * distance / denom
	alphay = y * distance / denom
	alphaz = z * distance / denom
	slopexy = alphax * alphay
	slopexz = alphax * alphaz
	slopeyz = alphay * alphaz
	return (slopexy, slopexz, slopeyz)

def float_compare_equal(f1, f2):
	# First check the two numbers
	# have the same sign
	s1 = math.copysign(1.0, f1)
	s2 = math.copysign(1.0, f2)
	if s1 == s2:
		# Now check that |f1| - |f2| < tolerance
		m1 = math.fabs(f1)
		m2 = math.fabs(f2)
		if math.fabs(m1 - m2) < config.getSlopeTolerance():
			return True
	else:
		return False

def slopes_compare_equal(s1, s2):
	for idx in range(3):
		if not float_compare_equal(s1[idx], s2[idx]):
			return False
	return True

def filter_cells_by_slope(cells):
	#print 'Before filtering:', len(cells), 'cells'
	data = [ ]
	for cell in cells:
		left = cell.getFirst()
		right = cell.getLast()
		slope = kahn_slope([left.getCoordinates(), right.getCoordinates()])
		length = right.distanceFrom(left)
		data.append((slope, length, cell))
	# sort by length, shortest --> longest
	data = sorted(data, key=operator.itemgetter(1))
	filtered = [ ]
	while len(data):
		d1 = data.pop()
		s1 = d1[0] # slope of d1
		use_d1 = True
		for d2 in data:
			s2 = d2[0] # slope of d2
			if slopes_compare_equal(s1, s2):
				# d2 is a shorter version of the same cell
				use_d1 = False
				break
		if use_d1:
			filtered.append(d1[2])
	#print 'After filtering:', len(filtered), 'cells'
	return filtered

class CellGenerator3D(object):
	"""Class to generate cells in 3D.
	Used by the 3D CA algorithm
	"""
	def __init__(self, hits, propagation_coordinate):
		# Store things passed in
		self.hits = hits
		
		# build coordinate system
		coords = [0, 1, 2]
		coords.remove(propagation_coordinate)
		self.x = propagation_coordinate
		self.y = coords[0]
		self.z = coords[1]

		# Store internal state
		self.radius = config.getCellRadius3D()
		self.radius_max = config.getCellRadius3DMax()
		self.radius_inc = config.getRadiusIncrement()
		self.min_neighbours = config.getMinLeftwardNeighbours()
		self.tree = self._buildTree()
		self.cells = [ ]
		
		# Build cells
		self._makeCells()
	
	def getCells(self):
		return self.cells
	
	def _buildTree(self):
		# Put data into a format suitable for building a kdtree
		points = [h.getCoordinates() for h in self.hits]
		return KDTree(points)

	def _pointFromIndex(self, index):
		return self.tree.data[index]

	def _pointsFromIndices(self, indices):
		points = [ ]
		for index in indices:
			points.append(self.tree.data[index])
		return points

	def _simpleLR(self, pt, point):
		return pt[self.x] < point[self.x]

	def _fullLR(self, pt, point):
		return ((pt[self.x] < point[self.x])
				or (pt[self.x] == point[self.x] and pt[self.y] < point[self.y])
				or (pt[self.x] == point[self.x] and pt[self.y] == point[self.y] and pt[self.z] < point[self.z]))

	def _pt_leftof(self, pt, point):
		# simple L-R relationship
		return self._simpleLR(pt, point)
		# OR: more complicated L-R relationship
		# return self._fullLR(pt, point)

	def _leftwardNeighbours(self, point, neighbours):
		"""Filter a list of neighbours to leave only
		those left of the given point"""
		left_neighbours = [ ]
		for pt in neighbours:
			if self._pt_leftof(pt, point):
				left_neighbours.append(pt)
		return left_neighbours

	def _expandSearch(self, point):
		left_neighbours = [ ]
		radius = self.radius
		while len(left_neighbours) < self.min_neighbours:
			radius += self.radius_inc
			if radius > self.radius_max:
				break
			indices = self.tree.query_ball_point(point, radius)
			neighbours = self._pointsFromIndices(indices)
			left_neighbours = self._leftwardNeighbours(point, neighbours)
		return left_neighbours

	def _hitFromPoint(self, pt):
		tmp_hit = Hit(pt)
		for h in self.hits:
			if tmp_hit.spatially_equals(h):
				return h
		return None

	def _buildCellsFromNeighbours(self, point, neighbours):
		right = self._hitFromPoint(point)
		cells = [ ]
		for neighbour in neighbours:
			left = self._hitFromPoint(neighbour)
			if left != right:
				c = Cell(left, right)
				cells.append(c)
		return cells

	def _makeCells(self):
		# 1. Query tree for all hits in some radius of each hit
		#  This is the "fast" way; for the most part this should
		#  be correct, but we'll adjust the search for some special
		#  cases, below.
		results = self.tree.query_ball_tree(self.tree, self.radius)
		
		# 2. Iterate over results
		for index in range(len(results)):
			point = self._pointFromIndex(index)
			# 3. Get leftward neighbours
			neighbours = self._pointsFromIndices(results[index])
			left_neighbours = self._leftwardNeighbours(point, neighbours)
			# 4. If we have too few, expand the search
			if len(left_neighbours) < self.min_neighbours:
				# enter expanding adaptive search
				left_neighbours = self._expandSearch(point)
			# 5. Build cells
			cells = self._buildCellsFromNeighbours(point, left_neighbours)
			# 6. If any cells share similar slopes, pick only the smallest
			if len(cells) > 1:
				self.cells.extend(filter_cells_by_slope(cells))
			else:
				self.cells.extend(cells)

