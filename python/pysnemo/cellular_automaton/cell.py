# cell - Cell objects for use by the cellular automaton
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
cell
==================

This module implements a cell object for use in the cellular automaton.
A cell consists of two point-like objects and one weight (the cell value).

A cell provides information such as whether it is a neighbour of another
cell, the angle between cells, the endpoints of the cell, and the weight
of the cell.
"""

__all__ = ['Cell', 'CellSet']

import math

from pysnemo.cellular_automaton.configuration import config
from pysnemo.cellular_automaton.clustering import Unclustering
from pysnemo.cellular_automaton import utility

class Cell(object):
	"""Class representing a cell used in the cellular automaton.
	A cell consists of two points and a weight.
	"""
	def __init__(self, first, last):
		"""Instantiate a Cell with the given endpoints, and set
		its weight to 1, and set it not to update on the next
		step
		"""
		self.first = first
		self.last = last
		self.weight = 1
		self.update_on_next_step = False
		self.unit_vector = [0, 0, 0]
		self._compute_unit_vector()

	def equal_to(self, other):
		return self.first == other.first and self.last == other.last

	def __repr__(self):
		return '[' + repr(self.first) + '-->' + repr(self.last)  + ']'
	
	def _compute_unit_vector(self):
		"""Compute a unit vector along the line
		defined by the cell
		"""
		self._compute_unit_vector_pca()

	def _compute_unit_vector_pca(self):
		"""Internal method to compute the 
		unit vector along this cell by
		principal components analysis
		"""
		# First get all the unclustered hits
		u = Unclustering()
		hits = u.uncluster(self.first)
		hits.extend(u.uncluster(self.last))

		pca = utility.PCA(hits)
		self.unit_vector = pca.getUnitVector()
		# Make unit vector point from first to last
		direction_vector = self._compute_unit_vector_endpoints_worker()
		sign_direction_x = math.copysign(1.0, direction_vector[0])
		sign_vector_x = math.copysign(1.0, self.unit_vector[0])
		if sign_direction_x != sign_vector_x:
			self.unit_vector = [-x for x in self.unit_vector]

	def _compute_unit_vector_endpoints_worker(self):
		"""Internal method to compute the
		unit vector along this cell by
		taking the centroid coordinates
		of the first and last clusters
		and drawing a line between them
		"""
		q = self.last.getCoordinates()
		p = self.first.getCoordinates()
		diffs = [q[i] - p[i] for i in range(3)]
		len = math.sqrt(sum([d**2 for d in diffs]))
		vec = [d / len for d in diffs]
		return vec

	def _compute_unit_vector_endpoints(self):
		
		self.unit_vector = self._compute_unit_vector_endpoints_worker()
	
	def getFirst(self):
		"""Return the point-like object corresponding to the beginning
		of the cell
		"""
		return self.first
	
	def getLast(self):
		"""Return the point-like object corresponding to the end of
		the cell
		"""
		return self.last
	
	def getWeight(self):
		"""Return the weight of the cell
		"""
		return self.weight
	
	def setUpdate(self):
		"""Set the cell to update on the next CA step
		"""
		self.update_on_next_step = True
	
	def step(self):
		"""Perform the next CA step; updating the weight of this
		cell if required
		"""
		if self.update_on_next_step:
			self.weight += 1
			self.update_on_next_step = False

	def vector(self):
		"""Return a unit vector along the line between first and last
		"""
		return self.unit_vector
	
	def angleWith(self, cell, invert=False):
		"""Return the cosine of the angle between this cell and another cell

		Returns:
		--------
		cos_theta : float
			Result of scalar product between unit vectors along each
			cell
		"""
		u = self.vector()
		v = cell.vector()
		#print 'Vectors:', u, v
		cos_angle = utility.compute_cosine_between(u, v)
		if invert:
			return 180.0 - math.acos(cos_angle)
		else:
			return math.acos(cos_angle)

	def isLeftwardNeighbourOf(self, cell, propagation_coordinate):
		mid_left = self.getLast()
		mid_right = cell.getFirst()
		# Test if there is a shared point
		if mid_left == mid_right:
			far_left = self.getFirst()
			far_right = cell.getLast()
			# Test if the other two points differ
			if far_left != far_right:
				theta = math.fabs(self.angleWith(cell))
				theta_max = config.getMaxAngle()
				if theta <= theta_max:
					return True
		return False


class CellSet(object):
	"""Set of Cell objects. Has the usual properties of a 
	(mathematical) set
	"""
	def __init__(self):
		"""Construct an empty CellSet
		"""
		self.cells = [ ]
		self.maxcell = None
		self.maxweight = 0

	def __repr__(self):
		s = ', '.join((repr(cell) for cell in self.cells))
		return '[' + s + ']'

	def __len__(self):
		"""Special method required to implement the len(x) operation
		for a CellSet
		"""
		return len(self.cells)
	
	def __iter__(self):
		"""Special method producing a generator to iterate
		over the elements in the CellSet
		"""
		for cell in self.cells:
			yield cell

	def append(self, cell):
		"""Add a cell to the CellSet
		"""
		if cell not in self.cells:
			weight = cell.getWeight()
			if weight > self.maxweight:
				self.maxcell = cell
				self.maxweight = weight
			self.cells.append(cell)
	
	def extend(self, cells):
		"""Add cells to the CellSet
		"""
		for cell in cells:
			self.append(cell)

	def getCells(self):
		"""Return the number of cells in the set
		"""
		return self.cells
	
	def cellsNotIn(self, other):
		"""Return the set of cells in this set but not in the other set
		A convenience function for improving program readability

		Returns:
		--------
		remaining : CellSet
			New set containing cells from this set which are not in the
			other set
		"""
		new_set = CellSet()
		new_set.extend((cell for cell in self.cells if cell not in other.cells))
		return new_set

	def getMaxCell(self):
		return self.maxcell	

