# point - Various point-like objects
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
Point-like objects
==================

This module implements a number of point-like objects:
	- Hit: A single space-point taken to be centered on a unit voxel
	  in 3D
	- Cluster: A collection of hits or clusters whose position is the
	  centroid position of all the sub-points.
	- Cell: A pair of point-like objects representing a single cell
	  in the cellular automaton. In addition to spatial properties,
	  cells have a 'weight' parameter representing their current state.
"""

__all__ = ['Hit', 'Cluster']

import math

from pysnemo.cellular_automaton.configuration import config

def distance_between(p0, p1):
	"""Returns the Euclidean distance between two 3D point-like objects.

	Parameters:
	-----------
	p0 : point-like object
	p1 : point-like object

	Returns:
	--------
	distance : float
		Euclidean distance between the two points in 3D
	"""
	x0 = p0.getCoordinates()
	x1 = p1.getCoordinates()
	r2 = sum([(x1[i] - x0[i])**2 for i in range(3)])
	return math.sqrt(r2)

def distance_between_2d(p0, p1, x, y):
	"""Compute the Euclidean distance between two point-like objects
	in 2D, using the coordinate indices given by x and y.

	Parameters:
	-----------
	p0 : point-like object
	p1 : point-like object
	x : index in (0,1,2) representing first coordinate
	y : index in (0,1,2) representing second coordinate

	Returns:
	--------
	distance : float
		Euclidean distance between the two points in a 2D projection
	"""
	x0 = p0.getCoordinates()
	x1 = p1.getCoordinates()
	r2 = (x1[x] - x0[x])**2 + (x1[y] - x0[y])**2
	return math.sqrt(r2)

class Hit(object):
	"""Class representing a single voxel as a point-like object
	with coordinates specifying its centre.
	"""
	def __init__(self, coordinates, radius=1.0, errors=[]):
		"""Instantiate a Hit with the supplied coordinates
		"""
		self.coordinates = coordinates
		self.errors = errors
		self.radius = radius
		self.meta_data = None

	def getRadius(self):
		return self.radius

	def getErrors(self):
		return self.errors

	def setMeta_data(self, gginfo):
		self.meta_data = gginfo
	
	def getMeta_data(self):
		return self.meta_data
	
	def hasMeta_data(self):
		if self.meta_data is None:
			return False
		else:
			return True
	
	def __repr__(self):
		parts = ', '.join((str(x) for x in self.coordinates))
		return '('+parts+')'

	def getCoordinates(self):
		"""Return coordinates as 3-element list of float
		"""
		return self.coordinates
	
	def getPoints(self):
		"""Return list of all contained points; since a Hit is not
		clustered, returns a list with just one element
		"""
		return [self]
	
	def distanceFrom(self, point):
		"""Get Euclidean distance in 3D between this Hit and another
		point-like object
		"""
		return distance_between(self, point)
	
	def distanceFrom2D(self, point, x, y):
		"""Get Euclidean distance in 2D between this Hit and another
		point-like object, using only the coordinate indices given by
		x and y
		"""
		return distance_between_2d(self, point, x, y)

		
	def isClustered(self):
		"""Returns False, indicating that a Hit is not a clustered object
		"""
		return False
	
	def isUnclustered(self):
		"""Returns True, indicating that a Hit is not a clustered object
		Opposite of Hit.isClustered()
		"""
		return not self.isClustered()
	
	def spatially_equals(self, other):
		p = self.coordinates
		q = other.getCoordinates()
		if p[0] == q[0] and p[1] == q[1] and p[2] == q[2]:
			return True
		else:
			return False

		
class Cluster(object):
	"""Class representing a collection of point-like objects
	(hits or other clusters)
	
	The coordinates of a cluster are the centroid of the coordinates
	of the point-like objects contained within it.
	"""
	def __init__(self, scaled=False):
		"""Instantiate a Cluster, which initially holds no points,
		and has invalid centroid value (represented by None)
		"""
		self.points = [ ]
		self.accumulated = [0.0, 0.0, 0.0] # Always holds unscaled values
		self.charge = 0.0
		self.coordinates = None # holds charge weighted and possibly scaled values
		self.coords_valid = False
		self.scaled = scaled
	
	def __repr__(self):
		if not self.coords_valid:
			self.getCoordinates() # Ensure coordinates are valid
		parts = ', '.join((str(x) for x in self.coordinates))
		return '('+parts+')'

	def __len__(self):
		return len(self.points)

	def _compute_coords(self):
		if self.charge>0.0:
			self.coordinates = [x / self.charge for x in self.accumulated]
			if self.scaled:
				factor = float(config.getScaleSize())
				self.coordinates = [x / factor for x in self.coordinates]
		self.coords_valid = True

	def setScaled(self, scaled):
		self.scaled = scaled
		self.coords_valid = False
		self._compute_coords()
	
	def isScaled(self):
		return self.scaled

	def getCoordinates(self):
		"""Returns the coordinates of the centroid of the cluster.
		If the cluster is empty, return None
		"""
		if not self.coords_valid:
			self._compute_coords()
		return self.coordinates
	
	def getPoints(self):
		"""Return a list of points held by the cluster
		"""
		return self.points
	
	def getCharge(self):
		"""Return the total charge of this cluster
		"""
		return self.charge

	def distanceFrom(self, point):
		"""Return the Euclidean distance in 3D between this cluster
		and another point-like object
		"""
		return distance_between(self, point)
	
	def distanceFrom2D(self, point, x, y):
		"""Returm the Euclidean distance in 2D between this cluster
		and another point-like object, where x and y give the indices
		of the spatial coordinates for the projection to be used
		"""
		return distance_between_2d(self, point, x, y)
	
	def isClustered(self):
		"""Return True, indicating that this object is a cluster of
		point-like objects
		"""
		return True

	def isUnclustered(self):
		"""Return False, indicating that this object is a cluster of
		point-like objects (opposite of Cluster.isClustered())
		"""
		return not self.isClustered()
	
	def addPoint(self, point):
		"""Add a point to the cluster.
		If the point is already in the cluster, do nothing.
		"""
		for point2 in self.points:
			if point.spatially_equals(point2):
				break
		else: # point wasn't found in the current set of points
			self.coords_valid = True
			self.coordinates = point.getCoordinates()
			self.points.append(point)
			
	def addPoints(self, points):
		"""Add many points to the cluster.
		Use self.addPoint() so as to omit points which are
		already in the cluster.
		"""
		for point in points:
			self.addPoint(point)
				
	def getSize(self):
		"""Return the number of points in the cluster
		"""
		return len(self.points)

	def spatially_equals(self, other):
		p = self.coordinates
		q = other.getCoordinates()
		if p[0] == q[0] and p[1] == q[1] and p[2] == q[2]:
			return True
		else:
			return False
	
