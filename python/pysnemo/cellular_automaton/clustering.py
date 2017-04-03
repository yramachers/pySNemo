# clustering - Module for clustering hits into larger groups
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
clustering
==========

This module provides services for clustering hits together into
larger groups in one of two ways:

1. Scaling
	Clusters hits together according to some scale size (from the
	global configuration)
2. Orthogonal Clustering
	Clusters hits together orthogonal to the current propagation
	direction, subject to certain criteria
"""

import math
from pysnemo.cellular_automaton.point import Cluster
from pysnemo.cellular_automaton.configuration import config

__all__ = ['Scaling', 'Unclustering']

class Scaling(object):
	"""Class providing access to clustering by scaling up the voxel size
	and re-interpreting coordinates in terms of the scaled coordinate system

	e.g. a point at (4,3,5) becomes a member of a cluster with coordinates
	determined by the charge-weighted centroid of that cluster
	"""
	def __init__(self, points):
		"""Initialise a Scaling object with a list of point-like objects
		and clusters the points
		"""
		self.points = points
		self.cluster_map = { }
		self._do_clustering()
	
	def _do_clustering(self):
		"""Member function to perform the clustering on initialisation
		"""
		for point in self.points:
			coords = point.getCoordinates()
			scale_size = config.getScaleSize()
			scaled = tuple([math.floor(x / scale_size) for x in coords])
			if scaled not in self.cluster_map:
				cluster = Cluster(scaled=True)
				self.cluster_map[scaled] = cluster
			self.cluster_map[scaled].addPoint(point)
	
	def getClusters(self):
		"""Return a list of clusters resulting from scaling

		Returns:
		--------
		clusters : list of Cluster objects
			Clusters resulting from the scaling operation
		"""
		scaled = self.cluster_map.values()
		output = [p for p in scaled if p.getCoordinates() is not None]
		return output

		
class Unclustering(object):
	"""Class providing a service which turns a cluster into a list of top-level
	(primary) hits.
	"""
	def __init__(self):
		"""Create an unclustering service which can be used to uncluster one or
		more clusters of hits
		"""
		pass
	
	def uncluster(self, cluster):
		"""Uncluster hits by operating recursively on clusters and adding their
		hits to a list
		"""
		hits = [ ]
		for item in cluster.getPoints():
			if item.isClustered():
				hits.extend(self.uncluster(item))
			else:
				hits.append(item)
		return hits

