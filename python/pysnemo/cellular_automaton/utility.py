# utility - module to hold utility functions for the
# cellular automaton and provide a single access point
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
utility
==================

This module implements utility functions which may be used anywhere
in the cellular automaton package.

These utility functions include:

	float_compare_equal(f1, f2) for comparing two floats
	coordinates_compare_equal(p1, p2) for comparing two space points
	compute_cosine_between(u, v) for computing cos theta between two unit vectors
	class PCA for principal component analysis
"""

import math

import numpy
import scipy
from heapq import nlargest

from pysnemo.cellular_automaton.configuration import config

def float_compare_equal(f1, f2):
	# First check the two numbers
	# have the same sign
	copysign = math.copysign
	fabs = math.fabs
	s1 = copysign(1.0, f1)
	s2 = copysign(1.0, f2)
	if s1 == s2 and math.fabs(f1 - f2) < config.getComparisonTolerance():
		return True
	else:
		return False

def coordinates_compare_equal(p1, p2):
	for idx in xrange(len(p1)):
		if not float_compare_equal(p1[idx], p2[idx]):
			return False
	return True

def compute_cosine_between(u, v):
	# Compute angle between two unit vectors
	cos_theta = sum([u[i]*v[i] for i in range(3)]) 
	if cos_theta > 1.0:
		cos_theta = 1.0
	elif cos_theta < -1.0:
		cos_theta = -1.0
	return cos_theta

class PCA(object):
	def __init__(self, hits):
		transposed = zip( *map(lambda hit: hit.getCoordinates(), hits) )
		data_array = numpy.array(transposed)
		# Get eigenvalues and eigenvectors
		eigenval, eigenvec = numpy.linalg.eig(scipy.cov(data_array))
		# Transpose eigenvec to return to dataset
		etranspose = numpy.transpose(eigenvec)
		# find index of largest element
		indexes = [0, 1, 2]
		largest_index = nlargest(1, indexes, key=lambda i: eigenval[i])
		# Now the principal component of the data set lies along the vector
		# given by eigenvec[largest_index]
		vec = etranspose[largest_index][0]
		self.unit_vector = [numpy.real(vec[0]), numpy.real(vec[1]), numpy.real(vec[2])]
	
	def getUnitVector(self):
		return self.unit_vector

class PCATuple(object):
	def __init__(self, hits):
		data_array = numpy.transpose(numpy.array(hits))
		# Get eigenvalues and eigenvectors
		eigenval, eigenvec = numpy.linalg.eig(scipy.cov(data_array))
		# Transpose eigenvec to return to dataset
		etranspose = numpy.transpose(eigenvec)
		# find index of largest element
		indexes = [0, 1, 2]
		largest_index = nlargest(1, indexes, key=lambda i: eigenval[i])
		# Now the principal component of the data set lies along the vector
		# given by eigenvec[largest_index]
		vec = etranspose[largest_index][0]
		self.unit_vector = [numpy.real(vec[0]), numpy.real(vec[1]), numpy.real(vec[2])]
	
	def getUnitVector(self):
		return self.unit_vector
