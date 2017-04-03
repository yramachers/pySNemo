# configuration - class to hold configuration parameters for the
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
configuration
==================

This module implements a configuration class which holds all of the
user-configurable parameters available to the cellular automaton and
allows for access to them from a single point.

Note that all angles are stored and returned as the cosine of the angle,
and all distances are in voxel units (i.e. 1 mm)

In addition to the class 'Configuration', there exists a module instance
of this class, called 'config' whose job it is to store the global
configuration information for the cellular automaton. Modules may access
this with the following import statement
"""

import math

__all__ = ['Configuration', 'config']

class Configuration(object):
	"""Class to hold configuration information.
	"""
	def __init__(self):
		"""Initialise a Configuration object with default values
		"""
		self.theta = math.radians(25.0)  # float in range [0,1]
		self.scale_size = 3.0            # float in range [1, inf)
		self.cell_radius_3d = 2.0        # float in range [0, inf)
		self.cell_radius_3d_max = 7.0    # float in range [0, inf)
		self.cell_radius_3d_increment = 0.05 # float in range [0, cell_radius_3d_max)
		self.min_leftward_neighbours = 1 # int - min number of leftward neighbours to terminate adaptive search
		self.comparison_tolerance = 0.05 # float in range [0,inf)
		self.slope_tolerance = 0.10       # float in range [0, inf)
		self.mask_features = False       # Whether to mask out features such as vertices
		self.mask_radius = 10.0          # float in range [0,inf)
		self.mask_origin = False         # Whether to mask out hits around the origin
		self.min_track_length = 3        # int - min. number of hits for cluster to be output

	def getSlopeTolerance(self):
		return self.slope_tolerance

	def setSlopeTolerance(self, tol):
		self.slope_tolerance = float(tol)

	def getRadiusIncrement(self):
		return self.cell_radius_3d_increment

	def setRadiusIncrement(self, inc):
		self.cell_radius_3d_increment = float(inc)

	def getMinLeftwardNeighbours(self):
		"""Return the minimum number of leftward neighbours required
		to terminate an adaptive neighbour search before the maximum
		radius is hit
		"""
		return self.min_leftward_neighbours

	def setMinLeftwardNeighbours(self, nmin):
		self.min_leftward_neighbours = int(nmin)

	def getMinTrackLength(self):
		"""Return the minimum number of hits which must be present
		for a cluster to appear in the final output
		
		Returns:
		--------
		min_track_length : int
			Minimum number of hits for a cluster to be shown in the output
		"""
		return self.min_track_length

	def setMinTracklength(self, min_length):
		"""Backward compatibility function;
		original version had this function only,
		but it is incompatible with the naming scheme

		Just call the new one! It's possible that nothing
		uses this, but leave this in for now.
		"""
		self.setMinTrackLength(min_length)
	
	def setMinTrackLength(self, min_length):
		"""Set the minimum number of hits which must be present
		for a cluster to appear in the final output

		Parameters:
		-----------
		min_length : int
			Minimum number of hits present in a cluster for it to be
			included in the output
		"""
		self.min_track_length = int(min_length)
	
	def getCellRadiusIncrement(self):
		"""Return the incremental step size to be used for 3D cell
		generation
		
		Returns:
		--------
		step : float
			Increment used during adaptive search for cell generation
		"""
		return self.cell_radius_3d_increment

	def setCellRadiusIncrement(self, increment):
		"""Set the incremental step size to be used for 3D cell
		generation

		Parameters:
		-----------
		step : float
			Increment used during adaptive search for cell generation
		"""
		self.cell_radius_3d_increment = float(increment)

	def getCellRadius3DMax(self):
		"""Return the max radius to be used for 3D cell generation

		Returns:
		--------
		radius : float
			Max radius within which points must be located to make cells
		"""
		return self.cell_radius_3d_max

	def setCellRadius3DMax(self, radius):
		"""Set the max radius within which points must be located
		to make 3D cells

		Parameters:
		-----------
		radius : float
			Max radius used for 3D cell generation
		"""
		self.cell_radius_3d_max = float(radius)

	def getCellRadius3D(self):
		"""Return the radius to be used for 3D cell generation

		Returns:
		--------
		radius : float
			Radius within which points must be located to make cells
		"""
		return self.cell_radius_3d

	def setCellRadius3D(self, radius):
		"""Set the radius within which points must be located
		to make 3D cells

		Parameters:
		-----------
		radius : float
			Radius used for 3D cell generation
		"""
		self.cell_radius_3d = float(radius)

	def getMaskRadius(self):
		"""Return the radius to be used for feature masking

		Returns:
		--------
		radius : float
			Radius around a feature point which is to be masked out
		"""
		return self.mask_radius

	def setMaskRadius(self, radius):
		"""Set the radius to be used for feature masking

		Parameters:
		-----------
		radius : float
			Radius around a feature point which is to be masked out
		"""
		self.mask_radius = float(radius)

	def getMaskFeatures(self):
		"""Return whether or not to mask out features such as
		vertices.

		Returns:
		--------
		mask : boolean
			True if features such as vertices should be masked
			out, false otherwise
		"""
		return self.mask_features
	
	def setMaskFeatures(self, mask):
		"""Set whether or not to mask out features such as
		vertices.

		Parameters:
		-----------
		mask : boolean
			True if features such as vertices should be masked
			out, false otherwise
		"""
		self.mask_features = mask

	def getMaskOrigin(self):
		"""Return whether or not to mask out origin

		Returns:
		--------
		mask : boolean
			True if origin should be masked
			out, false otherwise
		"""
		return self.mask_origin
	
	def setMaskOrigin(self, mask):
		"""Set whether or not to mask out origin

		Parameters:
		-----------
		mask : boolean
			True if origin should be masked
			out, false otherwise
		"""
		self.mask_origin = mask

	def getMaxAngle(self):
		"""Return maximum angle allowed between two
		track segments

		Returns:
		--------
		theta : float
			Maximum angle allowed between two track
			segments ('breaking angle') IN RADIANS
		"""
		return self.theta
	
	def setMaxAngle(self, theta):
		"""Set the maximum angle allowed between two
		track segments

		Parameters:
		-----------
		theta : float
			Maximum angle allowed between two track
			segments ('breaking angle') IN RADIANS
		"""
		self.theta = float(theta)

	def getMaxAngleDeg(self):
		"""Return maximum angle allowed between two
		track segments

		Returns:
		--------
		theta : float
			Maximum angle allowed between two track
			segments ('breaking angle') IN DEGREES
		"""
		return math.degrees(self.theta)
	
	def setMaxAngleDeg(self, theta):
		"""Set the maximum angle allowed between two
		track segments

		Parameters:
		-----------
		theta : float
			Maximum angle allowed between two track
			segments ('breaking angle') IN DEGREES
		"""
		self.theta = math.radians(float(theta))

	def getScaleSize(self):
		"""Return the cluster size for scaling operations

		Returns:
		--------
		scale : float
			Scale value to be used for clustering by scaling
		"""
		return self.scale_size

	def setScaleSize(self, scale_size):
		"""Set the cluster size for scaling operations

		Parameters:
		-----------
		scale : float
			Scale value to be used for clustering by scaling
		"""
		self.scale_size = float(scale_size)
	
	def getComparisonTolerance(self):
		"""Return the tolerance for comparison of
		floating-point values

		Returns:
		--------
		tolerance : float
			Tolerance to be used for floating-point
			comparisons
		"""
		return self.comparison_tolerance
	
	def setComparisonTolerance(self, tolerance):
		"""Set the tolerance for comparison of 
		floating-point values

		Parameters:
		-----------
		tolerance : float
			Tolerance to be used for floating-point
			comparisons
		"""
		self.comparison_tolerance = float(tolerance)


# Make an object of type Configuration in the module
config = Configuration()
