# control - Pipeline services for Clustering
#
# Copyright (c) 2013, YR
# Copyright (c) 2012 The University of Warwick
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
cluster.control - Control for clustering algorithms
===========================================================

This module provides pipeline services for performing
density-based clustering using DBSCAN as well as the set-cluster idea, a 
point cloud service, filament finder and image segmentation.

"""

__all__ = ['FilamentService']

from pysnemo.cluster.filamentsearch import FilamentCluster
import logging


class FilamentService(object):
	"""
	Pipeline entry for filament clustering

	Uses the raw hits only and can serve as preparation clustering for 
	more detailed clusterers subsequently.
	"""
	def __init__(self, inboxkey, outboxkey):
		"""
		Initialise a FilamentService with the parameters given

		inboxkey : the key to a list of tuples as input hits
		            for processing with ImageSegmentation - None
			    means - get the raw hits.
			    Must be list of tracker_hit objects.
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.FilamentgService')
		self.inkey = inboxkey
		self.outkey = outboxkey

		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
	
	def __repr__(self):
		s = 'Filament Cluster Service'
		return s


	def __call__(self, event):
		"""
		Process an event with Filament Search and output the 
		clusters to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
				if not isinstance(data,list):
					print 'Error in pipeline: not a list type input data as expected for raw tracker hits.'
					return None
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			data = event.getTrackerData()

		if len(data)>0:
			cand = self.process(data) # here call the filament search algorithm

		else: # empty event
			event.setKeyValue(self.outkey, {})
			return event

		# Then set the outboxkey
		event.setKeyValue(self.outkey, cand)

		# Finally, return the event object
		return event


	def process(self,data):
		# Now run 
		filament = FilamentCluster(data)
		filament.run()
		return filament.getClusters() # come as dictionary of numbered clusters

