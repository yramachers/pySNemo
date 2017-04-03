# control - Pipeline services for Glue_Cluster
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
ca_glue.control - Control for GlueClusters
===========================================================

This module provides pipeline services for performing
merging of CA clusters if hits in cluster members are close enough.

"""

__all__ = ['GlueService','RoadService']

from pysnemo.ca_glue.glue import MergeXY, RoadSweeper
import logging

class GlueService(object):
	"""
	Pipeline entry for merging of clusters
	
	Uses the clusters from the CA clustering algorithm.
	"""
	def __init__(self, inbox1key, inbox2key, outboxkey):
		"""
		Initialise a GlueService with the parameters given
		inboxkey1 : the key to a list of tuples as input hits
		            for processing CA clusters with 
		            propagation in x
		inboxkey2 : the key to a list of tuples as input hits
		            for processing CA clusters with 
		            propagation in y
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.GlueService')
		self.inkey1 = inbox1key
		self.inkey2 = inbox2key
		self.outkey = outboxkey
		# write logging info
		self.logger.info('Input key 1: %s',self.inkey1)
		self.logger.info('Input key 2: %s',self.inkey2)
		self.logger.info('Output key: %s',self.outkey)
		
	def __repr__(self):
		s = 'CA Glue Service'
		return s


	def __call__(self, event):
		"""
		Process an event with GlueClusters and output the new 
		clusters to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey1 is not None and self.inkey2 is not None):
			if (event.hasKey(self.inkey1) and event.hasKey(self.inkey2)):
				data1 = event.getKeyValue(self.inkey1)
				data2 = event.getKeyValue(self.inkey2)
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			print 'Error in pipeline: Need clustered data as input, not None!'
			return None
		
		# Run the glueing
		if len(data1)>0 and len(data2)>0:
			output = MergeXY(data1, data2).getClusters()
		else:
			if len(data1)>0:
				output = data1
			elif len(data2)>0:
				output = data2
			else:
				output = {}
		
		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)
		
		# Finally, return the event object
		return event



class RoadService(object):
	"""
	Pipeline entry for road sweeping of clusters
	
	Uses the clusters from the CA clustering algorithm.
	"""
	def __init__(self, inboxkey, outboxkey):
		"""
		Initialise a RoadService with the parameters given
		inboxkey : the key to the set of CA clusters.
		outboxkey : output and noise key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.RoadService')
		self.inkey = inboxkey
		self.outkey = outboxkey
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
		
	def __repr__(self):
		s = 'CA Road Sweeping Service'
		return s


	def __call__(self, event):
		"""
		Process an event with RoadSweeper and output the new 
		clusters to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			print 'Error in pipeline: Need clustered data as input, not None!'
			return None
		
		# Run the sweeper
		raw_data = event.getTrackerData()
		s = RoadSweeper(data, raw_data)
		output, noise = s.getClusters()

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)
		
		# noise store in event
		nkey  = self.outkey+'_noise'
		event.setKeyValue(nkey, noise)

		# Finally, return the event object
		return event
