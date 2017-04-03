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
point cloud service and image segmentation.

"""

__all__ = ['DBSCANService','NNSetClusterService','FlattenNN','PointCloudService','ImageSegService']

from pysnemo.cluster.dbscan import dbscan, dbprocess
from pysnemo.cluster.nnsetcluster import NNSetCluster
from pysnemo.cluster.pointcloud import PointCloud
from pysnemo.cluster.imagesegmentation import ImageSegmentation
import logging

class DBSCANService(object):
	"""
	Pipeline entry for density-based clustering using DBSCAN

	Uses the raw hits 
	which is to be in the form list of hits, i.e. [(x,y,z,Q), ...]

	See the dbscan.py file for details of output format as a 
	dictionary.
	"""
	def __init__(self, inboxkey, outboxkey, epsilon, min_points):
		"""
		Initialise a DBSCANService with the parameters given
		epsilon : neighbourhood radius, max is set to 12.6 in 2D
		min_points : int min number of points within epsilon radius
		inboxkey : the key to a list of tuples as input hits
		            for processing with dbscan - None
			    means - get the raw hits.
			    Must be list of geigercylinder objects.
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.DBSCANService')
		self.epsilon = epsilon
		self.min_points = min_points
		self.inkey = inboxkey
		self.outkey = outboxkey
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
		self.logger.info('Epsilon: %f',self.epsilon)
		self.logger.info('Minimum points: %d',self.min_points)
	
	def __repr__(self):
		s = 'DBScan Service'
		return s


	def __call__(self, event):
		"""
		Process an event with DBSCAN and output the clustered
		hits to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			data = event.getTrackerData()

		# Now run DBSCAN, first 2D only
		output2D, noise2D = self.adaptive2D(data)
		print 'Found %d clusters in adaptive2D.'%len(output2D)
		if (len(output2D)<1): # only if null cluster, try again
			output,noise = self.adaptive3D(data)
		else:
			output = output2D
			noise = noise2D

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)

		# store in event
		nkey  = self.outkey+'_noise'
		event.setKeyValue(nkey, noise)
	
		# Finally, return the event object
		if (len(output)>0):
			return event
		else:
			return None # stops the pipeline


	def _complement(self,l,length):
		A = set(l)
		Bl = [i for i in range(length)]
		B = set(Bl)
		return list(B.difference(A))


	def adaptive2D(self,data):
		# First, get wire coordinates for each gcylinder, no z value
		eps = self.epsilon
		spatial_only = [(gg.xwire,gg.ywire,0.0) for gg in data] 
		clusters = {}
		while (len(clusters)<1 and eps < 12.6):
			clusters = dbscan(spatial_only, eps, self.min_points)
			output = { }
			noise = []
			print 'In adaptive2D, eps=',eps
			if (len(clusters)):
				# build the dictionary of clusters for the event storage
				for cluster_id, indices in clusters.iteritems():
					out = []
					for ind in list(indices[0]):
						out.append(data[ind])
					output[cluster_id] = out

				# Noise - not loosing any wires
				outind = []
				for cluster_id, indices in clusters.iteritems():
					outind.extend(list(indices[0])) # collect all
				compl = self._complement(outind,len(data)) 
				for ind in compl:   # not in outind 
					noise.append(data[ind])
			else:
				eps += 0.5 # next iteration
		return output, noise


	def adaptive3D(self, data):
		# First, get wire coordinates for each gcylinder
		eps = self.epsilon
		spatial_only = [(gg.xwire,gg.ywire,gg.zwire) for gg in data] 
		clusters = {}
		while (len(clusters)<1 and eps < 100.0): # relax z cond. 3D
			clusters = dbscan(spatial_only, eps, self.min_points)
			output = { }
			noise = []
			print 'In adaptive3D, eps=',eps
			if (len(clusters)):
				# build the dictionary of clusters for the event storage
				for cluster_id, indices in clusters.iteritems():
					out = []
					for ind in list(indices[0]):
						out.append(data[ind])
					output[cluster_id] = out

				# Noise - not loosing any wires
				outind = []
				for cluster_id, indices in clusters.iteritems():
					outind.extend(list(indices[0])) # collect all
				compl = self._complement(outind,len(data)) 
				for ind in compl:   # not in outind 
					noise.append(data[ind])
			else:
				eps += 1.0 # next iteration
		return output, noise


		



class LegendreService(object):
	"""
	Pipeline entry for Legendre clustering

	Uses the raw hits 
	which is to be in the form list of hits, i.e. [(x,y,z,Q), ...]

	See the leg_cluster.py file for details of output format as a 
	dictionary. Should contain line_candidate objects, defined in 
	leg_cluster.py.
	"""
	def __init__(self, inboxkey, outboxkey, magnet=False):
		"""
		Initialise a LegendreService with the parameters given
		inboxkey : the key to a list of tuples as input hits
		            for processing - None
			    means - get the raw hits.
			    Must be list og geigercylinder objects
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.LegendreService')
		self.clusterflag = False # flag cluster data
		self.inkey = inboxkey
		self.outkey = outboxkey
		self.flag = magnet # default to lines
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
		self.logger.info('Magnet flag: %s',self.flag)
	
	def __repr__(self):
		s = 'Legendre Cluster Service'
		return s


	def __call__(self, event):
		"""
		Process an event with leg_cluster and output the clustered
		hits to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
				if isinstance(data,dict):
					self.clusterflag = True
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			data = event.getTrackerData()

		if self.clusterflag:
			cand = {}
			noise = {}
			for k,v in data.iteritems():
				c, n = self.process(v) # returns lists
				cand[k] = c
				noise[k] = n
		else:
			cand, noise = self.process(data)

		output = { }
		noiseout = []
		# build the dictionary of clusters for the event storage
		# full unfolding of clusters in dictionaries

		# output holds list of TrackCandidate objects
		# unwrap into wire_lists = list of geiger cylinders
		# noise holds list of geiger cylinders
		if self.clusterflag:
			counter = 1
			for k,v in cand.iteritems():
				for entry in v:
					for trc in entry:
						for gglist in trc.wire_list:
							output[counter] = gglist
							counter += 1
			counter = 1
			for k,v in noise.iteritems():
				for entry in v:
					noiseout.append(entry)
		else:
			for i,entry in enumerate(cand):
				for gglist in entry.wire_list:
					output[i+1] = gglist

			for entry in noise:
				noiseout.append(entry) 

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)

		# store noise in event
		nkey  = self.outkey+'_noise'
		event.setKeyValue(nkey, noiseout)

		# Finally, return the event object
		return event



	def process(self, fulldata):
		if self.flag:
			ll = Legendre_circles(fulldata)
		else:
			ll = Legendre_lines(fulldata)
		ll.run()
		return ll.getTracks(), ll.getNoise()




class NNSetClusterService(object):
	"""
	Pipeline entry for NNSet clustering

	Uses the raw hits 
	which is to be in the form list of hits, i.e. [(x,y,z,Q), ...]

	See the nnsetcluster.py file for details of output format as a 
	dictionary. Should contain line_candidate objects, defined in 
	nnsetcluster.py.
	"""
	def __init__(self, inboxkey, outboxkey, demonstrator_flag=False):
		"""
		Initialise a NNSetClusterService with the parameters given
		inboxkey : the key to a list of tuples as input hits
		            for processing - None
			    means - get the raw hits.

		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.NNSetClusterService')
		self.clusterflag = False # flag cluster data
		self.inkey = inboxkey
		self.outkey = outboxkey
		self.dflag = demonstrator_flag
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
		self.logger.info('Demonstrator flag: %d',self.dflag)
	
	def __repr__(self):
		s = 'NNSet Cluster Service'
		return s


	def __call__(self, event):
		"""
		Process an event with nnsetcluster and output the clustered
		hits to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
				if isinstance(data,dict):
					self.clusterflag = True
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			data = event.getTrackerData()

		if self.clusterflag:
			cand = {}
			noise = {}
			for k,v in data.iteritems():
				c, n = self.process(v) # returns lists
				cand[k] = c
				noise[k] = n
		else:
			if len(data)>0:
				cand, noise = self.process(data)
			else:
				event.setKeyValue(self.outkey, {})
				nkey  = self.outkey+'_noise'
				event.setKeyValue(nkey, [])
				return event

		output = { }
		noiseout = []
		# build the dictionary of clusters for the event storage
		# full unfolding of clusters in dictionaries

		# output dict holds list of dictionaries
		# noise holds list of tracker_hit objects
		if self.clusterflag:
			counter = 1
			for k,v in cand.iteritems():
				for entry in v:
					output[counter] = entry
					counter += 1
			counter = 1
			for k,v in noise.iteritems():
				for entry in v:
					noiseout.append(entry)
		else:
			for i,entry in enumerate(cand):
				output[i+1] = entry

			for entry in noise:
				noiseout.append(entry) 

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)

		# store noise in event
		nkey  = self.outkey+'_noise'
		event.setKeyValue(nkey, noiseout)

		# Finally, return the event object
		return event



	def process(self, fulldata):
		nn = NNSetCluster(fulldata,self.dflag)
		nn.run()
		return nn.getClusters(), nn.getNoise()




class FlattenNN(object):
	"""
	Pipeline entry for flattening a NNSet cluster to plain
	cluster format.

	Uses only NNSetClusterService output. Can't process anything else!
	"""
	def __init__(self, inboxkey, outboxkey):
		"""
		Initialise a FlattenNN with the parameters given
		inboxkey : the key to a list of tuples as input hits
		            for processing - 
			    must get NNcluster output as input
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.FlattenNN')
		self.inkey = inboxkey
		self.outkey = outboxkey
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
	
	def __repr__(self):
		s = 'NNSet flattening Service'
		return s


	def __call__(self, event):
		"""
		Process an  nnsetcluster output and flatten the data
		to plain cluster format to store in outbox key.
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			print 'Error in pipeline: Input must be a NNSet output'
			return None # stop pipeline

		output = { }
		# build the dictionary of clusters for the event storage
		# full unfolding of clusters in dictionaries

		counter = 1
		for k,v in data.iteritems():
			for i,entry in v.iteritems():
				if i>=0:
					for l in entry: # list of track_hits
						if len(l):
							output[counter] = l
							counter += 1
				else:
					if len(entry):
						output[counter] = entry
						counter += 1

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)

		# Finally, return the event object
		return event



class PointCloudService(object):
	"""
	Pipeline entry for preparing tangent point data for clustering

	Uses the raw hits or clustered data.
	"""
	def __init__(self, inboxkey, outboxkey):
		"""
		Initialise a PointCloudService with the parameters given

		inboxkey : the key to a list of tuples as input hits
		            for processing with dbscan - None
			    means - get the raw hits.
			    Must be list of tracker_hit objects.
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.PointCloudService')
		self.inkey = inboxkey
		self.outkey = outboxkey
		self.clusterflag = False
		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
	
	def __repr__(self):
		s = 'PointCloud Service'
		return s


	def __call__(self, event):
		"""
		Process an event with PointCloud and output the 
		tangent points to the outbox key
		"""
		# Get input data according to keystring
		if (self.inkey is not None):
			if (event.hasKey(self.inkey)):
				data = event.getKeyValue(self.inkey)
				if isinstance(data,dict):
					self.clusterflag = True
			else:
				print 'Error in pipeline: key {0} not in event'.format(self.inkey)
				return None
		else:
			data = event.getTrackerData()

		if self.clusterflag:
			cand = {}
			for k,v in data.iteritems():
				c = self.process(v) # returns lists
				cand[k] = c
		else:
			if len(data)>0:
				cand = self.process(data)
			else:
				event.setKeyValue(self.outkey, {})
				return event

		output = {}
		# output dict holds list of dictionaries
		if self.clusterflag:
			counter = 1
			for k,v in cand.iteritems():
				output[counter] = v
				counter += 1
		else:
			for i,entry in enumerate(cand):
				output[i+1] = entry

		# Then set the outboxkey
		event.setKeyValue(self.outkey, output)

		# Finally, return the event object
		return event


	def process(self,data):
		# Now run PointCloud
		pc = PointCloud(data)
		out = pc.run()
		return out





class ImageSegService(object):
	"""
	Pipeline entry for image segmentation clustering

	Uses the raw hits only and can serve as preparation clustering for 
	more detailed clusterers subsequently.
	"""
	def __init__(self, inboxkey, outboxkey):
		"""
		Initialise an ImageSegService with the parameters given

		inboxkey : the key to a list of tuples as input hits
		            for processing with ImageSegmentation - None
			    means - get the raw hits.
			    Must be list of tracker_hit objects.
		outboxkey : output key in the event dictionary
		"""
		self.logger = logging.getLogger('eventloop.ImageSegService')
		self.inkey = inboxkey
		self.outkey = outboxkey

		# write logging info
		self.logger.info('Input key: %s',self.inkey)
		self.logger.info('Output key: %s',self.outkey)
	
	def __repr__(self):
		s = 'Image Segmentation Service'
		return s


	def __call__(self, event):
		"""
		Process an event with Image Segmentation and output the 
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
			cand = self.process(data) # here call the image segmentation algorithm
		else: # empty event
			event.setKeyValue(self.outkey, {})
			return event

		# Then set the outboxkey
		event.setKeyValue(self.outkey, cand)

		# Finally, return the event object
		return event


	def process(self,data):
		# Now run 
		imageseg = ImageSegmentation(data)
		imageseg.run()
		return imageseg.getClusters() # come as dictionary of numbered clusters

