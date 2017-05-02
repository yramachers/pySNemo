# clustering - Module for clustering hits into larger groups
#
# Copyright (c) 2013 YR
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
glueing
-------

This module provides services for clustering hits together by merging
smaller clusters to the main cluster candidates:

"""
from numpy import array, where
from scipy.spatial import KDTree
from pysnemo.io.edm import tracker_hit

__all__ = ['MergeXY','RoadSweeper']

class MergeXY(object):
	'''
	Merging of largest clusters. One from CA x propagation, one
	from CA y propagation.
	'''
	def __init__(self, data1, data2):
		"""Initialise a MergeXY object with two(!) 
		dictionaries of clusters of wire coordinates as input. One 
		originating from the CA propagation in x direction, 
		the other from propagation in y direction.
		"""
		self.clsx = data1
		self.clsy = data2
		self.out = {}
		self._merge()


	def getClusters(self):
		return self.out
	


	def _remove_duplicates(self,temp):
		metacl = {}
		for k, v in temp.iteritems():
			vshort = []
			for entry in v:
				vshort.append(entry.meta_info) # meta info only
			metacl[k] = vshort
		remove = []
		ks = metacl.keys()
		for i in range(len(ks)): # for all clusters
			s1 = set(metacl[ks[i]])
			for j in range(i+1,len(ks)):
				s2 = set(metacl[ks[j]])
				sd1 = s1.difference(s2)
				sd2 = s2.difference(s1)
				if len(sd1)<1: # no diff
					# s1 is a sub-set
					remove.append(ks[i]) # [0]=key
				elif len(sd2)<1: # no diff
					# s2 is a sub-set
					remove.append(ks[j]) # [0]=key
		#print 'remove has: ',remove
		temp_out = {}
		counter = 1
		for k, v in temp.iteritems():
			if k not in remove:
				temp_out[counter] = v
				counter += 1
		for k, v in temp_out.iteritems():
			self.out[k] = v

				

	def _maxkeys(self,data,l):
		if len(data)>0:
			key1 = data.keys()[0]
			val1 = data.values()[0]
			nmax = len(val1)
			l.append(key1)
			for i in range(1,len(data)):
				k = data.keys()[i]
				v = data.values()[i]
				if len(v)>=nmax-1:
					l.append(k)
		return l


	def _merge(self):
		temp = {}
		klist = []
		#print 'MergeXY: clusters in x: ',len(self.clsx)
		#print 'clusters in y: ',len(self.clsy)
		klist = self._maxkeys(self.clsx,klist)
		#print 'clusters in x max keylist: ',klist

		if len(klist)>0:
			for i,entry in enumerate(klist,start=1):
				temp[i] = self.clsx[entry]
		klist = []
		klist = self._maxkeys(self.clsy,klist)
		#print 'clusters in y max keylist: ',klist

		if len(klist)>0:
			for j,entry in enumerate(klist,start=i+1):
				temp[j] = self.clsy[entry]
		self._remove_duplicates(temp)



class RoadSweeperOld(object):
	'''
	Road following along clusters and collecting left-over wire hits,
	within one wire distance to clusters.
	'''
	def __init__(self, data, raw_data):
		"""Initialise a RoadSweeper object with  
		dictionary of clusters of wire coordinates as input as well
		as the full raw data of wires for sweeping. Deltaz is the 
		max tolerance to count as a neighbour in z.
		"""
		self.cls = data
		self.ggdata = raw_data
		self.wirecoords = []
		for gg in self.ggdata:
			if isinstance(gg,tracker_hit):
				self.wirecoords.append((gg.x,gg.y,gg.z))
			elif isinstance(gg.tangentpoint_hit):
				self.wirecoords.append((gg.x,gg.y,gg.z))
		#print 'wire coords: ',self.wirecoords
		self.out = {}
		if len(data)<1:
			print 'No data to sweep in RoadSweeper!'
		else:
			self._sweep()



	def find_geiger(self, metainfo):
		for entry in self.ggdata:
			# check id's in place [0] in meta info
			if entry.meta_info[0] == metainfo[0]:
				gg = entry
		return gg


	def getClusters(self):
		# translate into gcylinder objects
		geiger = {}
		if len(self.out)<1:
			return {},[]
		for k,v in self.out.iteritems():
			gglist = []
			for data in v:
				meta = data.meta_info #meta info
				gglist.append(self.find_geiger(meta))
			geiger[k] = gglist
		noise = []
		idlist = []
		for k,gglist in geiger.iteritems():
			for gg in gglist:
				idlist.append(gg.meta_info[0]) # the id number
		nidlist = self._complement(idlist,len(self.ggdata))
		if len(nidlist):
			for id in nidlist:
				noise.append(self.find_geiger([id]))
		return geiger, noise
	

	def _complement(self,l,length):
		A = set(l)
		Bl = [i for i in range(length)]
		B = set(Bl)
		return list(B.difference(A))


	def _pair_check(self, wire):
 		#print 'in pair check wire is: ',wire.x,wire.y
		warr = array(self.wirecoords)
		# check on next in x 
		indxup = where(abs(wire.x-warr[:,0])<46.0)[0]
		sx = set(indxup)
 		#print 'in pair check sx: ',sx
		# check on next in y
		indyup = where(abs(wire.y-warr[:,1])<46.0)[0]
		sy = set(indyup)
 		#print 'in pair check sy: ',sy
		# check on next in z
		#indzup = where(abs(wire.z-warr[:,2])<self.deltaz)[0]
		#sz = set(indzup)
		# check on all conditions, logical 'and'
		#ind = sxu.intersection(syu.intersection(szu)) 
		ind = sx.intersection(sy)
		#print 'set intersection: ',ind
		return ind
	

	def _isinside(self,meta,data):
		tempmeta = []
		for entry in data:
			tempmeta.append(entry.meta_info)
		return (meta in tempmeta)


	def _remove_duplicates(self,temp):
		out = {}
		for k,v in temp.iteritems():
			temp_out = []
			remove = []
			for i,entry in enumerate(v):
				meta = entry.meta_info
				vrest = v[i+1:]
				if self._isinside(meta,vrest):
					remove.append(i)
			for i,entry in enumerate(v):
				if i not in remove:
					temp_out.append(entry)
			out[k] = temp_out
		return out



	def _remove_duplicate_clusters(self,temp):
		metacl = {}
		for k, v in temp.iteritems():
			vshort = []
			for entry in v:
				vshort.append(entry.meta_info) # meta info only
			metacl[k] = vshort
		remove = []
		ks = metacl.keys()
		for i in range(len(ks)): # for all clusters
			s1 = set(metacl[ks[i]])
			for j in range(i+1,len(ks)):
				s2 = set(metacl[ks[j]])
				sd1 = s1.difference(s2)
				sd2 = s2.difference(s1)
				if len(sd1)<2: # diff of 1 wire max
					# s1 is a sub-set
					remove.append(ks[i]) # [0]=key
				elif len(sd2)<2: # diff of 1 wire max
					# s2 is a sub-set
					remove.append(ks[j]) # [0]=key
		#print 'sweep remove has: ',remove
		temp_out = {}
		counter = 1
		for k, v in temp.iteritems():
			if k not in remove:
				temp_out[counter] = v
				counter += 1
		for k, v in temp_out.iteritems():
			self.out[k] = v

				
	def _sweep(self):
		temp = {}
		#print 'SWEEP: got %d clusters.'%len(self.cls)
		for k,cl in self.cls.iteritems():
			indlist = []
			for entry in cl:
				ind = self._pair_check(entry)
				if len(ind)>0:
					indlist.append(list(ind))
			for entry in indlist:
				for i in entry:
					cl.append(self.ggdata[i])
			temp[k] = cl

		out = self._remove_duplicates(temp)
		self._remove_duplicate_clusters(out) # sets self.out





class RoadSweeper(object):
	'''
	Road following along clusters and collecting left-over wire hits,
	within one wire distance to clusters.
	'''
	def __init__(self, data, raw_data):
		"""Initialise a RoadSweeper object with  
		dictionary of clusters of wire coordinates as input as well
		as the full raw data of wires for sweeping. Deltaz is the 
		max tolerance to count as a neighbour in z.
		"""
		self.cls = data
		self.ggdata = raw_data
		self.ggdatami = [trh.meta_info for trh in raw_data]
		if len(data)<1:
			print 'No data to sweep in RoadSweeper!'
		else:
			self._sweep()



	def getClusters(self):
		return self.out, self.noise
	

	def _remove_duplicate_clusters(self,temp):
		metacl = {}
		for k, v in temp.iteritems():
			vshort = [trh.meta_info for trh in v]
			metacl[k] = vshort
		remove = []
		ks = metacl.keys()
		for i in range(len(ks)): # for all clusters
			s1 = set(metacl[ks[i]])
			for j in range(i+1,len(ks)): # test all others
				s2 = set(metacl[ks[j]])
				sd1 = s1.difference(s2)
				sd2 = s2.difference(s1)
				if len(sd1)<2: # diff of 1 wire max
					# s1 is a sub-set
					remove.append(ks[i]) # [0]=key
				elif len(sd2)<2: # diff of 1 wire max
					# s2 is a sub-set
					remove.append(ks[j]) # [0]=key
		#print 'duplicate remove remove has: ',remove
		temp_out = {}
		counter = 1
		for k, v in temp.iteritems():
			if k not in remove:
				temp_out[counter] = v
				counter += 1
		return temp_out

				
	def _sweep(self):
		# find all  non clustered hits
		self.out = { }
		self.noise = []
		temp = { }
		miset = set()
		indexset = set()
		for mi in self.ggdatami: # check meta info
			miset.add(self.ggdatami.index(mi)) # all hit indices
			for k,cl in self.cls.iteritems():
				clsmi = [trh.meta_info for trh in cl]
				if mi in clsmi:
					indexset.add(self.ggdatami.index(mi))
		noiseidx = miset.difference(indexset) # hits not in a cluster
		
		swept = []
		for k,cl in self.cls.iteritems():
			side = cl[0].meta_info[3] # same for all cluster members
			clsmi = [trh.meta_info[4:] for trh in cl] # row, col
			kdt = KDTree(clsmi)
			indexlist = list(indexset)
			for idx in noiseidx:
				nmeta =  self.ggdatami[idx]
				if side == nmeta[3]: # only for same side noise
					pairs = kdt.query_ball_point(nmeta[4:],1.5)
					if len(pairs)>0: # found a nearest neighbour
						swept.append(idx) # to remove from noise
						cl.append(self.ggdata[idx]) # sweep up
		for k,cl in self.cls.iteritems():
			temp[k] = cl # output with refreshed cls
		for idx in noiseidx:
			if idx not in swept:
				self.noise.append(self.ggdata[idx])
		self.out = self._remove_duplicate_clusters(temp)

