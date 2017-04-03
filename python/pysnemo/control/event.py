# event - Event object and event loop
# Feature update Aug. 2012 YR
#
# Copyright (c) 2012 Andrew J. Bennieston <A.J.Bennieston@warwick.ac.uk>
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
control.event
=======================

Event: An event object which acts
as a wrapper to the dictionary generated by
the file readers. As reconstruction proceeds
more keys are added by processing modules.
These are requested by modules as input data
by a compulsory keyword string argument in 
all constructors where input data is needed.


EventLoop: An object which handles reading data from
a file and processing each event in turn. EventLoop
instances can apply a set of global cuts, and process
a pipeline of services to reconstruct an event.

"""

__all__ = ['Event', 'EventLoop'] # Note: feature change to 'no event wrappers'

import pysnemo.io.decorators
import logging

class Event(object):
	def __init__(self, filename, event_id, event_dict):
		"""
		Create an Event object from a dictionary
		created by a file reader

		"""
		self.filename = filename
		self.event_id = event_id
		self.event_dict = event_dict
		self.tracker_data = None
		self.calo_hits = None

	def __str__(self):
		'''
		Event object should know how to print itself
		'''
		s = 'Event from file: %s\n' % self.filename
		s += 'with event ID = %d\n' % self.event_id
		return s

	def getFilename(self):
		return self.filename

	def getEventID(self):
		return self.event_id
	
	def getTrackerData(self): # interface change to tracker data
		"""
		This method returns the tracker data from file
		(a list of hits)

		This is a list: [gcylinder, ...]
		of geiger cylinder objects, see edm.py. By default, it is
		all the hits in the event.
		"""
		if self.tracker_data is None:
			if 'gg' in self.event_dict['raw']:
				# Build the tracker data
				self.tracker_data = self.event_dict['raw']['gg']
			else:
				self.tracker_data = []
		return self.tracker_data


	def getCaloHits(self): # interface change to tracker data
		"""
		This method returns the calo hits from file
		(a list of hits)

		This is a list: [calo_hit, ...]
		of calo_hit objects, see edm.py. By default, it is
		all the hits in the event.
		"""
		if 'calo_hits' in self.event_dict['raw']:
			# Build the tracker data
			self.calo_hits = self.event_dict['raw']['calo_hits']
			return self.calo_hits
		else:
			print "NOT a SNG4 file - no calo hits"
			return []


	def getTrueParticle(self): # interface to truth data
		'''
		list of trueparticle objects, see edm.py
		'''
		return self.event_dict['raw']['trueparticle']


	def getTrueVertex(self): # interface to truth data
		'''
		list of truevertex objects, see edm.py
		'''
		if 'vertex' in self.event_dict['raw']:
			return self.event_dict['raw']['vertex']
		else:
			print "NOT a SNG4 file - vertex"
			return []
		


	def getKeys(self):
		"""Return a list of keys provided by this event.

		"""
		return self.event_dict.keys()

	def hasKey(self, key):
		"""
		Returns True if the key exists, False otherwise.

		"""
		if key in self.event_dict:
			return True
		else:
			return False
	
	def setKeyValue(self, key, value):
		"""
		Sets the key to the value specified.
		Should be the main module processed data storage method.

		"""
		self.event_dict[key] = value
	
	def addToKeyValue(self, key, value):
		"""Adds the value given to the value currently in the key.
		(!)Requires that += works for the relevant type(s).
		If the key is not in the dictionary, it is created and
		initialised to the value passed as an argument.

		"""
		if self.hasKey(key):
			self.event_dict[key] += value
		else:
			self.event_dict[key] = value

	def getKeyValue(self, key):
		"""
		Gets the value associated with the key given.
		Should be the main event data access method for modules.

		"""
		if self.hasKey(key):
			return self.event_dict[key]
		else:
			return None
	


class EventLoop(object):
	def __init__(self, filename, operation=None, global_cuts=None, failed_operation=None, first_event=None, last_event=None):
		"""Construct an EventLoop object 
		
		Parameters:
		  filename		   : file to read events from
		  - .root  : SNG4 root file
		  - .tsim  : toysimulation file
		  - .evt   : processed event file

		  operation		   : Python callable which takes an
		  event-like object as an argument

		  global_cuts	   : List of boolean predicates which
		  take event-like object arguments

		  failed_operation : Operation to perform on events
		  that fail the global cuts
		  (default: do nothing)

		  first_event      : Event number to start processing at

		  last_event       : Event number to end processing at (inclusive)
	
		The run() method runs operation for each event from
		filename for which all global_cuts return true.
		"""
		self.logger = logging.getLogger('eventloop.EventLoop')
		# Set member variables
		self.filename = filename
		self.operation = operation
		self.global_cuts = global_cuts
		self.failed_operation = failed_operation
		self.first_event = first_event
		self.last_event = last_event
		
		# Open file reader and determine number of events
		if isinstance(self.filename, list):
			self.reader = pysnemo.io.decorators.create_reader_instance(self.filename)
			self.num_events = self.reader[1].number_of_events() # processed file numbers
		else:
			self.reader = pysnemo.io.decorators.create_reader_instance(self.filename)
			self.num_events = self.reader.number_of_events()

		if self.first_event is None:
			self.first_event = 0
		if self.last_event is None:
			self.last_event = self.num_events - 1
		if self.last_event >= self.num_events:
			self.last_event = self.num_events - 1
		if (self.first_event >= self.num_events or self.first_event>self.last_event):
			self.first_event = 0
		# write in log file
		self.logger.info('Input: %s',self.filename)
		self.logger.info('From event %d to event %d.', self.first_event,self.last_event)

	
	def getNumEvents(self):
		"""
		Return the number of events contained in the Lamu file
		"""
		return self.num_events

	def run(self):
		"""Run over each event in turn.
		1. Read event from Lamu file
		2. Apply global cuts
		3. Perform operation if all global cuts return True

		operation must produce any output required, so if it is a
		pipeline, the last operation must produce output, otherwise
		all will be for nothing!
		"""
		# Keep some statistics
		num_processed = 0
		num_failed = 0
		for event_number in xrange(self.first_event, self.last_event + 1):
			if isinstance(self.reader, list):
				evt = self.reader[0].get_event(event_number) # raw key
				dummy = self.reader[1].get_event(event_number) # processed
				evt[dummy.keys()[0]] = dummy.values()[0] # add processed key to Event dictionary
			else:
				evt = self.reader.get_event(event_number)
			current_event = Event(self.filename, event_number, evt)
			# if cuts exist, apply them
			if self.global_cuts is not None:
				return_values = [cut(current_event) for cut in self.global_cuts]
			else:
				return_values = [True]
			if False in return_values:
				num_failed += 1
				if self.failed_operation is not None:
					self.failed_operation(current_event)
				continue
			else:
				# Now we've passed all the cuts, process this event
				num_processed += 1
				self.operation(current_event)
		# Return the statistics
		return (num_processed, num_failed)

