# pipeline - Pipelines for event processing
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
pysnemo.control.pipeline
=======================

Provides classes which can be used to chain
operations together in order to process a
single event.

These operations can do just about anything,
but each must be a Python callable which takes
a single event object as its argument.
"""

__all__ = ['Pipeline']
import logging

class Pipeline(object):
	def __init__(self, *operations):
		"""
		Construct a Pipeline from one or more operations.
		Each will be run in order.

		An operation is a Python callable which takes a
		single argument, an event-like object.
		
		The operation must return either an event-like
		object (which may be the original event it was
		passed), or None.

		If the operation returns None it is assumed that it
		does not wrap or unwrap the event, and the previous
		event object is used for the next operation.

		If it returns anything else, that value is used for
		the next operation's event argument; so it had better
		be an event-like object!
		"""
		self.logger = logging.getLogger('eventloop.Pipeline')
		self.operations = operations
		self.logger.info('Pipeline : %s',self.operations)
	
	def __call__(self, event):
		"""
		Function call operator. Apply each operation in turn.
		"""
		for op in self.operations:
			event_ = op(event)
			if event_ is not None:
				event = event_

