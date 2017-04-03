# ca - A cellular automaton
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
ca : Cellular Automaton
=======================

This module implements a cellular automaton for tracking in liquid
Argon. The module incorporates cell generation, as well as forward
and reverse passes of a CA algorithm.

CellGenerator : class
	Generates cells from point-like objects; this should be used
	to make a list of cells for the Automaton

Automaton : class
	Represents a cellular automaton and its state; used to generate
	clusters representing tracks from cells generated using point-like
	objects as input.
"""

__all__ = ['CellularAutomaton']

import math
import operator

from pysnemo.cellular_automaton.configuration import config
from pysnemo.cellular_automaton import utility
from pysnemo.cellular_automaton.cell import Cell, CellSet
from pysnemo.cellular_automaton.point import Cluster
from pysnemo.cellular_automaton.cellgen import CellGenerator3D
from scipy.spatial import KDTree

####                       ####
#### 3D Cellular Automaton ####
####                       ####

class Automaton3D(object):
	"""Class representing a cellular automaton in three dimensions,
	including code to update the state in both forward and reverse
	passes (see documentation and later comments for details).
	"""
	def __init__(self, cells, propagation_coordinate):
		self.cells = cells
		self.propagation_coordinate = propagation_coordinate
		self.forward_step_count = 0
		self.reverse_step_count = 0
		self.neighbourhood = { }
		self.track_candidates = [ ]
		self._buildNeighbourhood()
	
	def _buildNeighbourhood(self):
		"""Build neighbourhood map of cell -> list of leftward neighbours
		"""
		self.neighbourhood = { }
		for cell1 in self.cells:
			for cell2 in self.cells:
				if cell1 is cell2:
					pass
				elif cell2.isLeftwardNeighbourOf(cell1, self.propagation_coordinate):
					# Add neighbour
					if cell1 not in self.neighbourhood:
						self.neighbourhood[cell1] = [ ]
					if cell2 not in self.neighbourhood[cell1]:
						self.neighbourhood[cell1].append(cell2)

	def _doForwardStep(self):
		"""Perform a single step in the forward pass of the algorithm
		"""
		self.forward_step_count += 1
		#print 'Forward step:', self.forward_step_count
		to_update = [ ]
		for cell1, neighbours in self.neighbourhood.iteritems():
			for cell2 in neighbours:
				if cell1.getWeight() == cell2.getWeight():
					#print cell1, cell1.getWeight(), ':::', cell2, cell2.getWeight()
					cell1.setUpdate()
					to_update.append(cell1)
					break
		for cell in to_update:
			cell.step()
		if len(to_update):
			#print 'Cells updated:', len(to_update)
			return True
		else:
			return False
	
	def forwardRun(self):
		run = True
		while run:
			run = self._doForwardStep()
		return self.forward_step_count

	def _doReverseStep(self):
		self.reverse_step_count += 1
		
		# Find cell with max weight
		current_cell = self.remaining_cells.getMaxCell()
		#print 'Highest cell weight:', current_cell.getWeight(), current_cell

		# Append cell to new track candidate
		track = CellSet()
		track.append(current_cell)
		self.track_candidates.append(track)

		# Loop over finding a neighbour until we hit the
		# termination condition
		# N. B. We can use the neighbourhood map here too!
		while current_cell.getWeight() > 1:
			next_weight = current_cell.getWeight() - 1
			neighbours = [n for n in self.neighbourhood[current_cell] if n.getWeight() == next_weight]

			#print 'neighb. list: ', self.neighbourhood[current_cell]
			#print 'nei. weight: ', [n.getWeight() for n in self.neighbourhood[current_cell]]
			#print 'expected weight next: ', next_weight
			if len(neighbours) == 1:
				track.append(neighbours[0])
				current_cell = neighbours[0]
				#print ' Added cell of weight', current_cell.getWeight(), current_cell
			elif len(neighbours) > 1:
				# Find neighbour with smallest angle
				angles = [ ]
				for cell in neighbours:
					angle = cell.angleWith(current_cell)
					angles.append((angle, cell))
				angles = sorted(angles, key=operator.itemgetter(0))
				cell = angles[0][1]
				track.append(cell)
				current_cell = cell
				#print 'neigh > 1:  Added cell of weight', current_cell.getWeight(), current_cell
			else:
				#print 'Stopping. Current cell weight', current_cell.getWeight(), current_cell
				break

		# Now we're out of the loop, reprocess the list of remaining cells
		self.remaining_cells = self.remaining_cells.cellsNotIn(track)

	def reverseRun(self):
		# Build set of cells
		self.remaining_cells = CellSet()
		for cell in self.cells:
			self.remaining_cells.append(cell)

		# Perform reverse steps until no cells remain
		while len(self.remaining_cells):
			self._doReverseStep()
		return self.reverse_step_count

	def getTrackCandidates(self):
		return self.track_candidates

class CellularAutomaton3D(object):
	"""Wrapper around the Automaton3D and CellGenerator3D objects
	to make a set of cells, run an automaton over them, and turn
	the cells back into clusters
	"""
	def __init__(self, hits, propagation_coordinate):
		"""Construct a CellularAutomaton3D object for the given
		hits, running along the specified propagation coordinate.
		"""
		self.hits = hits
		self.propagation_coordinate = propagation_coordinate
		self.tracks = [ ]
	
	def run(self):
		"""Run the cellular automaton to produce a list of track candidates
		"""
		# Generate cells
		cellgen = CellGenerator3D(self.hits, self.propagation_coordinate)
		cells = cellgen.getCells()
		#print 'Number of cells:', len(cells)
		#print cells

		# Run the cellular automaton to obtain tracks
		ca = Automaton3D(cells, self.propagation_coordinate)
		fwd_run_count = ca.forwardRun()
		#print 'Forward run count:', fwd_run_count
		rev_run_count = ca.reverseRun()
		#print 'Reverse run count:', rev_run_count
		ca_tracks = ca.getTrackCandidates()
		# Filter out tracks containing just one cell
		filtered_tracks = [track for track in ca_tracks if len(track) > 1]
		# Turn cells back into clusters
		for track in ca_tracks:
			expanded = Cluster()
			self.tracks.append(expanded)
			for cell in track:
				expanded.addPoints((cell.getFirst(), cell.getLast()))

	def getTracks(self):
		"""Return the list of track candidates
		"""
		return self.tracks

