# __init__ - top level package initialization for pySNemo
#
# Copyright (c) 2013 YR
# Copyright (c) 2013 by The University of Warwick
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
pySNemo: A package for reducing and analysing SuperNEMO data
============================================================

Documentation is provided through the docstrings.

Contents
--------
pySNemo provides the following subpackages:

Subpackages
-----------
::
    ca_glue              --- merge two cellular automaton propagation direction clusters
    cellular_automaton   --- cellular automaton clustering algorithm
    fitter               --- helix and line fitting routines, produces trajectories
    reconstruction       --- final step routines to end at Particle Track Data objects for analysis
    control              --- Control mechanisms for event processing.
    io                   --- Interfaces to persistant data formats.
    utility              --- Widely useful classes and functions.
    cluster              --- clustering for reconstruction
    graphtrack           --- convert raw data rings to hit candidates and produce possible paths


Utility tools
-------------
::

    __version__    --- pySNemo version string

"""

__version__ = '0.1.2'
