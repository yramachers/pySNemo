# unique - Uniquely allocate hits to one cluster
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
pysnemo.utility.unique
====================

Utility to allocate hits uniquely to the longest track in which
they appear.
"""

__all__ = ['unique_hit_assignment']

def unique_hit_assignment(tracks):
    '''Perform unique allocation of hits to the
    longest track in which they appear.
    This is done by sorting the tracks by length,
    then filtering the hits in the shortest, leaving only
    those which are not in any other track.
    '''
    sorted_tracks = sorted(tracks, key=lambda x: len(x),
            reverse=True) # Reversed, since pop() operates on the end of a list
    new_tracks = [ ]
    while len(sorted_tracks):
        current_track = sorted_tracks.pop()
        new_track = [ ]
        for hit in current_track:
            if not _hit_in_one_of(hit, sorted_tracks):
                # Have unique hit
                new_track.append(hit)
        if len(new_track):
            # Have something in this track
            new_tracks.append(new_track)
    return new_tracks

def _hit_in_one_of(hit, tracks):
    '''Return True if the hit is in one
    of the tracks given, otherwise return
    False.
    '''
    for track in tracks:
        if hit in track:
            return True
    return False

