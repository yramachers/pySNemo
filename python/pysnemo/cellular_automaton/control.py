# control - Pipeline service for Cellular Automaton
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
control
=======

This module implements the pipeline services for running the
Cellular Automaton through the pySNemo.control mechanism.

"""

__all__ = ['CAService']

from pysnemo.cellular_automaton.configuration import config
from pysnemo.cellular_automaton.clustering import Scaling, Unclustering
from pysnemo.cellular_automaton.ca import CellularAutomaton3D
from pysnemo.cellular_automaton.point import Hit, Cluster
from pysnemo.utility.unique import unique_hit_assignment
from pysnemo.io.edm import tracker_hit
from numpy import asarray, delete, append
from math import sqrt
import logging

class CAService(object):
    def __init__(self, inboxkey, outboxkey, propagation_coordinate=0, theta=30.0, scale_size=3.0, flat=True):
        """
        Initialise the CA service for use in a pipeline.

        The service will then process events one at a time
        as they are passed through the pipeline.

        The inboxkey selection will be taken as the input,
        and output will be to the new key outboxkey:
        
        output in format of dictionary {cluster ID: list of tuples}
            representing the clusters found by the CA

        Some common options can be set in the constructor, others
        must be set through the configure method, which takes
        keyword arguments.

        Parameters permissible in constructor:
        --------------------------------------
        propagation_coordinate : 0, 1, 2
            Coordinate to take as "beam direction"
            Fixed to 0 for SuperNEMO, i.e. x = perpendicular to foil
        theta : float
            Max angle in DEGREES between two segments
            considered part of the same track
            Default: 25.0
        scale_size : float
            Scale size to use when rescaling the voxels
            Default: 3.0
        flat: boolean for switching 2D treatment, True or not, False
        """
        self.logger = logging.getLogger('eventloop.CAService')
        config.setMaxAngleDeg(theta)
        config.setScaleSize(scale_size)
        self.propagation_coordinate = propagation_coordinate
        self.flatbool = flat
        self.inkey = inboxkey
        self.outkey = outboxkey
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)
        self.logger.info('Maximum angle: %f',theta)
    
    def __repr__(self):
        s = 'CA Service'
        return s


    def configure(self, **kwargs):
        """
        Configure the Cellular Automaton service.
        Uses keyword arguments; any arguments passed
        will be set; anything else will retain previous
        or default value.

        Valid keyword arguments:
        ------------------------
        theta : float
            Max angle in DEGREES between two segments
            considered part of the same track

        scale_size : float
            Scale size to use when rescaling the voxels

        cell_radius : float
            Initial search radius to use when generating
            cells

        cell_radius_max : float
            Maximum search radius for cell generation

        cell_radius_increment : float
            Step size between min and max radius for
            cell generation

        min_neighbours : int
            Minimum number of leftward neighbours
            required within search radius

        comparison_tolerance : float
            Tolerance used on floating-point comparisons

        slope_tolerance : float
            Tolerance used on slope comparisons in
            adaptive search for cell generation

        min_track_length : int
            Minimum number of hits in a cluster for it to
            appear in the output list of clusters

        For default values see the pysnemo.cellular_automaton.configuration module.
        Feature masking options are not exported here because masking should be
        done earlier in the pipeline.
        """
        # Map the keyword to the function used to set it
        # so we can iterate through and test kwargs for
        # each value in turn
        arg_map = {
                'theta' : config.setMaxAngleDeg,
                'scale_size' : config.setScaleSize,
                'cell_radius' : config.setCellRadius3D,
                'cell_radius_max' : config.setCellRadius3DMax,
                'cell_radius_increment' : config.setCellRadiusIncrement,
                'min_neighbours' : config.setMinLeftwardNeighbours,
                'comparison_tolerance' : config.setComparisonTolerance,
                'slope_tolerance' : config.setSlopeTolerance,
                'min_track_length' : config.setMinTrackLength
                }
        for keyword, function in arg_map.iteritems():
            if keyword in kwargs:
                function(kwargs[keyword])
    
    def __call__(self, event):
        """
        Process the event given.
        Assumes any charge weighting, feature masking, etc.
        has been done earlier in the pipeline.
        """
        # 0. Get hit selection as list of tuples
        # Get input data according to keystring
        flag = False
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                geiger = event.getKeyValue(self.inkey)
                if isinstance(geiger,list):
                    flag = True
                elif isinstance(geiger,dict):
                    flag = False
                else: 
                    print "Error in pipeline: unknown type to CA Service"
                    return None
            else: 
                print "Error in pipeline: unknown key to CA Service"
                return None
        else:
            geiger = event.getTrackerData()

        # First, get wire coordinates for each gcylinder
        if flag:
            sorted_tracks = self._process(geiger)
            # 7. Store the results
            mtl = config.getMinTrackLength()
            output = {} # dictionary format for output clusters
            for i,entry in enumerate(sorted_tracks):
                if len(entry)>mtl:
                    output[i+1] = entry # numbering from 1 to N cluster
            
        else:
            counter = 1
            output = {} # dictionary format for output clusters
            for k,v in geiger.iteritems():
                sorted_tracks = self._process(v)
                
                # 7. Store the results
                mtl = config.getMinTrackLength()
                for entry in sorted_tracks:
                    if len(entry)>mtl:
                        output[counter] = entry # numbering from 1 to N cluster
                        counter += 1

        event.setKeyValue(self.outkey, output)

        # 8. Return the event object for the next entry in the pipeline
        return event
    

    def _process(self, data):
        hits = self._transform_hits_in(data)
        #print 'in process: transform hits in:',hits
        if len(hits)>0:
            # 1. Scale the hits
            scaled_hits = Scaling(hits).getClusters()
            #print 'in process: scaled hits in:',scaled_hits

            # 2. Perform CA run
            ca = CellularAutomaton3D(scaled_hits, self.propagation_coordinate)
            ca.run()

            # 3. Uncluster CA output
            tracks = self._uncluster(ca.getTracks())
            
            # 4. Transform back to list of plain hits
            output_tracks = self._transform_hits_out(tracks)
            
            # 5. Uniquely assign hits to a cluster
            unique_tracks = unique_hit_assignment(output_tracks)
            
            # 6. Sort by track length, largest to smallest
            sorted_tracks = sorted(unique_tracks, key=lambda x: len(x), reverse=True)
            return sorted_tracks
        else:
            return [] # input data wrong, output will be empty


    ### INTERNAL HELPER FUNCTIONS FOLLOW ###

    def _transform_hits_in(self, hits):
        """
        Convert a list of hits of the form (x, y, z, Q)
        to a list of CA Hit objects
        """
        cluster = Cluster()
        for hit in hits:
            if isinstance(hit,tracker_hit):
                if self.flatbool:
                    wireset = (hit.x,hit.y,0.0,hit.r)
                else:
                    wireset = (hit.x,hit.y,hit.z,hit.r)
                err = [hit.sigmar,hit.sigmaz]
                h = Hit(wireset[:3], wireset[3], err)
                h.setMeta_data(hit.meta_info) # tracker_hit meta data
                cluster.addPoint(h)
            else:
                print 'data type not recognized, not a tracker_hit'
                return []
        return cluster.getPoints()

    def _transform_hits_out(self, tracks):
        """Convert a list of form [list of CA hits, ...] to a list of
        form [list of tracker_hits, ...]
        """
        output_tracks = [ ]
        for track in tracks:
            output = [ ]
            for ca_hit in track:
                coords = ca_hit.getCoordinates()
                err = ca_hit.getErrors()
                plain_hit = tracker_hit(coords[0], coords[1], coords[2], err[1],ca_hit.getRadius(),err[0])
                mi = ca_hit.getMeta_data()
                plain_hit.set_info(mi[0],mi[1],mi[2],mi[3],mi[4],mi[5])
                output.append(plain_hit)
            output_tracks.append(output)
        return output_tracks
    
    def _uncluster(self, tracks):
        """
        Undo the clustering that was performed by the Scaling
        operation
        """
        result = [ ]
        service = Unclustering()
        for track in tracks:
            result.append(service.uncluster(track))
        return result
