# control - Pipeline services for pysnemo.graphtrack
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
graphtrack.control - Control for tracking
===========================================================

This module provides pipeline services for performing
track finding using the graph algorithm on ring data

"""

__all__ = ['TrackingService']

from pysnemo.graphtrack.graphpathfinder import cleangraph
from pysnemo.utility.rawgeigerhits import raw_hits
import logging

class TrackingService(object):
    """
    Tracking Service object - running the tracking reconstruction!
    Input: lattice constant fixed default at 4.4 cm
    """
    def __init__(self, inboxkey, outboxkey, sigma=1.0, lattice=44.0):
        """
        Initialise a Tracking Service with the parameters given
        inboxkey  : the key to a list of tuples as input hits
                   for processing.
                   Input must be a dictionaries !
        outboxkey : output key in the event dictionary
        lattice   : wire grid lattice constant - geometry default = 4.4cm
        """
        self.logger = logging.getLogger('eventloop.TrackingService')
        self.initsigma = sigma
        self.lattice = lattice

        self.inkey = inboxkey
        self.outkey = outboxkey
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)
        self.logger.info('tolerance on tangent points: %f',self.initsigma)
        self.logger.info('Lattice constant: %f',self.lattice)

        
    def __repr__(self):
        s = 'Tracking Service'
        return s


    def __call__(self, event):
        """
        process the data from file, find all tangent points and 
        finish with the list of shortest paths + path lengths, 
        which are by definition the track candidates.

        Output: in key outbox a dictionary is added to event
                containing keys id cluster candidates
                each holding a dictionary of track candidates
                with lists of 
                tuples of path points (x,y,z,dx,dy,dz) and their errors.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                data = event.getKeyValue(self.inkey) # must be dict
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in Tracking Service: Will not work on raw data.'
            return None
        
        
        # number crunshing get candidates as 
        # numbered (key) list of tuples of points and errors in a dictionary
        cand = {}
        for k,v in data.iteritems(): # each cluster gets a
            if len(v)>3: # have at least 4 hits
                cand[k] = self.get_paths(v) # dict of paths

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event



    def duplicatecheck(self,path):
        tol = 1.e-6
        ind = []
        for count,tup1 in enumerate(path):
            for m in xrange(count+1,len(path)):
                if abs(tup1[0]-path[m][0])<tol and abs(tup1[1]-path[m][1])<tol:
                    ind.append(m)
        if len(ind)>0:
            newpath = []
            for n,p in enumerate(path):
                if n not in ind:
                    newpath.append(p)
            return newpath
        else:
            return path
        

    def tracker_side_check(self,path):
        sidelist = []
        for p in path:
            mi = p[-1]
            sidelist.append(mi[3]) # side info
        sideset = set(sidelist) # removes all multiple entries
        return list(sideset) 
        

    def get_paths(self, gdata):
        '''
        Using raw_hits object from utility.
        Provides: make_hitdictionary() method with
        Output:   dict of actual hits for each ring: format is
        key = id number taken from input list sequences
        value is a list of triplets: 2 cylinder objects to make a pair
        and one resulting tp_controller object holding the 
        actual filtered hits on both rings in that pair (and errors)
        '''
        self.sigma = self.initsigma # reset running self.sigma
        candidates = {}
        while (len(candidates)<1 and self.sigma<7): # give up at sigma=7
#            print 'sigma prefactor try: ',self.sigma

            rh = raw_hits(self.sigma, self.lattice, gdata)
            hitsdict = rh.make_hitdictionary()
            
            if (len(hitsdict)>0):
                
                start, end = rh.start_end_points()
                #print 'list of start and end points: '
                #print start
                #print end
                # find out connections between any two rings
                gg=cleangraph(hitsdict,self.sigma) # sigma tolerance
                nd=gg.graph.nodes()
                #if self.sigma==self.initsigma:
                    #print 'all nodes: ',nd
                # find nearest match in nodes to start,end points
                # generous check
                startcopy=[]
                ind=[]
                for t in start:
                    tx=t[0]
                    ty=t[1]
                    for m in xrange(0,len(nd)):
                        if ((abs(nd[m][0]-tx)<6.3) and (abs(nd[m][1]-ty)<6.3)):
                            if (m not in ind):
                                ind.append(m)
                for entry in ind:
                    startcopy.append(nd[entry])
                #print "Start:"
                #print startcopy
                
                endcopy=[]
                ind=[]
                for t in end:
                    tx=t[0] #x, y comparison sufficient
                    ty=t[1]
                    for m in xrange(0,len(nd)):
                        if ((abs(nd[m][0]-tx)<6.3) and (abs(nd[m][1]-ty)<6.3)):
                            if (m not in ind):
                                ind.append(m)
                for entry in ind:
                    endcopy.append(nd[entry])
                #print "Ends:"
                #print endcopy
                
                pl = []
                paths = []
                for t in startcopy: # tuples in the list are the starting nodes
                    for te in endcopy:
                        p,l = gg.find_path(t,te)
                        paths.append(p)
                        pl.append(l)
                        #print 'found path: ',p
            
                ll = [l for l in pl if l>0] # collect all valid paths
#                print "Graph tracking: list of pathlengths:"
#                print ll
                counter = 1
                for p in paths:
                    cleanpath = self.duplicatecheck(p)
                    if len(cleanpath):
                        side = self.tracker_side_check(cleanpath)
                        if len(side)>1:
                            continue # next, path with mixed sides
                        if side[0]==0: # fix key sign for one tracker side
                            counter = - abs(counter)
                        else:  # opposite sign for opposite tracker side
                            counter = abs(counter)
                        candidates[counter] = cleanpath#(x,y,z,ex,ey,ez,ids,mi)
                        counter = abs(counter) + 1 # always increase no matter what sign

            self.sigma += 1
        return candidates
