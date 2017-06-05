# control - Pipeline services for pysnemo.reconstruction
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
reconstruction.control - Control for reconstruction
===========================================================

This module provides pipeline services for performing
reconstruction and filtering of entire events, i.e.
pulling together reconstructed tracker data and 
calorimeter data, filtering outputs for the most 
likely model explaining the data and providing therefore
useful output for analysis.
"""

__all__ = ['FitExtrapolatorService', 'ConsolidatorService']

import pysnemo.reconstruction.extrapolator as extra
import pysnemo.reconstruction.alt_extrapolator as aext
from pysnemo.io.edm import Particle
import logging


class FitExtrapolatorService(object):
    """
    Fit Extrapolator Service object - from fitted track candidates
    create intersection ellipses at the foil location and 
    calorimeter location and test for forming a Particle
    object, i.e. finding an associated calohit.
    Input: must find dictionary of fitted track candidates
    
    Output: dictionary - keys of calorimeter hits meta_info
            containing either a list of Particle objects
            or a calohit object for a non-associated calohit.
    """
    def __init__(self, inboxkey, outboxkey, alternative=False):
        """
        Initialise a Fit Extrapolation Service with the parameters given
        inboxkey  : the key to a list of tuples as input hits
                   for processing.
                   alternative switch for switching geometry, default demonstrator 
        outboxkey : output key in the event dictionary

        """
        self.logger = logging.getLogger('eventloop.FitExtrapolatorService')
        self.inkey = inboxkey
        self.outkey = outboxkey
        self.alternative = alternative
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)

        
    def __repr__(self):
        s = 'Multiple Scatter Extrapolator Service in 3D'
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a dictionary 
        of track candidates from fitting
        Output: in key outbox a dictionary is added to event
                containing keys of clusters (int) and a tuple of 
                tuples of extrapolation intersection ellipses
                for foil and calorimter.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                fitdata = event.getKeyValue(self.inkey)
                d = event.getKeyValue('raw')
                if 'calo_hits' in d:
                    calohits = d['calo_hits']
                else:
                    print 'Error in pipeline: calo hits not in event'
                    return None
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None; not permitted.'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...
        out = self.process(fitdata, calohits)

        # Then set the outboxkey
        event.setKeyValue(self.outkey, out)

        # Finally, return the event object
        return event


    def process(self, data, calohits):
        if self.alternative: # catch geometry switch early
            return self.process_alternative(data, calohits)

        # unravel the dictonaries to get to the full fit data
        output = {}
        for calohit in calohits:
            calo_associated = False
            paths = []
            for tup in data: # three fitters in tuple
                results = []
                for cand in tup: # FitResults objects
                    for l,fitresult in cand.fitterpaths.iteritems(): # extract Fit object
                        counter = 0
                        for entry in fitresult:
                            if isinstance(entry,tuple): # is a line fit
                                counter += 1
                        if counter < 2: # not a broken line fit
                            for item in fitresult: # either one or two fits in there
                                if isinstance(item,tuple) and item[2]>=0: # unbroken line fit with valid chi2
                                    #print 'found line tuple ',item
                                    out = extra.extrapolate_line([item],cand.side,calohit)
                                    if len(out): # successful extrapolation
                                        #print 'line extrapolation done.'
                                        particle = Particle((cand.id,l,0)) # tuple as id
                                        particle.set_vertex_foil(out[0])
                                        particle.set_vertex_calo(out[1])
                                        particle.set_calo_hit(calohit)
                                        particle.set_fitter(item)
                                        calo_associated = True
                                        results.append(particle) 
                                elif not isinstance(item,tuple): # is a HelixFit object
                                    out = extra.extrapolate_helix(item,cand.side,calohit)
                                    if len(out): # successful extrapolation
                                        #print 'made a helix extrapolation.'
                                        particle = Particle((cand.id,l,1)) # tuple as id
                                        particle.set_vertex_foil(out[0])
                                        particle.set_vertex_calo(out[1])
                                        particle.set_calo_hit(calohit)
                                        particle.set_fitter(item)
                                        calo_associated = True
                                        results.append(particle)

                        else: # broken line fit
                            grouptuples = []
                            for item in fitresult: # several tuples for broken fit
                                if isinstance(item,tuple): # exclude possible helix fit
                                    grouptuples.append(item) # ignore helix when broken line fitted

                            start = grouptuples[1] # [0] was the initial fit
                            end = grouptuples[-1]

                            if len(start[0])==2: # short fit 2D
                                next = grouptuples[2] # second short fit
                                b = start[0]+next[0] # concatenate
                                e = start[1]+next[1]
                                left = (b,e,start[2],start[3],start[4])
                            else:
                                left = start
                        #print 'left tuple: ',left
                            if len(end[0])==2:
                                previous = grouptuples[-2] # second short fit
                                b = previous[0]+end[0] # concatenate
                                e = previous[1]+end[1]
                                right = (b,e,end[2],end[3],end[4])
                            else:
                                right = end
                        #print 'right tuple: ',right
                            out = extra.extrapolate_line([left,right],cand.side,calohit)
                        #print 'Extrapolate Line result 2: ',out
                            if len(out): # successful extrapolation
                                #print 'made a broken line extrapolation.'
                                particle = Particle((cand.id,l,2)) # tuple as id
                                particle.set_vertex_foil(out[0])
                                particle.set_vertex_calo(out[1])
                                particle.set_calo_hit(calohit)
                                particle.set_fitter(item)
                                particle.set_kink(out[2])
                                particle.set_angles(out[3])                        
                                calo_associated = True
                                results.append(particle)
                paths.append(results)
            if calo_associated: # store listof particles with calohit
                output[calohit.meta_info] = paths
            else: # no particle, only calohit
                output[calohit.meta_info] = [calohit]
        return output




    def process_alternative(self, data, calohits):
        # unravel the dictonaries to get to the full fit data
        output = {}
        for calohit in calohits:
            calo_associated = False
            paths = []
            for cand in data: # FitResults objects
                results = []
                for l,fitresult in cand.fitterpaths.iteritems(): # extract Fit object
                    for item in fitresult: # should be only one, a line fit tuple
                        out = aext.extrapolate_line([item],cand.side,calohit) # alternative extrapolator
                        if len(out): # successful extrapolation
                            particle = Particle((cand.id,l,0)) # tuple as id
                            particle.set_vertex_foil(out[0])
                            particle.set_vertex_calo(out[1])
                            particle.set_calo_hit(calohit)
                            particle.set_fitter(item)
                            calo_associated = True
                            results.append(particle) 
                        
                paths.append(results)
            if calo_associated: # store listof particles with calohit
                output[calohit.meta_info] = paths
            else: # no particle, only calohit
                output[calohit.meta_info] = [calohit]
        return output




class ConsolidatorService(object):
    """
    Consolidator Service collects all reconstruction results and orders
    them ready for the final file output as 'PTD' Particle Track Data
    similar to final Falaise output.
    Input: must find dictionary of extrapolated results and a cut off integer
           stating how many fitted and extrapolated path solutions 
           should maximally be kept, irrespective of the fit being a line
           or a helix. Distinction will merely be the charge measurement.
    
    Output: info collection ready to output to file similar to Falaise
            PTD output in Root file for analysis.
    """
    def __init__(self, inboxkey, outboxkey, cut = 1):
        """
        Initialise a Consolidator Service with keys as below and a cut off integer
        inboxkey  : the key to a dictionary as input for processing.
        outboxkey : output key in the event dictionary

        """
        self.logger = logging.getLogger('eventloop.Chi2FilterService')
        self.inkey = inboxkey
        self.outkey = outboxkey
        self.cut = cut
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)
        self.logger.info('Cut off number: %s',self.cut)

        
    def __repr__(self):
        s = 'Consolidator Service, restricting the number of'
        s+= 'fitted and extrapolated path solutions with associated calorimeter\n'
        s += 'to: '
        s += self.cut
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a dictionary 
        of with lists of lists of Particles for each associated 
        calorimeter hit or a calo_hit object for each isolated calo hit.
        Output: in key outbox a dictionary is added to event
                containing a reduced number, maximally to the cut off,
                of Particle objects in case of associated calo hit
                of calo_hit objects for non-associated calo hits. Can then 
                be written to disk as final reconstruction output.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                data = event.getKeyValue(self.inkey)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None; not permitted.'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...
        out = self.process(data)

        # Then set the outboxkey
        event.setKeyValue(self.outkey, out)

        # Finally, return the event object
        return event
    

    
    def select_particles(self, particles):
        # list of Particles objects for one calorimeter hit.
        sparticles = sorted(particles, key = lambda p: p.chi2)
        return sparticles[:self.cut] # the best N=self.cut in sorted list of chi2



    def clean_results(self, results):
        # make unique path to calo associations 
        # for each available hit cluster depending on chi2
        clids = []
        coll = []
        for entry in results:
            if isinstance(entry,list): # found particles
                for particle in entry:
                    if not (particle.id[0] in clids): # get cluster id if unique
                        coll.append(particle)
                        clids.append(particle.id[0])
                    else:
                        idx = clids.index(particle.id[0]) # where is it if not unique
                        if coll[idx].chi2 > particle.chi2: # replace if better
                            coll.pop(idx)
                            clids.pop(idx)
                            coll.append(particle)
                            clids.append(particle.id[0])
            else:
                coll.append(entry)
                clids.append(-1)
        return coll
                
                

    def process(self, data):
        # unravel the dictonaries to get to the full fit data
        result = []
	for k,calohit in data.iteritems():
            for paths in calohit:
                if isinstance(paths,list): # found particles
                    if len(paths)>0:
                        collection = self.select_particles(paths)
                        #print 'particle selected ',collection[0]
                        #print 'fitter with that particle: ',collection[0].fitter
                        result.append(collection)
                else:
                    result.append(paths) # is a single calo_hit object
        out = self.clean_results(result)
        return out




