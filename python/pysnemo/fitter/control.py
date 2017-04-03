# control - Pipeline services for ROOT Fitter
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
fitter.control
--------------

This module provides pipeline services for performing
track fitting.
"""

__all__ = ['MigradFittingService', 'MS_FittingService2D', 'MS_FittingService3D',
           'Chi2FilterService', 'FitExtrapolatorService']

from ROOT import TGraph2DErrors
import logging
import pysnemo.fitter.pyroot3dfitterV2 as pfit
import pysnemo.fitter.mslinearfitter as msfit
import pysnemo.fitter.ms3dlinefitter as ms3d
from pysnemo.io.edm import FitResults
import numpy as np

class MigradFittingService(object):
    """
    Fitting Service object - fit track candidates
    Input: must find list of ClusterPaths objects
    
    Output: ROOT fitter object

    """
    def __init__(self, inboxkey, outboxkey, Bfield = 25.0):
        """
        Initialise a Fitting Service with the parameters given
        inboxkey  : the key to a list of ClusterPaths objects
        outboxkey : output key in the event dictionary

        """
        self.logger = logging.getLogger('eventloop.FittingService')
        self.inkey = inboxkey
        self.outkey = outboxkey
        self.bfield = Bfield
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)

        
    def __repr__(self):
        s = 'Fitting Service'
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a dictionary 
        of track candidate points with errors.
        Output: in key outbox key a dictionary is added to event
                containing keys of candidate numbers (int) and 
                a list of root virtual fitter result objects.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                data = event.getKeyValue(self.inkey)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...
        cand = self.process(data)

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event



    def process(self, data):
        candidates = []
        for cpaths in data:
            fr = FitResults(cpaths.id)
#            results = {}
            for pnumber, path in cpaths.paths.iteritems():
                hits = []
                errors = []
                for p in path:
                    h = tuple(p[:3]) #(x,y,z) split path tuples
                    err = tuple(p[3:6]) # (ex,ey,ez) split path tuples
                    hits.append(h)
                    errors.append(err)
                # Now construct the data object for fitting
                hist = TGraph2DErrors(len(path))
                hist.SetName("data")
            
                for i,(point,err) in enumerate(zip(hits,errors)):
                    hist.SetPoint(i,point[0],point[1],point[2])
                    hist.SetPointError(i,err[0],err[1],err[2])
                
                    # and fit as line with output for line extrapolator
                if (len(hits)>0):
#                    results[pnumber] = []
                    lf = pfit.linefitter(hist)
                    lfitter = lf.lfitter
                    if (lfitter): # line fit, store a tuple
#                        print lf 
                        fr.add_fitter(pnumber,(lf.linepar,lf.fit_errors,lf.chi2,[],[]))
#                        results[pnumber].append((lf.linepar,lf.fit_errors,lf.chi2,[],[]))
                    else:
                        print "Failed MIGRAD line fit"
                    if self.bfield > 0.0: # want a decent bfield for this, units ??Tesla??
                        hf = pfit.helixfitter(hist,self.bfield)
                        fitresult = hf.helix_result
                        if (fitresult is not None): # store HelixFit object
#                            print fitresult
                            fr.add_fitter(pnumber,fitresult) # added HelixFit
#                            results[pnumber].append(fitresult) 
                        else:
                            print "Failed MIGRAD helix fit"
                hist.Clear()
            candidates.append(fr)
#            candidates[cpaths.id] = results
        return candidates



class MS_FittingService2D(object):
    """
    Fitting Service object - fit track candidates
    Input: must find dictionary of track candidates
    
    Output: dictionary - keys of candidate numbers (int) and 
            a tuple of a best fit list, error list and chisq.

    """
    def __init__(self, inboxkey, outboxkey):
        """
        Initialise a Fitting Service with the parameters given
        inboxkey  : the key to a list of tuples as input hits
                   for processing.
        outboxkey : output key in the event dictionary

        """
        self.logger = logging.getLogger('eventloop.MS_FittingService2D')
        self.inkey = inboxkey
        self.outkey = outboxkey
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)

        
    def __repr__(self):
        s = 'Multiple Scatter Fitting Service in 2D'
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a dictionary 
        of track candidate points with errors.
        Output: in key outbox key a dictionary is added to event
                containing keys of candidate numbers (int) and 
                a tuple of a best fit list, error list and chisq.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                fulldata = event.getKeyValue(self.inkey)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...
        cand = self.process(fulldata)

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event


    def process(self, data):
        cluster = {}
        for cpaths in data:
            results = {}
            for l, path in cpaths.paths.iteritems():
                hits = path

                # Now construct the data object for fitting
                ic, sl = self.start_values(hits)
                #print 'Start values: ic=%f, sl=%f'%(ic,sl)
                bestfit, error, chi, beta_angles, beta_errors = msfit.fitter(hits,ic,sl)
                #print 'Beta Values: ',beta_angles
                #print 'Beta Errors list: ',beta_errors
                if chi>=0: # valid fit
                    bins, angles = msfit.kink_finder(beta_angles, beta_errors)
                    results[l] = []
                    results[l].append((bestfit,error,chi,bins,angles))
                    if len(bins)>0: # split at kinks
                        # now split the data and fit only innermost and 
                        # outermost segment in order to extrapolate better to 
                        # foil and calorimeter
                        bins.sort()
                        bin1 = bins[0]
                        bin2 = bins[-1] # would be identical for single split

                        # split the data
                        if bin1<3:
                            hits1 = hits[:bin1+2] # have at least 3 points
                            #print 'short data1: ',hits1
                        else:
                            hits1 = hits[:bin1]
                            #print 'Hits to left: ',hits1
                        if bin2>len(hits)-3:
                            hits2 = hits[bin2-3:] # have at least 3 points
                            #print 'short data2: ',hits2
                        else:
                            hits2 = hits[bin2:]
                            #print 'Hits to right: ',hits2
                        ic1, sl1 = self.start_values(hits1)
                        #print 'Left part: Start values: ic=%f, sl=%f'%(ic1,sl1)
                        bestfit, error, chi = msfit.segment_fit(hits1,ic1,sl1)
                        results[l].append((bestfit,error,chi,bins,angles))

                        ic2, sl2 = self.start_values(hits2)
                        #print 'Right part: Start values: ic=%f, sl=%f'%(ic2,sl2)
                        bestfit, error, chi = msfit.segment_fit(hits2,ic2,sl2)
                        results[l].append((bestfit,error,chi,bins,angles))
                        
                else: # initial fit failed
                    results[l] = [(bestfit,error,chi,[],[])]
            cluster[cpaths.id] = results
        return cluster


    def start_values(self,hits):
        xtemp=[]
        ytemp=[]
        for hit in hits:
            xtemp.append(hit[0])
            ytemp.append(hit[1])
        xa=np.array(xtemp)
        ya=np.array(ytemp)
        z=np.polyfit(xa,ya,1) # linear quick fit
        sl = z[0]
        ic = z[1]
        return ic, sl



class MS_FittingService3D(object):
    """
    Fitting Service object - fit track candidates
    Input: must find dictionary of track candidates
    
    Output: dictionary - keys of candidate numbers (int) and 
            a tuple of a best fit list, error list and chisq.

    """
    def __init__(self, inboxkey, outboxkey):
        """
        Initialise a Fitting Service with the parameters given
        inboxkey  : the key to a list of tuples as input hits
                   for processing.
        outboxkey : output key in the event dictionary

        """
        self.logger = logging.getLogger('eventloop.MS_FittingService3D')
        self.inkey = inboxkey
        self.outkey = outboxkey
        # write logging info
        self.logger.info('Input key: %s',self.inkey)
        self.logger.info('Output key: %s',self.outkey)

        
    def __repr__(self):
        s = 'Multiple Scatter Fitting Service in 3D'
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a dictionary 
        of track candidate points with errors.
        Output: in key outbox key a dictionary is added to event
                containing keys of candidate numbers (int) and 
                a tuple of a best fit list, error list and chisq.
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                fulldata = event.getKeyValue(self.inkey)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...
        cand = self.process(fulldata)

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event


    def process(self, data):
        cluster = {}
        for cpaths in data:
            results = {}
            for l, path in cpaths.paths.iteritems():
                hits = path

                # Now construct the data object for fitting
                icxy, slxy = self.start_values(hits,0)
                icxz, slxz = self.start_values(hits,1)
                #print 'Initial Start values: icxy=%f, slxy=%f'%(icxy,slxy)
                #print 'Initial Start values: icxz=%f, slxz=%f'%(icxz,slxz)

                stvals = [icxy,slxy,icxz,slxz]
                bestfit, error, chi, beta_angles, beta_errors = ms3d.fitter(hits,stvals)
                #print 'Beta Values: ',beta_angles
                #print 'Beta Errors list: ',beta_errors
                if chi>=0: # valid fit
                    bins, angles = ms3d.kink_finder(beta_angles, beta_errors)
                    #print 'KinkFinder result: ',bins,angles
                    results[l] = []
                    results[l].append((bestfit,error,chi,bins,angles))
                    if len(bins)>0: # split at kinks
                        # now split the data and fit only innermost and 
                        # outermost segment in order to extrapolate better to 
                        # foil and calorimeter
                        bins.sort()
                        bin1 = bins[0]
                        bin2 = bins[-1] # would be identical for single split

                        # split the data
                        if bin1<5:
                            hits1 = hits[:bin1+4] # have at least 5 points
                            #print 'short data1: ',hits1
                            b1, e1, c1, b2, e2, c2 = self.short_fit(hits1)
                            results[l].append((b1,e1,c1,bins,angles))
                            results[l].append((b2,e2,c2,bins,angles))
                        else:
                            hits1 = hits[:bin1]
                            #print 'Hits to left: ',hits1
                            icxy, slxy = self.start_values(hits1,0)
                            icxz, slxz = self.start_values(hits1,1)
                            #print 'Start values: icxy=%f, slxy=%f'%(icxy,slxy)
                            #print 'Start values: icxz=%f, slxz=%f'%(icxz,slxz)
                            stvals = [icxy,slxy,icxz,slxz]
                            bestfit, error, chi, d1,d2 = ms3d.fitter(hits1,stvals)
                            results[l].append((bestfit,error,chi,bins,angles))

                        if bin2>len(hits)-6:
                            hits2 = hits[bin2-5:] # have at least 5 points
                            #print 'short data2: ',hits2
                            b1, e1, c1, b2, e2, c2 = self.short_fit(hits2)
                            results[l].append((b1,e1,c1,bins,angles))
                            results[l].append((b2,e2,c2,bins,angles))
                        else:
                            hits2 = hits[bin2:]
                            #print 'Hits to right: ',hits2
                            icxy, slxy = self.start_values(hits2,0)
                            icxz, slxz = self.start_values(hits2,1)
                            #print 'Start values: icxy=%f, slxy=%f'%(icxy,slxy)
                            #print 'Start values: icxz=%f, slxz=%f'%(icxz,slxz)
                            stvals = [icxy,slxy,icxz,slxz]
                            bestfit, error, chi, d1,d2 = ms3d.fitter(hits2,stvals)
                            results[l].append((bestfit,error,chi,bins,angles))

                else: # initial fit failed
                    results[l] = [(bestfit,error,chi,[],[])]
            cluster[cpaths.id] = results
        return cluster


    def start_values(self,hits,which):
        if which==0: # do the xy projection
            xtemp=[]
            ytemp=[]
            for hit in hits:
                xtemp.append(hit[0])
                ytemp.append(hit[1])
            xa=np.array(xtemp)
            ya=np.array(ytemp)
            z=np.polyfit(xa,ya,1) # linear quick fit
            sl = z[0]
            ic = z[1]
        else:   # do the xz projection
            xtemp=[]
            ytemp=[]
            for hit in hits:
                xtemp.append(hit[0])
                ytemp.append(hit[2])
            xa=np.array(xtemp)
            ya=np.array(ytemp)
            z=np.polyfit(xa,ya,1) # linear quick fit
            sl = z[0]
            ic = z[1]
        return ic, sl


    def short_fit(self,hits):
        datax = []
        datay = []
        dataz = []
        errx = []
        erry = []
        errz = []
        for d in hits:
            datax.append(d[0])
            datay.append(d[1])
            dataz.append(d[2])
            errx.append(d[3])
            erry.append(d[4])
            errz.append(d[5])
        icxy, slxy = self.start_values(hits,0)
        icxz, slxz = self.start_values(hits,1)
        #print 'short fit, Start values: icxy=%f, slxy=%f'%(icxy,slxy)
        #print 'short fit, Start values: icxz=%f, slxz=%f'%(icxz,slxz)
        bestfit1, error1, chi1 = ms3d.segment_fit2D(datax,datay,errx,erry,icxy,slxy)
        bestfit2, error2, chi2 = ms3d.segment_fit2D(datax,dataz,errx,errz,icxz,slxz)
        return bestfit1, error1, chi1, bestfit2, error2, chi2




