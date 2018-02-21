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

__all__ = ['FittingService', 'MigradFittingService', 'MS_FittingService2D', 'MS_FittingService3D']

from ROOT import TGraph2DErrors
import logging
import pysnemo.fitter.pyroot3dfitterV2 as pfit
import pysnemo.fitter.mslinearfitter as msfit
import pysnemo.fitter.ms3dlinefitter as ms3d
from pysnemo.io.edm import FitResults, ClusterPaths
import numpy as np



class FittingService(object):
    """
    Fitting Service object - fit track candidates
    Attempt several different fit models in decreasing degree of difficulty:
    From helix to broken line to simple line. Results if present are handed back

    Input: must find list of ClusterPaths objects or a dictionary from 
           direct processing of the path finding module.
    
    Output: Tuple of lists of FitResults objects (see EDM)

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
        process the data from event, which must find a list of ClusterPaths 
        objects.
        Output: Tuple of lists of FitResults objects
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                data = event.getKeyValue(self.inkey)
                if isinstance(data,list): # list of clusterpath objects from reading file
                    cand = self.process(data)
                elif isinstance(data,dict): # direct tracker output
                    clpdata = self.transfer(data)
                    cand = self.process(clpdata)
                else:
                    print 'Error in pipeline: unknown input type ',type(data)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event



    def transfer(self, datain):
        clpdata = []
        for k,val in datain.iteritems():
            clp = ClusterPaths(k)
            for pid, entry in val.iteritems():
                clp.add_path(pid,entry)
            clpdata.append(clp)
        return clpdata



    def process(self, data):
        '''
        Use the separate fitting services stand-alone, separately from
        event processing.
        '''
        if self.bfield>0:
            migradhelix = MigradFittingService("in","out",self.bfield) # inbox and outbox keys 
            # Helix first
            helixfits = migradhelix.process(data)
        else:
            helixfits = None

        migradline = MigradFittingService("in","out",0.0)
        brline = MS_FittingService3D("in2","out2") # are irrelevant for internal use like here

        # Broken lines next
        blfits = brline.process(data)

        # Line next
        linefits = migradline.process(data)      

        if helixfits is not None:
            return (helixfits, blfits, linefits)
        else:
            return (blfits, linefits)


class MigradFittingService(object):
    """
    Fitting Service object - fit track candidates
    Input: must find list of ClusterPaths objects
    
    Output: List of FitResults objects

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
        s = 'Migrad Fitting Service'
        return s


    def __call__(self, event):
        """
        process the data from event, which must be a list of 
        ClusterPath objects.
        Output: List of FitResults objects
        """
        # Get input data according to keystring
        if (self.inkey is not None):
            if (event.hasKey(self.inkey)):
                data = event.getKeyValue(self.inkey)
                if isinstance(data,list): # list of clusterpath objects from reading file
                    cand = self.process(data)
                elif isinstance(data,dict): # direct tracker output
                    clpdata = self.transfer(data)
                    cand = self.process(clpdata)
                else:
                    print 'Error in pipeline: unknown input type ',type(data)
            else:
                print 'Error in pipeline: key %s not in event'%self.inkey
                return None
        else:
            print 'Error in pipeline: key is None'
            return None

        
        # number crunshing get candidates as 
        # numbered (key) list of ...

        # Then set the outboxkey
        event.setKeyValue(self.outkey, cand)

        # Finally, return the event object
        return event



    def transfer(self, datain):
        clpdata = []
        for k,val in datain.iteritems():
            clp = ClusterPaths(k)
            for pid, entry in val.iteritems():
                clp.add_path(pid,entry)
            clpdata.append(clp)
        return clpdata



    def process(self, data):
        candidates = []
        for cpaths in data:
            fr = FitResults(cpaths.id)
            for pnumber, path in cpaths.paths.iteritems():
                hits = []
                errors = []
                for p in path:
                    h = tuple(p[:3]) # (x,y,z) split path tuples
                    err = tuple(p[3:6]) # (ex,ey,ez) split path tuples
                    hits.append(h)
                    errors.append(err)
                # Now construct the data object for fitting
                hist = TGraph2DErrors(len(path))
                hist.SetName("data")
                # print 'path number %d'%pnumber
                # and fit as line with output for line extrapolator
                if (len(hits)>0):
                    if self.bfield <= 0.0: # no bfield, line fit
                        for i,(point,err) in enumerate(zip(hits,errors)): 
                            hist.SetPoint(i,point[0]*1.0e-3,point[1]*1.0e-3,point[2]*1.0e-3) # in [m]
                            hist.SetPointError(i,err[0]*1.0e-3,err[1]*1.0e-3,err[2]*1.0e-3)
                            # print 'data in: ',point
                        lf = pfit.linefitter(hist)
                        lfitter = lf.lfitter
                        if (lfitter): # line fit, store a tuple
                            fr.add_fitter(pnumber,(lf.linepar,lf.fit_errors,lf.chi2,[],[]))
                        #else:
                            #print "Failed MIGRAD line fit"
                    elif self.bfield > 0.0: # want a decent bfield for this, units Tesla
                        for i,(point,err) in enumerate(zip(hits,errors)):  # more generous errors for stiff helices
                            hist.SetPoint(i,point[0]*1.0e-3,point[1]*1.0e-3,point[2]*1.0e-3) # in [m]
                            hist.SetPointError(i, 3*err[0]*1.0e-3, 3*err[1]*1.0e-3, 3*err[2]*1.0e-3)
                        fitresult = None
                        step = 0.1
                        counter = 0
                        while fitresult is None and counter < 20:
                            startvalue = -1.04 + counter * step
                            hf = pfit.helixfitter(hist,self.bfield, startvalue) # with inverse radius start value
                            fitresult = hf.helix_result # sensitive to inverse radius start value!
                            counter += 1
                            if fitresult is not None and abs(fitresult.par[3]) > 3.0: # failed fit convergence
                                fitresult = None
                        if (fitresult is not None): # store HelixFit object
                            fr.add_fitter(pnumber,fitresult) # added HelixFit
                            #print 'Found with start value %f'%startvalue
#                            print fitresult
                        else:
                            hist.Clear() # even more generous errors
                            for i,(point,err) in enumerate(zip(hits,errors)):  # more generous 10% error for stiff helices
                                hist.SetPoint(i,point[0]*1.0e-3,point[1]*1.0e-3,point[2]*1.0e-3) # in [m]
                                hist.SetPointError(i, 0.1*point[0]*1.0e-3, 0.1*point[1]*1.0e-3, 3*err[2]*1.0e-3)
                            counter = 0
                            while fitresult is None and counter < 20:
                                startvalue = -1.06 + counter * step
                                hf = pfit.helixfitter(hist,self.bfield, startvalue) # with inverse radius start value
                                fitresult = hf.helix_result # sensitive to inverse radius start value!
                                counter += 1
                                if fitresult is not None and abs(fitresult.par[3]) > 3.0: # failed fit convergence
                                    fitresult = None
                            if (fitresult is not None): # store HelixFit object
                                fr.add_fitter(pnumber,fitresult) # added HelixFit
                                #print 'Found 2nd time with start value %f'%startvalue
                            #else:
                                #print "Failed MIGRAD helix fit"
                hist.Clear()
            candidates.append(fr)
        return candidates



class MS_FittingService2D(object):
    """
    Fitting Service object - fit track candidates
    Input: must find dictionary of track candidates
    
    Output: List of FitResults objects

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
        Output: List of FitResults objects
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
        candidates = []
        for cpaths in data:
            fr = FitResults(cpaths.id)
            # results = {}
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
                    # results[l] = []
                    fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                    # results[l].append((bestfit,error,chi,bins,angles))
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
                        fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                        # results[l].append((bestfit,error,chi,bins,angles))

                        ic2, sl2 = self.start_values(hits2)
                        #print 'Right part: Start values: ic=%f, sl=%f'%(ic2,sl2)
                        bestfit, error, chi = msfit.segment_fit(hits2,ic2,sl2)
                        fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                        # results[l].append((bestfit,error,chi,bins,angles))
                        
                else: # initial fit failed
                    fr.add_fitter(l, (bestfit,error,chi,[],[]))
                    # results[l] = [(bestfit,error,chi,[],[])]
            candidates.append(fr)
            # cluster[cpaths.id] = results
        return candidates


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
    
    Output: List of FitResults objects

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
        Output: List of FitResults objects
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
        candidates = []
        for cpaths in data:
            fr = FitResults(cpaths.id)
            # results = {}
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
                    # results[l] = []
                    fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                    # results[l].append((bestfit,error,chi,bins,angles))
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
                            if c1>0 and c2>0:
                                # concatenate
                                bestshort = b1+b2
                                errshort = e1+e2
                                tchi = c1+c2
                                fr.add_fitter(l, (bestshort,errshort,tchi,bins,angles))
                            #fr.add_fitter(l, (b2,e2,c2,bins,angles))
                            # results[l].append((b1,e1,c1,bins,angles))
                            # results[l].append((b2,e2,c2,bins,angles))
                        else:
                            hits1 = hits[:bin1]
                            #print 'Hits to left: ',hits1
                            icxy, slxy = self.start_values(hits1,0)
                            icxz, slxz = self.start_values(hits1,1)
                            #print 'Start values: icxy=%f, slxy=%f'%(icxy,slxy)
                            #print 'Start values: icxz=%f, slxz=%f'%(icxz,slxz)
                            stvals = [icxy,slxy,icxz,slxz]
                            bestfit, error, chi, d1,d2 = ms3d.fitter(hits1,stvals)
                            if chi>0:
                                fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                            # results[l].append((bestfit,error,chi,bins,angles))

                        if bin2>len(hits)-6:
                            hits2 = hits[bin2-5:] # have at least 5 points
                            #print 'short data2: ',hits2
                            b1, e1, c1, b2, e2, c2 = self.short_fit(hits2)
                            if c1>0 and c2>0:
                                # concatenate
                                bestshort = b1+b2
                                errshort = e1+e2
                                tchi = c1+c2
                                fr.add_fitter(l, (bestshort,errshort,tchi,bins,angles))
                            #fr.add_fitter(l, (b1,e1,c1,bins,angles))
                            #fr.add_fitter(l, (b2,e2,c2,bins,angles))
                            # results[l].append((b1,e1,c1,bins,angles))
                            # results[l].append((b2,e2,c2,bins,angles))
                        else:
                            hits2 = hits[bin2:]
                            #print 'Hits to right: ',hits2
                            icxy, slxy = self.start_values(hits2,0)
                            icxz, slxz = self.start_values(hits2,1)
                            #print 'Start values: icxy=%f, slxy=%f'%(icxy,slxy)
                            #print 'Start values: icxz=%f, slxz=%f'%(icxz,slxz)
                            stvals = [icxy,slxy,icxz,slxz]
                            bestfit, error, chi, d1,d2 = ms3d.fitter(hits2,stvals)
                            if chi>0:
                                fr.add_fitter(l, (bestfit,error,chi,bins,angles))
                            # results[l].append((bestfit,error,chi,bins,angles))

                    #else: # initial fit failed
                    #fr.add_fitter(l, (bestfit,error,chi,[],[]))
                    # results[l] = [(bestfit,error,chi,[],[])]
                    candidates.append(fr)
            # cluster[cpaths.id] = fr
        return candidates


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




