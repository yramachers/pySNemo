import ROOT
from array import array
from scipy.optimize import curve_fit
import numpy as np
from math import sqrt

class Leastsquare(ROOT.TPyMultiGenFunction):
    ''' 
    2D line fit function including multiple scattering term for 
    broken line algorithm, V.Blobel, NIM A566 (2006) 14
    Expand by allowing to choose model function as member function in class.
    So far only 2D line implemented - easy to expand this class 
    for choice.
    '''
    def __init__(self):
        ROOT.TPyMultiGenFunction.__init__(self,self)

    def NDim(self):
        return 2 # number of fit parameters

    def DoEval(self, par):
        try:
            x=self.datax[0]
        except IndexError:
            print 'Must set the data first'
            
        sum = 0.0

        # least square term first with errors on x and y!
        for x,y,ey,ex in zip(self.datax,self.datay,self.erry,self.errx):
            u = self.linefunction(x,par)
            sum += (y-u)**2/(ey**2 + par[1]**2 * ex**2)
 
        # additional term
        # build the error first from errx, erry - make self.errbeta
        self.calculate_angle_error(par)
        # build the angles - make self.beta
        self.calculate_beta(par)
        for i in range(1,len(self.datax)-1):
            sum += self.beta[i-1]*self.beta[i-1]/(self.errbeta[i-1])

        return sum

    def linefunction(self,x,par):
        # model function,
        # par list containing N model parameter
        # x the independent function parameter
        # here assume a 2D line with par[1] = slope, par[0] = intercept
        return par[1] * x + par[0]


    def calculate_beta(self,par):
        self.beta = [] # clear list
        for i in range(1,len(self.datax)-1): # N-2 slopes
            # residuals
            uim1 = self.datay[i-1]-self.linefunction(self.datax[i-1],par)
            ui = self.datay[i]-self.linefunction(self.datax[i],par)
            uip1 = self.datay[i+1]-self.linefunction(self.datax[i+1],par)

            thetai = (uip1-ui)/(self.datax[i+1]-self.datax[i])
            thetaim1 = (ui-uim1)/(self.datax[i]-self.datax[i-1])
            self.beta.append(thetai - thetaim1)
         

    def calculate_angle_error(self,par):
        self.errbeta = [] # clear list
        for i in range(1,len(self.erry)-1): # N-2 errors
            ui = self.datay[i]-self.linefunction(self.datax[i],par)
            uip1 = self.datay[i+1]-self.linefunction(self.datax[i+1],par)
            uim1 = self.datay[i-1]-self.linefunction(self.datax[i-1],par)

            value = self.datax[i+1] - self.datax[i]
            dvsq = self.erry[i+1]**2 + self.erry[i]**2
            dtheta2 = dvsq / value**2
            value2 = uip1 - ui
            dvsq = self.errx[i+1]**2 + self.errx[i]**2
            dtheta2 += value2**2 * dvsq / value**4

            value = self.datax[i] - self.datax[i-1]
            dvsq = self.erry[i]**2 + self.erry[i-1]**2
            dtheta1 = dvsq / value**2
            value2 = ui - uim1
            dvsq = self.errx[i]**2 + self.errx[i-1]**2
            dtheta1 += value2**2 * dvsq / value**4


            #value = uip1 - ui
            #dvsq = self.errx[i+1]**2 + self.errx[i]**2
            #dtheta2  = dvsq * value**2
            #value = self.datax[i+1] - self.datax[i]
            #dvsq = self.erry[i+1]**2 + self.erry[i]**2
            #dtheta2 += dvsq * value**2

            #value = ui - uim1
            #dvsq = self.errx[i]**2 + self.errx[i-1]**2
            #dtheta1  = dvsq * value**2
            #value = self.datax[i] - self.datax[i-1]
            #dvsq = self.erry[i]**2 + self.erry[i-1]**2
            #dtheta1 += dvsq * value**2

            self.errbeta.append(dtheta1+dtheta2)
        

    def SetData(self, datax, datay):
        self.datax = datax
        self.datay = datay


    def SetErrors(self, errx, erry):
        self.errx = errx
        self.erry = erry

        
        
def fitter(data,ic,sl):
    '''
    Input: data for one track and 
           start values: intercept and slope
           Expand by allowing Errors - needs more detailed error calc
           in Leastsquare object to get to errors on beta from x,y.
    '''
    # data preparation into x and y
    datax = []
    datay = []
    errx = []
    erry = []
    for d in data:
        datax.append(d[0])
        datay.append(d[1])
        errx.append(d[3])
        erry.append(d[4])

    # try using root fitting
    min = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    min.SetMaxFunctionCalls(1000000)
    min.SetTolerance(0.001)
    min.SetPrintLevel(0) # verbose level

    lstsq = Leastsquare()
    lstsq.SetData(datax,datay)
    lstsq.SetErrors(errx,erry) 

    step = [0.01,0.01]
    variable = [ic,sl] # start values
    
    min.SetFunction(lstsq)
    
    # Set the free variables to be minimized!
    min.SetVariable(0,"intercept",variable[0], step[0])
    min.SetVariable(1,"slope",variable[1], step[1])
    if min.Minimize():
        # get some results from minimizer
        best = array('f',[0]) # pyroot pointer helper
        best = min.X()
        err = array('f',[0]) # pyroot pointer helper
        err = min.Errors()
        val = min.MinValue() / (len(data)-lstsq.NDim()-1)

        # angle,error list after fit finished, i.e. final lists
        beta_angles = lstsq.beta
        beta_errors = lstsq.errbeta

        blist = []
        errlist = []
        for i in range(lstsq.NDim()):
            #print 'Fit parameter %d: %f +- %f'%(i,best[i],err[i])
            blist.append(best[i])
            errlist.append(err[i])
        return blist, errlist, val, beta_angles, beta_errors
    else:
        print 'Minuit minimizer failed!'
        return [], [], -1, [], []



def kink_finder(betaangles, berrors):
    '''
    Attempts to locate kinks (usually just one) using the error list 
    of angles along fitted multiple scattering model. Peaks in those 
    errors along fitted model point at kink locations.
    '''
    n = len(betaangles)
    hist = ROOT.TH1D("hist","title",n,1,n)
    for i,(val,err) in enumerate(zip(betaangles,berrors)):
        hist.SetBinContent(i+1,abs(val))
        hist.SetBinError(i+1,sqrt(err))
    sp = ROOT.TSpectrum()
    npeaks = sp.Search(hist,1,"goff",0.8) # 80% threshold for second peak
    print "found npeaks = ",npeaks
    bins = []
    angles = []
    if npeaks>0:
        peakpos = array('f',[0]) # pyroot pointer helper
        peakpos = sp.GetPositionX()
        for i in range(npeaks):
            bin = hist.FindBin(peakpos[i])
            print 'peak at bin: ',bin
            print 'and angle: ',hist.GetBinContent(bin)
            bins.append(bin)
            angles.append(hist.GetBinContent(bin))
    return bins, angles



def segment_fit(data,ic,sl):
    '''
    Input: data for one track and 
           start values: intercept and slope
    '''
    # data preparation into x and y
    datax = []
    datay = []
    errx = []
    erry = []
    for d in data:
        datax.append(d[0])
        datay.append(d[1])
        errx.append(d[3])
        erry.append(d[4])
    x = np.array(datax)
    y = np.array(datay)
    ex = np.array(errx)
    ey = np.array(erry) 
    # only one error as weight permitted - sum squared
    err = np.sqrt(np.square(ex)+np.square(ey)) 
    # start values
    p0 = np.array([ic,sl])
    # fit
    popt, pcov = curve_fit(model_func, x, y, p0, err)
    if not isinstance(pcov,float):
        # errors
        errlist = [sqrt(pcov[0][0]), sqrt(pcov[1][1])] # diagonal
        # chisquare
        sum=0.0
        for i,(xv,yv) in enumerate(zip(x,y)):
            nom = model_func(xv,popt[0],popt[1]) - yv
            denom = ey[i]**2 + popt[1]**2 * ex[i]**2
            sum += nom*nom/denom
        sum /= float(len(x)-2) # for a line take two dof out
    else:
        print 'failed fit'
        errlist = [0.0,0.0]
        sum = -1.0
    return list(popt), errlist, sum


def model_func(x, ic, sl):
    return sl * x + ic
