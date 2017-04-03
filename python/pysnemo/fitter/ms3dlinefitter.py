import ROOT
from array import array
from scipy.optimize import curve_fit
import numpy as np
from math import sqrt, pi, sin, asin
from pysnemo.utility.euclid import Vector3

class Leastsquare(ROOT.TPyMultiGenFunction):
    ''' 
    3D line fit function including multiple scattering term for 
    broken line algorithm, V.Blobel, NIM A566 (2006) 14
    Expand by allowing to choose model function as member function in class.
    '''
    def __init__(self):
        ROOT.TPyMultiGenFunction.__init__(self,self)

    def NDim(self):
        return 4 # number of fit parameters

    def DoEval(self, par):
        try:
            x=self.datax[1]
        except IndexError:
            print 'Must set the data first'
            
        sum = 0.0

        # least square term first with errors on x and y!
        for x,y,z,ex,ey,ez in zip(self.datax,self.datay,self.dataz,self.errx,self.erry,self.errz):
            point=[x,y,z]
            error=[ex,ey,ez]
            dsqr, wsqr, dummy1, dummy2 = self.linefunction(point,error,par)
            sum += dsqr/wsqr

        #print 'first sum = ',sum
        # additional term
        # build the angles - make self.beta
        self.calculate_beta(par)
        # build the error - make self.errbeta
        self.calculate_angle_error(par)
        for i in range(1,len(self.datax)-1):
            if self.errbeta[i-1]>0.0:
                sum += self.beta[i-1]*self.beta[i-1]/(self.errbeta[i-1])

        #print 'second sum = ',sum
        return sum

    def linefunction(self,point,error,par):
        ''' model function,
        par list containing N model parameter
        x the independent function parameter
        here assume a 3D line 
        par[0]: xy intercept
        par[1]: xy slope
        par[2]: xz ic
        par[3]: xz sl
        '''
        # assume never parallel to z,y plane = foil plane
        xp = Vector3(point[0],point[1],point[2])

        x0 = Vector3(0.0, par[0], par[2])
        x1 = Vector3(1.0, par[0]+par[1], par[2] + par[3])
        u = (x1-x0).normalized()
        # residual
        dvec = ((xp-x0).cross(u))

        # errors from xp in cross product
        weight = [(error[1]*u.z)**2+(error[2]*u.y)**2,(error[0]*u.z)**2+(error[2]*u.x)**2,(error[0]*u.y)**2+(error[1]*u.x)**2]
        model = xp + dvec
    
        dsqr = dvec.magnitude_squared()
        wsqr = sum(weight)
        
        return dsqr, wsqr, model, weight


    def calculate_beta(self,par):
        self.beta = [] # clear list
        for i in range(1,len(self.datax)-1): # N-2 angles
            point1=[self.datax[i-1],self.datay[i-1],self.dataz[i-1]]
            point2=[self.datax[i],self.datay[i],self.dataz[i]]
            point3=[self.datax[i+1],self.datay[i+1],self.dataz[i+1]]
            error=[0.0,0.0,0.0] # not needed here
            # residuals
            dummy1, dummy2, um1, dw = self.linefunction(point1,error,par)
            dummy1, dummy2, ui, dw  = self.linefunction(point2,error,par)
            dummy1, dummy2, up1, dw = self.linefunction(point3,error,par)

            v1 = ui-um1
            v2 = up1-ui

            angle = v2.angle(v1)
            if angle>pi/2.0: # backward to forward hemisphere
                angle -= pi
            #if not angle>=0.0:
            #    print 'v1= ',v1
            #    print 'v2= ',v2
            #    print 'points: '
            #    print point1,point2,point3
            self.beta.append(angle)


    def calculate_angle_error(self,par):
        self.errbeta = [] # clear list
        for i in range(1,len(self.erry)-1): # N-2 angle errors
            point1=[self.datax[i-1],self.datay[i-1],self.dataz[i-1]]
            point2=[self.datax[i],self.datay[i],self.dataz[i]]
            point3=[self.datax[i+1],self.datay[i+1],self.dataz[i+1]]
            error1=[self.errx[i-1],self.erry[i-1],self.errz[i-1]]
            error2=[self.errx[i],self.erry[i],self.errz[i]]
            error3=[self.errx[i+1],self.erry[i+1],self.errz[i+1]]

            # residual and error square (as a list) on connection vector
            dsqr1, wsq1, um1, w1 = self.linefunction(point1,error1,par)
            dsqr2, wsq2, ui , w2 = self.linefunction(point2,error2,par)
            dsqr3, wsq3, up1, w3 = self.linefunction(point3,error3,par)

            # the error square on connection vectors v1 and v2
            dv1 = [w1[0]+w2[0],w1[1]+w2[1],w1[2]+w2[2]]
            dv2 = [w3[0]+w2[0],w3[1]+w2[1],w3[2]+w2[2]]

            # prepare the error on the angle (from acos(v1.v2/|v1||v2|)
            v1 = ui-um1
            v2 = up1-ui
            angle = v2.angle(v1)
            if angle>pi/2.0:
                angle -= pi
            const = abs(sin(angle))

            # denominators
            denom1 = v1.magnitude()*v2.magnitude()
            denom2 = (v1.magnitude())**3 * v2.magnitude()
            denom3 = v1.magnitude() * (v2.magnitude())**3

            # derivative terms, for v1 and also v2
            tv1x = v2.x/denom1 - v1.x*(v1.dot(v2))/denom2
            tv1y = v2.y/denom1 - v1.y*(v1.dot(v2))/denom2
            tv1z = v2.z/denom1 - v1.z*(v1.dot(v2))/denom2
            tv2x = v1.x/denom1 - v2.x*(v1.dot(v2))/denom3
            tv2y = v1.y/denom1 - v2.y*(v1.dot(v2))/denom3
            tv2z = v1.z/denom1 - v2.z*(v1.dot(v2))/denom3

            # component error square terms
            cerrx = const * (dv1[0]*tv1x**2 + dv2[0]*tv2x**2)
            cerry = const * (dv1[1]*tv1y**2 + dv2[1]*tv2y**2)
            cerrz = const * (dv1[2]*tv1z**2 + dv2[2]*tv2z**2)
            sinerr = sqrt(cerrx + cerry + cerrz)
            if sinerr>1.0:
                sinerr = 1.0
            angle_error = asin(sinerr)
            if angle_error<1.0e-5: # cut tiny errors
                angle_error = 1.0e-5
            self.errbeta.append(angle_error)


    def SetData(self, datax, datay, dataz):
        self.datax = datax
        self.datay = datay
        self.dataz = dataz
        #print 'Data: '
        #for x,y,z in zip(self.datax,self.datay,self.dataz):
        #    print x,y,z


    def SetErrors(self, errx, erry, errz):
        self.errx = errx
        self.erry = erry
        self.errz = errz

        
        
def fitter(data,par):
    '''
    Input: data for one track and 
           start values: intercepts and slopes in xy and xz projection
    '''
    # data preparation into x and y
    datax = []
    datay = []
    dataz = []
    errx = []
    erry = []
    errz = []
    for d in data:
        datax.append(d[0])
        datay.append(d[1])
        dataz.append(d[2])
        errx.append(d[3])
        erry.append(d[4])
        errz.append(d[5])

    # try using root fitting
    min = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    min.SetMaxFunctionCalls(1000000)
    min.SetTolerance(0.001)
    min.SetPrintLevel(0) # verbose level

    lstsq = Leastsquare()
    lstsq.SetData(datax,datay,dataz)
    lstsq.SetErrors(errx,erry,errz) 

    step = [0.01,0.01,0.01,0.01]
    icxy = par[0]
    slxy = par[1]
    icxz = par[2]
    slxz = par[3]
    variable = [icxy,slxy,icxz,slxz] # start values
    
    min.SetFunction(lstsq)
    
    # Set the free variables to be minimized!
    min.SetVariable(0,"intercept xy",variable[0], step[0])
    min.SetVariable(1,"slope xy",variable[1], step[1])
    min.SetVariable(2,"intercept xz",variable[2], step[2])
    min.SetVariable(3,"slope xz",variable[3], step[3])
    if min.Minimize():
        # get some results from minimizer
        best = array('f',[0]) # pyroot pointer helper
        best = min.X()
        err = array('f',[0]) # pyroot pointer helper
        err = min.Errors()
        val = min.MinValue() / (len(data)-lstsq.NDim())

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
    npeaks = sp.Search(hist,1,"goff",0.9) # 90% threshold for second peak
    print "found npeaks = ",npeaks
    if npeaks==0: # try one more
        npeaks = sp.Search(hist,2.0,"goff",0.9) # 90% threshold for second peak
        print "2nd attempt npeaks = ",npeaks
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



def segment_fit2D(datax,datay,errx,erry,ic,sl):
    '''
    Input: data for one track and 
           start values: intercept and slope
    '''
    # data preparation into x and y
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
