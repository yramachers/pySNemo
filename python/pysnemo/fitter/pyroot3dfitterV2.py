import ROOT as root
import pysnemo.utility.helix as hhelix
root.gSystem.Load("/home/epp/phsdaq/Code/pySNemo/python/pysnemo/fitter/fit3dVers2_cxx.so")

class linefitter(object):
    '''
    Input: Data in form of a ROOT object, TGraph2DErrors, holding
           3D data, x, y, z and associated errors ex, ey, ez.
           Only used as a convenient container here
    Process: Fit a straight line to the data in 3D.
    Output: best fit = linepar, a list
            associated errors from fit = fit_errors, a list
            *** Can implement more here if need be ***
    '''
    def __init__(self,datagraph=None):
        self.linepar = [] # solutions
        self.fit_errors = [] # corresponding errors
        self.lfitter = None
        self.line_fit(datagraph)


    def __str__(self):
        if len(self.linepar)<1:
            return "Apply line_fit method first: results are empty!"
        else:
            s = "Fit results:\n"
            s += "Intercept in x-y plane = %f +- %f\n" % (self.linepar[0],self.fit_errors[0])
            s += "    Slope in x-y plane = %f +- %f\n" % (self.linepar[1],self.fit_errors[1])
            s += "Intercept in x-z plane = %f +- %f\n" % (self.linepar[2],self.fit_errors[2])
            s += "    Slope in x-z plane = %f +- %f\n" % (self.linepar[3],self.fit_errors[3])
            s += "Goodness of fit chi2 = %f\n" %self.chi2
            return s

    def line_fit(self,data):
        '''
        Fill the attributes linepar and fit_errors as results of a fit.
        '''
        lf = root.line_fit(data)
        if (lf):
            parfit = []
            parerr = []
            for i in range(4):
                parfit.append(lf.GetParameter(i))
                parerr.append(lf.GetParError(i))

            self.linepar = parfit
            self.fit_errors = parerr
            chi = root.Double()
            d1 = root.Double()
            d2 = root.Double()
            i1 = root.Long()
            i2 = root.Long()
            lf.GetStats(chi,d1,d2,i1,i2)
            self.chi2 = chi
            self.lfitter = lf

        else:
            print "Warning: Fit failed with MIGRAD"
            self.linepar = []
            self.fit_errors = []
            self.chi2 = -1.0



class helixfitter(object):
    '''
    Input: Data in form of a ROOT object, TGraph2DErrors, holding
           3D data, x, y, z and associated errors ex, ey, ez.
           Also Bfield is offered as second argument
    Process: Fit a helix to the data in 3D.
    Output: Access by helix_result member which is a HelixFit object.
    '''
    def __init__(self, datagraph, bfield = 25.0):
        self.helix_result = None
        self.bfield = bfield # magnetic field strength in [Gauss]
        self.helix_fit(datagraph)
        

    def __str__(self):
        if self.helix_result is not None:
            return "Apply helix_fit method first: results are empty!"
        else:
            s = "Print HelixFit object in helix_result for output.\n"
            return s

    def helix_fit(self,data):
        '''
        Fill the attributes linepar and fit_errors as results of a fit.
        '''
        lf = root.helix_fit(data, self.bfield)
        if (lf):
            test = hhelix.HelixFit(lf)
            if test.hel is not None:
                self.helix_result = test
        else:
            print "Warning: Fit failed with MIGRAD"
