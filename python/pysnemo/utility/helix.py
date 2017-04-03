import ROOT as root
from math import sqrt, sin, cos, atan, atan2, pi, isnan
class helix(object):
    '''
    Helix object, from
    http://www.lcsim.org/software/pandora/tag-1.22/doxygen/html/classpandora_1_1Helix.html
    Changes for correct Helix calculations by YR, 02/08/2011
    Preferred helix object for event generators.
    '''
    def __init__(self,position=(0,0,0),momentum=(1,0,0),charge=1.0,bfield=0.5):
        self.referencePoint = position # length unit metre
        self.momentum = momentum # [MeV/c]
        self.charge = charge # elementary charge units
        self.Bfield = bfield # unit: Tesla
        self.unit_constant = 2.99792458e2 # for momentum units in MeV/c
        self.buildHelix()

    def __str__(self):
        s = "Helix with reference point: (%f,%f,%f) [m]\n" % self.referencePoint
        s += "Momentum: (%f,%f,%f) [MeV/c]\n" % self.momentum
        s += "Charge = %f [e] in a %f Tesla field" % (self.charge,self.Bfield)
        return s

    def buildHelix(self):
        '''
        from Pandora Helix constructor
        setting all the required internal variables for
        Distance to Point function below.
        '''
        if self.Bfield==0.0:
            print "No B-field - no helix"
            return 0
        kappa = self.unit_constant
        half_pi = pi/2.0
        sign = self.Bfield/abs(self.Bfield)

        (px,py) = (self.momentum[0],self.momentum[1])
        self.pxy = sign*sqrt(px*px + py*py)
        self.radius = self.pxy/(kappa*self.Bfield) # [m]
        self.omega = self.charge / self.radius
        self.tanlambda = self.momentum[2] / self.pxy
        self.phiMomRefPoint = atan2(py,px)

        (x,y) = (self.referencePoint[0],self.referencePoint[1])
        self.xCentre = x + self.radius*cos(self.phiMomRefPoint - half_pi*self.charge)
        self.yCentre = y + self.radius*sin(self.phiMomRefPoint - half_pi*self.charge)
        if self.charge>0:
            self.d0 = self.charge*self.radius - sqrt(self.xCentre*self.xCentre + self.yCentre*self.yCentre)
        else:
            self.d0 = self.charge*self.radius + sqrt(self.xCentre*self.xCentre + self.yCentre*self.yCentre)

        self.phiRefPoint = atan2(y - self.yCentre,x - self.xCentre)
        self.phiAtPCA = atan2(-self.yCentre,-self.xCentre)
        self.phi0 = -half_pi * self.charge + self.phiAtPCA
        while self.phi0<0:
            self.phi0 += 2.0*pi
        while self.phi0>=2.0*pi:
            self.phi0 -= 2.0*pi

        self.xAtPCA = self.xCentre + self.radius * cos(self.phiAtPCA)
        self.yAtPCA = self.yCentre + self.radius * sin(self.phiAtPCA)
        self.pxAtPCA = self.pxy*cos(self.phi0)
        self.pyAtPCA = self.pxy*sin(self.phi0)

        deltaphi = self.phiRefPoint - self.phiAtPCA
        if abs((self.radius*self.tanlambda - deltaphi))<1.0e-6:
            xCircles = 0
        else:
            xCircles = (-self.referencePoint[2]*self.charge) / (2.0*pi*(self.radius*self.tanlambda - deltaphi))
        if xCircles>=0.0:
            n1 = int(xCircles)
            n2 = n1+1
        else:
            n1 = int(xCircles) - 1
            n2 = n1+1

        if (abs(n1-xCircles)<abs(n2-xCircles)):
            nCircles = n1
        else:
            nCircles = n2

        self.z0 = self.referencePoint[2] + self.radius*self.tanlambda*(deltaphi + 2.0*pi*nCircles)


    def PointOnHelix(self, t=1.0):
        '''
        Get a point (triplet) on the parametrised helix; 
        parameter t=0 when initial point is at phase phi0.
        '''
        sign = self.Bfield/abs(self.Bfield)
        phi0 = self.phiRefPoint
        if (self.charge*sign)>0:
            phi = phi0 + t*(2.0*pi)
            dphi = phi - phi0
        else:
            phi = pi + phi0 - t*(2.0*pi)
            dphi = phi - phi0 - pi

        x = self.xCentre + sign*self.charge*self.radius*(cos(phi))
        y = self.yCentre + sign*self.charge*self.radius*(sin(phi))
        z = self.referencePoint[2] + sign*self.charge*self.radius*self.tanlambda*dphi
        return (x,y,z)


    def GetDistanceToPoint(self,point=(0,0,0)):
        '''
        Input: Point tuple (x,y,z)
        Output: Distance tuple
                x component: distance in R-Phi plane 
                y-component: distance along Z axis 
                z-component: 3D distance magnitude
        '''
        sign = self.Bfield/abs(self.Bfield)
        phi0 = atan2(self.referencePoint[1] - self.yCentre, self.referencePoint[0] - self.xCentre)

        nCircles = 0
        if (self.charge*sign)>0:
            phi = atan2(point[1] - self.yCentre, point[0] - self.xCentre)
            offset = (pi + phi0)/(2.0*pi)
        else:
            phi = atan2(point[1] - self.yCentre, point[0] - self.xCentre) + pi
            offset = (phi0)/(2.0*pi) - 0.5

        xCircles = offset + (sign*self.charge*(point[2]+self.referencePoint[2])) / (2.0*pi*(self.tanlambda*self.radius))
        if (xCircles==1.0):
            nCircles = int(xCircles-1)
        else:
            nCircles = int(xCircles)

        if (self.charge*sign)>0:
            dPhi = 2.0*pi*float(nCircles) + phi - phi0
        else:
            dPhi = 2.0*pi*float(nCircles) + phi - phi0 - pi

        zOnHelix = self.referencePoint[2] + sign*self.charge*self.radius*self.tanlambda*dPhi
        distZ = abs(zOnHelix - point[2])
        distXY = sqrt((self.xCentre - point[0])**2.0 + (self.yCentre - point[1])**2.0)
        distXY = abs(distXY - self.radius)

        distance = (distXY, distZ, sqrt(distXY*distXY + distZ*distZ))
        return distance


    def intersectionXY(self,plane=(0.0,1.0,0.0,1.0)):
        '''
        Get the point of intersection with a plane parallel to z-axis. In 
        case of no intersection, returns 'None'
        Input: Plane parameter, 2 points in plane and normal vector components
               (x0,y0,ax,ay)
        '''
        sign = self.Bfield/abs(self.Bfield)
        if self.tanlambda==0.0:
            return None
        else:
            signz = self.tanlambda/abs(self.tanlambda)
        (x0,y0,ay,ax) = plane # swapped ax, ay for correct surface normal
        AA = sqrt(ax*ax+ay*ay)
        BB = (ax*(x0 - self.xCentre) + ay*(y0 - self.yCentre)) / AA
        CC = ((x0-self.xCentre)**2.0 + (y0 - self.yCentre)**2.0 - self.radius*self.radius) / AA
        det = BB*BB - CC
        if (det<0.0):
            return None

        tt1 = -BB + sqrt(det)
        tt2 = -BB - sqrt(det)
        xx1 = x0 + tt1*ax
        yy1 = y0 + tt1*ay
        xx2 = x0 + tt2*ax
        yy2 = y0 + tt2*ay
        
        phi1 = atan2(yy1-self.yCentre, xx1-self.xCentre)
        phi2 = atan2(yy2-self.yCentre, xx2-self.xCentre)
        phi0 = atan2(self.referencePoint[1]-self.yCentre, self.referencePoint[0] - self.xCentre)
        dphi1 = phi1 - phi0
        dphi2 = phi2 - phi0

        if (dphi1<0 and self.charge*sign>0):
            dphi1 = dphi1 + 2.0*pi
        elif (dphi1>0 and self.charge*sign<0):
            dphi1 = dphi1 - 2.0*pi
        if (dphi2<0 and self.charge*sign>0):
            dphi2 = dphi2 + 2.0*pi
        elif (dphi2>0 and self.charge*sign<0):
            dphi2 = dphi2 - 2.0*pi
        
        tt1 = -self.charge * dphi1 * self.radius / self.pxy
        tt2 = -self.charge * dphi2 * self.radius / self.pxy

        if (tt1<0 or tt2<0):
            tt1 = abs(tt1)
            tt2 = abs(tt2)
        
        if (tt1<tt2):
            time = tt1
            point = (xx1,yy1,self.referencePoint[2]+signz*time*abs(self.momentum[2]))
        else:
            time = tt2
            point = (xx2,yy2,self.referencePoint[2]+signz*time*abs(self.momentum[2]))

        return point


    def intersectionZ(self,zplane=0.0):
        '''
        Get the point of intersection with a plane perpendicular to z-axis. In 
        case of no intersection or no z-momentum, returns 'None'
        Input: Z-coordinate of plane
        '''
        if (self.momentum[2]==0.0):
            return None

        sign = self.Bfield/abs(self.Bfield)
        time = -(zplane - self.referencePoint[2]) / self.momentum[2]
        phi0 = atan2(self.referencePoint[1]-self.yCentre, self.referencePoint[0] - self.xCentre)
        phi = phi0 - sign*self.charge*self.pxy*time / self.radius
        point = (self.xCentre+self.radius*cos(phi), self.yCentre+self.radius*sin(phi), zplane)
        return point


class Par5Helix(object):
    '''
    5 parameter constructor for a helix - 
    example of parametrized helix function
    '''
    def __init__(self,par=(0.0,0.0,0.0,0.0,0.0), Bfield=1.0):
        '''
        My parametrization as opposed to LEP parameterization:
        swap phi0 and d0 for x,y of start point
        '''
        unit_constant = 2.99792458e2 # for momentum units in MeV/c
        # Parameter
        xatpca = par[0]
        yatpca = par[1]
        z0 = par[2]
        omega = par[3]
        tanlambda = par[4]

        if omega==0:
            print "No curvature omega, no helix"
            return 0

        charge = omega/abs(omega)
        radius = 1.0 / abs(omega)
        phi0 = atan(-xatpca/yatpca)
        ref_point = (xatpca,yatpca,z0)
        pxy = unit_constant * Bfield * radius
        momentum = (pxy*cos(phi0),pxy*sin(phi0),tanlambda*pxy)
        
        self.helix = helix(ref_point,momentum,charge,Bfield)

    def PointOnHelix(self, t=1.0):
        return self.helix.PointOnHelix(t)

    def GetDistanceToPoint(self,point=(0,0,0)):
        return self.helix.GetDistanceToPoint(point)

    def intersectionXY(self,plane=(0.0,1.0,0.0,1.0)):
        return self.helix.intersectionXY(plane)

    def intersectionZ(self,zplane=0.0):
        return self.helix.intersectionZ(zplane)


class HelixFit(object):
    '''
    Input: A ROOT fitter object (TVirtualFitter) containing the results 
           of a fit of data to a helix, i.e. 5 best fit parameter values
           as for the Par5Helix object and their errors!
    Attributes: Calculates all helix parameters and errors = 
                become accessible by name.
    '''
    def __init__(self,hfitter=None):
        p0 = hfitter.GetParameter(0)
        e0 = hfitter.GetParError(0)
        p1 = hfitter.GetParameter(1)
        e1 = hfitter.GetParError(1)
        p2 = hfitter.GetParameter(2)
        e2 = hfitter.GetParError(2)
        p3 = hfitter.GetParameter(3)
        e3 = hfitter.GetParError(3)
        p4 = hfitter.GetParameter(4)
        e4 = hfitter.GetParError(4)
        chi = root.Double()
        d1 = root.Double()
        d2 = root.Double()
        i1 = root.Long()
        i2 = root.Long()
        hfitter.GetStats(chi,d1,d2,i1,i2)
        self.chi2 = chi
        self.bf = hfitter.GetParameter(5)
        self.par = (p0,p1,p2,p3,p4)
        self.errors = (e0,e1,e2,e3,e4)
        valid = True
        for val,err in zip(self.par,self.errors):
            if isnan(val) or isnan(err):
                valid = False
        if valid:
            self.hel = Par5Helix(self.par,self.bf) # helix from best fit values
            self.observables()
        else:
            self.hel = None


    def __str__(self):
        s = "\nHelix Fit object creates all helix attributes from\n"
        s += "the 5 hit parameters:\n"
        s += "x-ref = (%f +- %f)\n" % (self.par[0],self.errors[0])
        s += "y-ref = (%f +- %f)\n" % (self.par[1],self.errors[1])
        s += "z-ref = (%f +- %f)\n" % (self.par[2],self.errors[2])
        s += "inverse radius and charge sign= (%f +- %f)\n" % (self.par[3],self.errors[3])
        s += "slope in z= (%f +- %f)\n" % (self.par[4],self.errors[4])
        s += "chisq = %f\n" % self.chi2
        s += "\n"
        s += "For a magnetic field of %f Tesla, the following observables\n" %self.bf
        s += "are also calculated and accessible directly by name:\n"
        s += "HelixFit_object.fitmomentum [MeV/c]\n"
        s += str(self.fitmomentum)
        s += "\nand its errors HelixFit_object.fitmomentum_errors\n"
        s += str(self.fitmomentum_errors)
        s += "\nParticle charge as HelixFit_object.fitcharge = %d [e]\n" % self.fitcharge
        s += "with error: %f [e]\n" % self.fitcharge_error
        return s

    def observables(self):
        '''
        Calculate momentum and its uncertainties as well as the fitted charge.
        Accessible by 'fitmomentum', 'fitmomentum_errors' and 'fitcharge'
        as the interesting observables after a helix fit to data.
        Intersections with planes to follow later.
        '''
        # simply get the values from the helix object constructed earlier
        self.fitmomentum = self.hel.helix.momentum
        self.fitcharge = self.hel.helix.charge
        # error calculations
        domega = self.errors[3]
        dtanl = self.errors[4]
        dx = self.errors[0]
        dy = self.errors[1]
        dz= self.errors[2]
        
        self.fitcharge_error = self.fitcharge * abs(domega/self.par[3])

        dpxy = self.hel.helix.pxy * (domega/self.par[3])
        arg = -self.par[0]/self.par[1]
        darg = arg * sqrt((dx/self.par[0])**2.0 + (dy/self.par[1])**2.0)
        df1 = -arg / (1.0+arg**2.0)**1.5
        df2 = 1.0 / (1.0+arg**2.0)**1.5
        d_px = sqrt((dpxy/self.hel.helix.pxy)**2.0 + (df1*df1 * darg*darg))
        d_py = sqrt((dpxy/self.hel.helix.pxy)**2.0 + (df2*df2 * darg*darg))
#        d_pz = sqrt((dpxy/self.hel.helix.pxy)**2.0 + (dtanl/self.par[4])**2.0)
        d_px *= self.fitmomentum[0]
        d_py *= self.fitmomentum[1]
#        d_pz *= self.fitmomentum[2]
        self.fitmomentum_errors = (abs(d_px), abs(d_py), abs(dtanl))
