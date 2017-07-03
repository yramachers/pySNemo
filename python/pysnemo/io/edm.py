from pysnemo.utility import euclid
class line(object):
    """
    simple line object holding slope and intercept
    """
    def __init__(self,slope=1.0,icept=0.0):
        self.slope = slope
        self.intercept = icept
        v2 = euclid.Vector2(1.0,self.slope)
        p1 = euclid.Point2(0.0,self.intercept)
        self.l = euclid.Line2(p1,v2)

    def __str__(self):
        return "line with slope=%f, intercept=%f" % (self.slope,self.intercept)

    def get_euclid_line(self):
        return self.l



class ClusterPaths(object):
    '''
    Hold cluster and associated paths: 
    - cluster id
    - the tracker side of the hits as sign of path id
    - potential paths.
    '''
    def __init__(self,clusterid):
        self.id = clusterid
        self.paths = {} # paths with id as key and points as list of tuples
        
        
    def add_path(self, id, pathlist):
        self.paths[id] = pathlist

    


class FitResults(object):
    '''
    Holds fit results objects
    - cluster id
    - the tracker side of the hits and 
    - list of fit objects like tuples for lines and HelixFit objects.
    '''
    def __init__(self,clusterid):
        self.id = clusterid
        self.fitterpaths = {} # paths with id as key and list of fitter objects
        
        
    def set_side(self, sideid):
        if sideid<0:
            self.side = 0 # tracker side negative x
        else:
            self.side = 1 # tracker side positive x


    def add_fitter(self, id, fitter):
        if not self.fitterpaths.has_key(id):
            self.fitterpaths[id] = [] # fresh list at key id
        self.set_side(id)
        self.fitterpaths[id].append(fitter)

    


class Particle(object):
    '''
    Hold final particle information: 
    - extrapolated vertices to foil and calo (with errors, all in one tuple)
    - the associated calo hit object
    - the fit object, i.e. either tuple(s) for a line
      or a HelixFit.
    - a tuple as id (cluster id, path id, type)
      where type is 0: line fit; 1: Helix fit; 2: broken line fit
    '''
    def __init__(self,particleid):
        self.id = particleid
        self.kink = False # for broken line fits
        self.angles = []
        
        
    def __str__(self):
        s = "Particle object, id = (%d, %d, %d)" % self.id
        s += "\nFoil vertex tuple: "+str(self.foilvertex)
        s += "\nCalorimeter vertex tuple: "+str(self.calovertex)
        s += "\nAssociated Calo: "+"\n"+str(self.calo_hit)
        s += "\nChi2: %f"%self.chi2
        return s

    def set_vertex_foil(self, foiltuple):
        self.foilvertex = foiltuple


    def set_vertex_calo(self, calotuple):
        self.calovertex = calotuple


    def set_calo_hit(self, calohit):
        self.calo_hit = calohit

    
    def set_fitter(self, fitobj):
        self.fitter = fitobj
        if isinstance(fitobj,tuple):
            self.chi2 = fitobj[2] # line fit
        else:
            self.chi2 = fitobj.chi2 # helix fit


    def add_fitter_for_BL(self, fitobj):
        dummylist = [self.fitter]
        dummylist.append(fitobj)
        self.fitter = dummylist # change type to list
        self.chi2 = 0.0
        for entry in self.fitter:
            if isinstance(entry,tuple): # get a total chi2 value
                self.chi2 += entry[2] # line fit
            else:
                self.chi2 += entry.chi2 # helix fit


    def set_kink(self, flag):
        self.kink = flag


    def set_angles(self,angles):
        self.angles = angles


class gcylinder(object):
    """
    The Geiger cylinder only contains data originating from measurements. That 
    data needs to be handed over when constructing the cylinder.
      xwire, ywire : geiger counter wire coordinates, known precisely
      radius       : measured radial distance of track to wire
      dr           : radial uncertainty
      zwire        : coordinate along the wire, measured, usually to far less 
                     precision than radial value
      dz           : error on z
      r_max        : maximum cylinder radius to avoid overlap; constant

    It also holds meta data as tuple of int in the container called 
    meta_info as (Id, Side, Layer and Column)
    """
    def __init__(self,x=0.0,y=0.0,z=0.0,dz=1.0,r=0.0,dr=1.0):
        self.xwire = x
        self.ywire = y
        self.zwire = z
        self.dz = dz
        self.checkRadius(r,dr)
        
    def __str__(self):
        s = "Geiger Cylinder of radius = %f and height = %f\n" % (self.radius,2*self.dz)
        s += "Uncertainty on radius dr = %f\n" % self.dr
        s += "around wire position (%f, %f, %f)" % (self.xwire,self.ywire,self.zwire)
        return s

    def __eq__(self,other):
        '''
        check on equality of two cylinders, self and other, tolerance set
        to 1 mm of the wire coordinates, assumed to be in units of cm
        '''
        tol = 0.1
        if (abs(self.xwire - other.xwire) < tol):
            if (abs(self.ywire - other.ywire) < tol):
                if (abs(self.zwire - other.zwire) < tol):
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False


    def checkRadius(self,r,dr):
        """
        The radius can be pointlike but has a finite error, i.e. check if 
        its too small for its error and set it to the 
        error boundary, i.e. avoid zero radii values.
        """
        if (r < dr):
            self.radius = dr # minimal radius
        else:
            self.radius = r

        self.dr = dr
    

    def set_info(self, ggid, id, module, side, layer, column):
        self.meta_info = (ggid, id, module, side, layer, column) # all int type

                 

class tracker_hit(object):
    """
    The tracker Geiger cylinder contains data after calibration.
    x,y,z,sigmaz : x,y,z wire coordinates and error on z
    r,sigmar     : r, sigmar ring radius and error

    It also holds meta data as tuple of int in the container called 
    meta_info as (Id, module, Side, Layer and Column)
    """
    def __init__(self,x,y,z,sz,r,sr):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.sigmar = sr
        self.sigmaz = sz
        
    def __str__(self):
        s = "Geiger Cylinder of radius = %f and height = %f\n" % (self.r,2*self.sigmaz)
        s += "Uncertainty on radius dr = %f\n" % self.sigmar
        s += "around wire position (%f, %f, %f)\n" % (self.x,self.y,self.z)
        s += "Meta Info: (%d,%d,%d,%d,%d,%d)"%(self.meta_info[0],self.meta_info[1],self.meta_info[2],self.meta_info[3],self.meta_info[4],self.meta_info[5])
        return s
    

    def set_info(self, ggid, trueid, module, side, layer, column):
        self.meta_info = (ggid, trueid, module, side, layer, column)

                 

class tangentpoint_hit(object):
    """
    The tracker Geiger cylinder tangent point has position and
    errors and meta information relating to its origin like 
    a cylinder meta_info block and cylinder ids involved 
    in its creation.
    x,y,z            : x,y,z point coordinates and error on z
    ex,ey,ez         : errors in coordinates

    It also holds meta data as tuple of tuple, tuple in the 
    container called meta_info as 
    ((creation cylinder ids),(on cylinder meta_info))
    """
    def __init__(self,x,y,z,ex,ey,ez):
        self.x = x
        self.y = y
        self.z = z
        self.ex = ex
        self.ey = ey
        self.ez = ez
        
    def __str__(self):
        s = "Tangent point coords: (%f,%f,%f)\n" % (self.x,self.y,self.z)
        s += "Uncertainties (%f,%f,%f)\n" % (self.ex,self.ey,self.ez)
        s += "Meta Info: "
        s += str(self.meta_info[1])
        s += str(self.meta_info[2])
        return s
    

    def set_info(self, trig, ids, mi):
        self.meta_info = (trig, ids, mi)

                 

class gcylinder_truth(object):
    """
    The truth Geiger cylinder only contains data originating 
    from truth monte carlo data.
    time   : hit time
    x,y,zstart : x,y,z step start in cell
    x,y,zstop  : x,y,z step stop in cell

    It also holds meta data as tuple of int in the container called 
    meta_info as (Id, module, Side, Layer and Column, trackid, parentid)
    """
    def __init__(self,time,xstart,ystart,zstart,xstop,ystop,zstop):
        self.time = time
        self.x0 = xstart
        self.y0 = ystart
        self.z0 = zstart
        self.x1 = xstop
        self.y1 = ystop
        self.z1 = zstop
        
    def __str__(self):
        s = "Geiger Truth Cylinder"
        s += "Start (x,y,z = (%f,%f,%f)\n" % (self.x0,self.y0,self.z0)
        s += "Stop (x,y,z = (%f,%f,%f)\n" % (self.x1,self.y1,self.z1)
        s += "Time = %f\n" % self.time
        return s
    

    def set_info(self, ggid, module, side, layer, column, trid, parentid):
        self.meta_info = (ggid, module, side, layer, column, trid, parentid)

                 

class calo_hit(object):
    """
    A calo_hit object holds simulated calorimeter data in
    time    : time of hit on calorimeter
    sigma_t : error on time
    energy  : deposited total energy from hit on calorimeter
    sigma_e : energy error

    It also holds meta data as tuple of int in the container called 
    meta_info as (Id, Type, Module, Side, Column, Row and Wall)
    """
    def __init__(self,time, dt, energy, de):
        self.time = time
        self.sigma_t = dt
        self.energy = energy
        self.sigma_e = de
        
    def __str__(self):
        s = "Calorimeter Hit at time = %f +- %f\n" % (self.time,self.sigma_t)
        s += "and energy = %f +- %f\n" % (self.energy,self.sigma_e)
        s += "meta info: id=%d, type=%d, module=%d, side=%d, column=%d, row=%d, wall=%d" % self.meta_info
        return s

    def set_info(self, id, type, module, side, column, row, wall):
        self.meta_info = (id, type, module, side, column, row, wall) # all int type


class calo_truth_hit(object):
    """
    A calo_truth_hit object holds simulated calorimeter truth data in
    time    : time of hit on calorimeter
    x, y ,z : hit coordinates on calo
    energy  : deposited total energy from hit on calorimeter

    It also holds meta data as tuple of int in the container called 
    meta_info as (Id, type, module, Side, Column, Row and Wall)
    """
    def __init__(self,time, x, y, z, energy):
        self.time = time
        self.x = x
        self.y = y
        self.z = z
        self.energy = energy
        
    def __str__(self):
        s = "Calorimeter Hit at time = %f\n" % self.time
        s += "and location (x,y,z) = (%f,%f,%f)\n" % (self.x,self.y,self.z)
        s += "and energy = %f\n" % self.energy
        s += "meta info: id=%d, type=%d, module=%d, side=%d, column=%d, row=%d, wall=%d" % self.meta_info
        return s

    def set_info(self, id, type, module, side, column, row, wall):
        self.meta_info = (id, type, module, side, column, row, wall) # all int


class truevertex(object):
    """True vertex storage from file. Holds
    time : production time
    xc: x coordinate of vertex
    yc: y -"-
    zc: z -"-
    """
    def __init__(self, x, y, z, time):
        self.time = time
        self.xc = x
        self.yc = y
        self.zc = z
        
    def __str__(self):
        s = "Vertex at = (%f,%f,%f)\n" % (self.xc,self.yc,self.zc)
        s += "and time = %f\n" % self.time
        return s



class trueparticle(object):
    """ True particle momentum storage from file. Holds
    time : production time
    px: x momentum of particle
    py: y -"-
    pz: z -"-
    """
    def __init__(self, id, type, px, py, pz, time, energy):
        self.id = id
        self.type = type
        self.px = px
        self.py = py
        self.pz = pz
        self.time = time
        self.kinenergy = energy
        
    def __str__(self):
        s = "True particle momentum at = (%f,%f,%f)\n" % (self.px,self.py,self.pz)
        s += "and time = %f\n" % self.time
        s += "and id, type = %d, %d\n" % (self.id,self.type)
        s += "and energy = %f\n" % self.kinenergy
        return s



class toytruth(object):
    """ Toy simulation truth storage from file. Holds 6 parameter
    and an id counter. The 6 numbers make up originally a euclid line3 object
    but can be used for other geometric structures too, like a helix.
    """
    def __init__(self, id, par, bplist):
        self.id = id
        self.bplist = bplist # contains break points (bpx,bpy,bplayer) if any
        vec = euclid.Vector3(par[0], par[1], par[2])
        pt = euclid.Point3(par[3], par[4], par[5])
        self.line = euclid.Line3(pt,vec) # convenience, not useful for helix truth
        self.parstore = par # in case the bare data is requested

    def __str__(self):
        s = "Toy simulation truth data\n"
        s += 'id = %d '%self.id
        s += 'parameter: (%f, %f, %f, %f, %f, %f, %d)\n'%self.parstore
        for i in range(len(self.bplist)): # only if bplist not empty
            s += 'break point x=%f, y=%f with angle=%f at layer %d.'%(self.bplist[i][0],self.bplist[i][1],self.bplist[i][2],self.bplist[i][3])
        return s

    def getLine(self):
        return self.line # as euclid Line3 object

    def getParameters(self):
        return self.parstore # as tuple

    def getCharge(self):
        return self.parstore[-1] # final entry

    def getbreakpointlist(self):
        return self.bplist
