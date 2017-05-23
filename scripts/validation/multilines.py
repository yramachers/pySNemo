from math import sqrt, tan, pi, atan
import numpy as np
import random
from pysnemo.utility import euclid
from pysnemo.utility import helix as HX


class demonstratorgrid(object):
    """
    This object should construct and hold the wire grid coordinates
    no input
    Data attribute: grid matrix of coordinate tuples
    Output: several hit tuples and truth values from hits() function
    """
    def __init__(self):
        self.rows = 113  # tracker geometry
        self.columns = 9  # tracker geometry
        self.d = 44.0  # [mm] sense wire distance
        self.offsetx = 53.0 # [mm] for wire plane 0 closest to foil
        self.offsety = -2464.0 # [mm] bottom wire row

        # have a 2d array of given dimensions
        self.grid_left = np.zeros((self.rows,self.columns),dtype=('f4,f4')) # tracker side = 0
        self.grid_right = np.zeros((self.rows,self.columns),dtype=('f4,f4')) # tracker side = 1
        self.wireinfo = []

        # with origin bottom left
        m = 0
        for entry in self.grid_left: # rows
            for n in range(len(entry)): # columns
                tup = (n*self.d + self.offsetx, m*self.d + self.offsety)
                self.grid_right[m][n] = tup

                tup2 = (-(n*self.d + self.offsetx), m*self.d + self.offsety)
                self.grid_left[m][n] = tup2
            m += 1


    def __str__(self):
        s = "grid with %d columns and %d rows\n" % (self.columns,self.rows)
        s += "in a lattice with lattice constant %f cm" % self.d
        return s
        
    
    def hits(self, line=None, side=0):
        """
        Input: A euclid line3 object
               The tracker side (=0 -> left tracker default)
        Process: Determine the nearest wires given a line throught the 
                 initial grid and cut on those wires closer than the 
                 lattice constant.
        Output: 2 Lists of wire coordinates and radii = shortest distance 
                wire to line
        """
        zvec = euclid.Vector3(0.0,0.0,1.0) # zaxis vector for 3D wire line
        darr = np.zeros((self.rows, self.columns))

        if side < 1:
            tracker = self.grid_left
        else:
            tracker = self.grid_right

        self.wireinfo = [] # clean start
        m = 0
        for entry in tracker: # rows
            for n in range(len(entry)): # columns
                if isinstance(line, euclid.Line3): # input was a Line3 object
                    pos = euclid.Point3(entry[n][0], entry[n][1], 0.0) # point of wire in grid
                    wire = euclid.Line3(pos,zvec) # wire into a vertical 3D line
                    darr[m][n] = wire.connect(line).length # shortest distance line to line in 3D

                elif isinstance(line, HX.helix): # input was a Helix object
                    distance = line.GetDistanceToPoint((entry[n][0]*1.0e-3, entry[n][1]*1.0e-3, 0.0)) # input in [m]
                    darr[m][n] = distance[0]*1.0e3 # [mm] shortest distance helix to wire point in 2D

                else:
                    print 'grid ERROR: input type unknown, not line nor helix'
                    return [], []
            m += 1
                
        cells = []
        radius = []
        iarr,jarr = np.where(darr <= (0.51*self.d)) #allow minimal overlap
        for ind,item in enumerate(iarr):
            jind = jarr[ind]
            self.wireinfo.append((side, jind, item)) # side, row, column
            cells.append(tracker[item][jind])
            radius.append(darr[item][jind])
  
        return cells, radius


    def multi_track_hits(self, lines = []):
        """
        Input: A list of tuples: (euclid line3 object, side)
               
        Process: Find all wires for each line on a specific tracker side. 
                 identical hits in separate clusters are fine.
                 
        Output: cluster dictionary, counter from one, each track a cluster
                containing (cells, radii, wire_info) as tuple of lists
        """
        cells = []
        radius = []
        cluster = { } # each line is a cluster

        for counter, (l, s) in enumerate(lines):
            cells, radius = self.hits(l,s)
            cluster[counter+1] = (cells, radius, self.wireinfo)
            
        if len(cluster)>1:
            return self.remove_doubles(cluster)
        else:
            return cluster


    def remove_doubles(self, cluster):
        for k in range(1,len(cluster)): # keys start at 1
            for pos, mi in enumerate(cluster[k][2]): # check info entries against all of next cluster entry
                nextinfo = np.array(cluster[k+1][2])  # from modified next cluster
                if len(nextinfo)<1:
                    break # nothing to compare to
                dublet = list(mi) # comes as tuple -> list convert
                nilist = nextinfo[:,-3:].tolist() # smart slicing, back to list

                if dublet in nilist: # found a double entry
                    indx = nilist.index(dublet) # first entry; should be only one
                    nextradius = cluster[k+1][1][indx] # choose according to radius
                    if cluster[k][1][pos] > nextradius: # smaller wins
                        cluster[k][1][pos] = nextradius # in-place change
                    cluster[k+1][0].pop(indx) # change cluster k+1
                    cluster[k+1][1].pop(indx)
                    cluster[k+1][2].pop(indx)

        return cluster




class track_generator(object):
    """
    simple toy track generator. Consider only 2D objects in z=0 plane
    no input
    Output: Euclid 3D line objects depending on chosen generator method
    """
    def __init__(self):
        self.lines = [] # line objects container
        

    def __str__(self):
        s = 'contains %d entries in lines container\n'%len(self.lines)
        for l in self.lines:
            s += str(l) # assumes line knows how to print itself
            s += '\n'
        return s


    def getLines(self):
        return self.lines


    def single_line_manual(self, slope, intercept = 0.0):
        # no container needed
        vec = euclid.Vector3(1.0, slope,0.0)
        pos = euclid.Point3(0.0, intercept, 0.0) # fix at x=0 plane -> foil plane
        return euclid.Line3(pos,vec)


    def single_line_random_slope(self, intercept = 0.0):
        # no container needed
        angle = random.uniform(-pi*0.5+0.17, pi*0.5-0.17) # taking vertical out
        sl = tan(angle)
        return self.single_line_manual(sl, intercept)


    def single_line_manual_atplane(self, slope, interceptx, intercepty):
        # no container needed
        vec = euclid.Vector3(1.0, slope,0.0)
        pos = euclid.Point3(interceptx, intercepty, 0.0) # fix at icx, iy point
        return euclid.Line3(pos,vec)


    def single_line_random_atplane(self, interceptx, intercepty):
        # no container needed
        angle = random.uniform(-pi*0.5+0.17, pi*0.5-0.17) # taking vertical out
        sl = tan(angle)
        return self.single_line_manual_atplane(sl, interceptx, intercepty)


    def double_manual_atvertex(self, sl1, sl2, intercept = 0.0):
        self.lines = [] # clear
        # v-shape double line from vertex
        self.lines.append(self.single_line_manual(sl1, intercept))
        self.lines.append(self.single_line_manual(sl2, intercept))
        

    def double_random_atvertex(self, intercept = 0.0):
        self.lines = [] # clear
        # v-shape double line from vertex
        self.lines.append(self.single_line_random_slope(intercept))
        self.lines.append(self.single_line_random_slope(intercept))
        

    def haystack(self, nlines=2):
        self.lines = [] # clear
        if nlines<1:
            nlines = 1 # minimum 1 line
        
        for i in range(nlines):
            intercept = random.uniform(-2464.0,2464.0) # limit from demonstrator y-axis
            self.lines.append(self.single_line_random_slope(intercept))



class helix_generator(object):
    """
    simple helix generator
    Input: None
    Output: Helix object from utility helix.py with generated properties, 
            to be used to get wire hits from grid object.
    """
    def __init__(self):
        self.helices = [] # helix container


    def single_random_momentum(self, intercept = 0.0, side=0):
        pos = (0.0, intercept * 1.0e-3, 0.0) # unit [m] for helix object

        px = random.uniform(0.3,1.4) # random momentum x, right
        py = random.uniform(-0.1,0.1) # random momentum y 
        pz = 0.0 # try only 2D on wire distance

        if side==1: # right tracker
            charge = -1.0 # unit [e]
        else:
            charge = 1.0 # unit [e]

        momentum = (px,py,pz) # unit [MeV/c]
        bfield = 2.5e-3 # 25 Gauss

        return HX.helix(pos,momentum,charge,bfield)
