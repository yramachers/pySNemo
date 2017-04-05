from math import sqrt, tan, pi, atan
import numpy as np
import random
from pysnemo.utility import euclid

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
        self.offsety = -2464.0 # [mm] bottom of foil

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

        m = 0
        for entry in tracker: # rows
            for n in range(len(entry)): # columns
                pos = euclid.Point3(entry[n][0], entry[n][1], 0.0) # point of wire in grid
                wire = euclid.Line3(pos,zvec) # wire into a vertical 3D line
                darr[m][n] = wire.connect(line).length # shortest distance line to line in 3D
            m += 1
                
        self.wireinfo = [] # clean start
        cells = []
        radius = []
        iarr,jarr = np.where(darr <= (0.51*self.d)) #allow minimal overlap
        for ind,item in enumerate(iarr):
            jind = jarr[ind]
            self.wireinfo.append((side, jind, item)) # side column, layer
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
        angle = pi/2.0
        while (angle>(0.5*pi-0.08) or angle<(-0.5*pi+0.08)):# cut vertical out
            angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
        sl = tan(angle)
        return self.single_line_manual(sl, intercept)


    def double_manual_atvertex(self, sl1, sl2, intercept = 0.0):
        self.lines = [] # clear
        self.lines.append(self.single_line_manual(sl1, intercept))
        self.lines.append(self.single_line_manual(sl2, intercept))
        

    def double_random_atvertex(self, intercept = 0.0):
        self.lines = [] # clear
        self.lines.append(self.single_line_random_slope(intercept))
        self.lines.append(self.single_line_random_slope(intercept))
        

    def haystack(self, nlines=2):
        self.lines = [] # clear
        if nlines<1:
            nlines = 1 # minimum 1 line
        
        for i in range(nlines):
            intercept = random.uniform(-2464.0,2464.0) # limit from demonstrator y-axis
            self.lines.append(self.single_line_random_slope(intercept))

