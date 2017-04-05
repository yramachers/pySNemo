from math import sqrt, tan, pi, atan
import numpy as np
import random
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


class line3D(object):
    '''
    simple 3D line object holding slope and intercept in 
    (x-y: p0 and p1) and 
    (x-z: p2 and p3) planes
    Output: triplet of cartesian coordinates of a point on the line
            defined by 4-tuple of parameter values
    '''
    def __init__(self,p=[1.0,1.0,1.0,1.0]):
        '''
        A parameteric line is defined from 6 parameters but 4 are independent.
        x0,y0,z0,x1,y1,z1 are the coordinates of two points on the line
        can choose x0 = 0 if line not parallel to y-z plane (=foil plane)
        and x1 = 1; 
        Input: p=4-tuplet of line parameter values
        Process: provides function to return a 
                 point_on_line(t=scalar line parameter)
                 as a cartesian triplet (x,y,z).
        '''
        self.par = p
        self.line = line(self.par[1],self.par[0]) # xy plane projection line
        
    def __str__(self):
        s = "line with xy slope=%f, intercept=%f\n" % (self.par[1],self.par[0])
        s +="line with xz slope=%f, intercept=%f" % (self.par[3],self.par[2])
        return s

    def point_on_line(self, t):
        x = t
        y = self.par[0] + self.par[1]*t
        z = self.par[2] + self.par[3]*t
        return x,y,z


class grid3D(object):
    """
    This object should construct and hold the wire grid coordinates
    Input: columns and rows of the grid and wire distance [cm]
    Data attribute: grid matrix of coordinate tuples
    Output: tuple of x,y wire coordinates 
    """
    def __init__(self,columns=9,rows=9,lattice=4.4):
        self._rows = rows
        self._columns = columns
        self.d = lattice
        # have a 2d array of given dimensions
        self.grid = np.zeros((rows,columns),dtype=('f4,f4,f4'))
        self.wireinfo = []
        # with origin bottom left
        m = 0
        for entry in self.grid:
            n = 0
            for pair in entry:
                pair['f0'] = n*lattice
                pair['f1'] = m*lattice
                self.grid[m][n] = (pair['f0'],pair['f1'],0) #triplet in grid
                n += 1
            m += 1

    def __str__(self):
        s = "grid with %d columns and %d rows\n" % (self._columns,self._rows)
        s += "in a lattice with lattice constant %f cm" % self.d
        return s
        
    def distance_to_line(self,line=None,point=(0.0,0.0,0.0)):
        """
        Input: a line3D object and a tuple of x,y,z coordinates
        Output: the distance of point to line in x,y-plane, 
                for now ignores the z-coordinate
        """
        l = line.line
        distance = abs(-l.intercept + point[1] - l.slope*point[0])/sqrt(l.slope*l.slope + 1.0)

        el = l.get_euclid_line()
        p1 = euclid.Point2(float(point[0]),float(point[1]))
        
        ls = euclid._connect_point2_line2(p1,el)
        #truthpoint = (ls.p.x+ls.v.x, ls.p.y+ls.v.y)
        xval = ls.p.x+ls.v.x
        ponline = line.point_on_line(xval)
        truthpoint = (xval, ls.p.y+ls.v.y, ponline[2])
        return distance, truthpoint


    def get_info(self):
        return self.wireinfo


    def hits(self,line=None):
        """
        Input: A line object
        Process: Determine the nearest wires given a line throught the 
                 initial grid and cut on those wires closer than the 
                 lattice constant.
        Output: 2 Lists of wire coordinates and radii = shortest distance 
                wire to line
        """
        #print line
        # have a 2d array of given dimensions
        trutharray = np.zeros((self._rows,self._columns),dtype=('f4,f4'))
        ll = np.zeros((self._rows,self._columns),dtype='f4')
        zval = {}
        m = 0
        for entry in self.grid:
            n = 0
            for pair in entry:
                distance,truth = self.distance_to_line(line,(pair['f0'],pair['f1'],0.0))
                ll[m][n] = distance
                trutharray[m][n] = (truth[0],truth[1])
                zval[(m,n)] = truth[2]
                n += 1
            m += 1
                
        cells = []
        truthpoints = []
        zlist = []
        radius = []
        self.wireinfo = [] # clean start
        iarr,jarr = np.where(ll <= (0.51*self.d)) #allow minimal overlap
        for ind,item in enumerate(iarr):
            jind = jarr[ind]
            zlist.append(zval[(item,jind)])
            cells.append(self.grid[item][jind])
            radius.append(ll[item][jind])
            truthpoints.append(trutharray[item][jind])
            self.wireinfo.append((0, jind, item)) # side column, layer
        # concatenate
        newcells = []
        for entry,z in zip(cells,zlist):
            x = entry[0]
            y = entry[1]
            #print 'cells: ',(x,y,z)
            newcells.append((x,y,z))
        return newcells,radius,truthpoints


class generator3D(object):
    """
    simple 3D line generator
    Input: Wire_grid; generates a single line with random slope and intercept
           in the uper right quadrant consistent with wire_grid object
    Output: Returns results from the hits method in wire_grid, i.e. 
            list of wire coordinates closest to line and list of radii
            as closest approach to wire.
    """
    def __init__(self, grid=None):
        # Initialisation only
#        self.centre = (0.5*grid._rows*grid.d+0.25*grid.d,0.5*grid._columns*grid.d+0.25*grid.d)
        self.centre = (0.0,0.5*grid._rows*grid.d+0.25*grid.d) # left border
        self.slope = []
        self.intercept = []
        self.angle = 0
        self._grid = grid
        random.seed() # system time as seed

    def __str__(self):
        s = ""
        for sl,ic in zip(self.slope,self.intercept):
            s += "line with slope=%f, intercept=%f\n" % (sl,ic)
        return s


    def generate(self):
        self.slope = []
        self.intercept = []
        par = []
        angle = pi/2.0
        while (angle>(0.5*pi-0.08) or angle<0.2):# cut vertical out
            angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
        sl = tan(angle)
        ic = self.centre[1] - sl*self.centre[0]
        self.slope.append(sl)
        self.intercept.append(ic)
        par.append(ic)
        par.append(sl)
        angle = pi/2.0
        while (angle>(0.5*pi-0.08) or angle<0.2):# cut vertical out
            angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
        sl = tan(angle) # xz slope
        par.append(ic)
        par.append(sl)
        self.slope.append(sl)
        self.intercept.append(ic)
        
        l3 = line3D(par)
        c, r, t = self._grid.hits(l3)
        return c, r, t


    def haystack(self, nlines = 2):
        self.slope = []
        self.intercept = []

        if nlines<1:
            nlines = 1 # no less than 1 line accepted.

        c, r, t = self.generate()  #first centre line
            
        # additional lines
        for i in range(1,nlines):
            par = []
            frac = random.uniform(0.2,0.8) # next starting point on foil
            self.centre = (0.0,frac*self._grid._rows*self._grid.d+0.25*self._grid.d) # new
            angle = pi/2.0
            while (angle>(0.5*pi-0.08) or angle<0.2):# cut vertical out
                angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
            sl = tan(angle) # xy slope
            ic = self.centre[1] - sl*self.centre[0] # common intercept
            self.slope.append(sl)
            self.intercept.append(ic)
            par.append(ic)
            par.append(sl)
            angle = pi/2.0
            while (angle>(0.5*pi-0.08) or angle<0.2):# cut vertical out
                angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
            sl = tan(angle) # xz slope
            par.append(ic)
            par.append(sl)
            self.slope.append(sl)
            self.intercept.append(ic)

            l3 = line3D(par)
            c2, r2, t2 = self._grid.hits(l3)

            c += c2 # appends list to previous list
            r += r2
            t += t2
        return c,r,t



class line_generator(object):
    """
    simple in-plane line generator
    Input: Wire_grid; generates a single line with random slope and intercept
           in the uper right quadrant consistent with wire_grid object
    Output: Returns results from the hits method in wire_grid, i.e. 
            list of wire coordinates closest to line and list of radii
            as closest approach to wire.
    """
    def __init__(self, grid=None):
        # Initialisation only
#        self.centre = (0.5*grid._rows*grid.d+0.25*grid.d,0.5*grid._columns*grid.d+0.25*grid.d)
        self.centre = (0.0,0.5*grid._rows*grid.d+0.25*grid.d) # left border
        self.slope = []
        self.intercept = []
        self.angle = 0
        self._grid = grid
        self.meta_info = None
        random.seed() # system time as seed

    def __str__(self):
        s = ""
        for sl,ic in zip(self.slope,self.intercept):
            s += "line with slope=%f, intercept=%f\n" % (sl,ic)
        return s


    def generate(self):
        self.slope = []
        self.intercept = []
        angle = pi/2.0
        while (angle>(0.5*pi-0.08) or angle<(-0.5*pi+0.08)):# cut vertical out
            angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
        sl = tan(angle)
        ic = self.centre[1] - sl*self.centre[0]
        self.angle = angle
        #print "angle: ",self.angle
        par = [ic,sl,ic,0.0] # no slope in xz plane
        self.slope.append(sl)
        self.intercept.append(ic)
        self.slope.append(0.0)
        self.intercept.append(ic)
        l = line3D(par)
        c, r, t = self._grid.hits(l)
        self.meta_info = self._grid.get_info()
        return c, r, t


    def vertex_double(self):
        self.slope = []
        self.intercept = []

        #first random line
        c, r, t = self.generate()

        # second line mirrored
        sl = tan(-self.angle)
        ic = self.intercept[0] # copy
        self.slope.append(sl)
        self.intercept.append(ic)
        self.slope.append(0.0)
        self.intercept.append(ic)
        par = [ic,sl,ic,0.0] # no slope in xz plane
        l2 = line3D(par)
        c2,r2, t2 = self._grid.hits(l2)
        meta_info2 = self._grid.get_info()
        self.meta_info += meta_info2 # concatenate
        
        # done
        c.reverse()
        r.reverse()
        t.reverse()
        c2.reverse()
        r2.reverse()
        t2.reverse()
        self.angle *= 2.0 # make it an opening angle
        return c+c2, r+r2, t+t2

    def haystack(self, nlines = 2):
        self.slope = []
        self.intercept = []

        if nlines<1:
            nlines = 1 # no less than 1 line accepted.

        c, r, t = self.generate()  #first centre line
            
        # additional lines
        for i in range(1,nlines):
            frac = random.uniform(0.2,0.8) # next starting point on foil
            self.centre = (0.0,frac*self._grid._rows*self._grid.d+0.25*self._grid.d) # new

            angle = pi/2.0
            while (angle>(0.5*pi-0.08) or angle<0.2):# cut vertical out
                angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
            sl = tan(angle)
            ic = self.centre[1] - sl*self.centre[0]
            self.slope.append(sl)
            self.intercept.append(ic)
            self.slope.append(0.0)
            self.intercept.append(ic)
            par = [ic,sl,ic,0.0] # no slope in xz plane
            l = line3D(sl,ic)
            c2,r2, t2 = self._grid.hits(l)
            c += c2 # appends list to previous list
            r += r2
            t += t2
        return c,r,t


    def manual_vertex_double(self, angle, vertexlocation):
        self.slope = []
        self.intercept = []
        
        sl = tan(angle)
        self.slope.append(sl)
        self.intercept.append(vertexlocation)
        self.slope.append(0.0)
        self.intercept.append(vertexlocation)
        par = [vertexlocation,sl,vertexlocation,0.0] # no slope in xz plane

        l = line3D(par)
        c, r, t = self._grid.hits(l)

        # second line
        sl = tan(-angle)
        self.slope.append(sl)
        self.intercept.append(vertexlocation)
        self.slope.append(0.0)
        self.intercept.append(vertexlocation)
        par = [vertexlocation,sl,vertexlocation,0.0] # no slope in xz plane

        l2 = line3D(par)
        c2, r2, t2 = self._grid.hits(l2)

        # done
        c.reverse()
        r.reverse()
        t.reverse()
        c2.reverse()
        r2.reverse()
        t2.reverse()
        self.angle = 2.0*angle # make it an opening angle
        return c+c2, r+r2, t+t2


    def manual_clock(self, angle, vertexlocation):
        self.slope = []
        self.intercept = []
        
        sl = tan(angle)
        self.slope.append(sl)
        self.intercept.append(vertexlocation)
        self.slope.append(0.0)
        self.intercept.append(vertexlocation)
        par = [vertexlocation,sl,vertexlocation,0.0] # no slope in xz plane

        l = line3D(par)
        c, r, t = self._grid.hits(l)

        # second line
        sl = 0.0  # fixed horizontal line
        self.slope.append(sl)
        self.intercept.append(vertexlocation)
        self.slope.append(0.0)
        self.intercept.append(vertexlocation)
        par = [vertexlocation,sl,vertexlocation,0.0] # no slope in xz plane

        l2 = line3D(par)
        c2, r2, t2 = self._grid.hits(l2)

        # done
        c.reverse()
        r.reverse()
        t.reverse()
        c2.reverse()
        r2.reverse()
        t2.reverse()
        self.angle = angle # make it an opening angle
        return c+c2, r+r2, t+t2


class ms_line_generator(object):
    """
    multiple scatter line generator.
    Produces its own wire_grid, 
    generates a single line with random multiple scatter events, 
    nscatter, as parameter and scatter angle normal std as second
    parameter.
    Output: Returns results from the hits method in wire_grid, i.e. 
            list of wire coordinates closest to line and list of radii
            as closest approach to wire.
    """
    def __init__(self, columns=9, rows=9, lattice=4.4):
        self.rows = rows
        self.columns = columns
        self.d = lattice
        self.slope = [] 
        self.intercept = []
       # have a 2d array of given dimensions
        self.grid = np.zeros((rows,columns),dtype=('f4,f4,f4'))
        # with origin bottom left
        m = 0
        for entry in self.grid:
            n = 0
            for pair in entry:
                pair['f0'] = n*lattice
                pair['f1'] = m*lattice
                self.grid[m][n] = (pair['f0'],pair['f1'],0) #triplet in grid
                n += 1
            m += 1

        self.centre = (0.0,0.5*self.rows*self.d+0.25*self.d) # left border
        random.seed() # system time as seed



    def distance_to_line(self, line, point):
        """
        Input: a euclid line object and a euclid point
        Output: the distance of point to line in x,y-plane, 
        """
        ls = euclid._connect_point2_line2(point,line)
        distance = ls.length
        truthpoint = (ls.p.x+ls.v.x, ls.p.y+ls.v.y)
        return distance, truthpoint


    def hits(self,line=None):
        """
        Input: A euclid line object
        Process: Determine the nearest wires given a line throught the 
                 initial grid and cut on those wires closer than the 
                 lattice constant.
        Output: 2 Lists of wire coordinates and radii = shortest distance 
                wire to line
        """
        trutharray = np.zeros((self.rows,self.columns),dtype=('f4,f4'))
        ll = np.zeros((self.rows,self.columns),dtype='f4')
        m = 0
        for entry in self.grid:
            n = 0
            for pair in entry:
                p = euclid.Point2(pair['f0'],pair['f1'])
                distance,truth = self.distance_to_line(line,p)
                ll[m][n] = distance
                trutharray[m][n] = (truth[0],truth[1]) #dublet in array
                n += 1
            m += 1
                
        cells = []
        truthpoints = []
        radius = []
        iarr,jarr = np.where(ll <= (0.51*self.d)) #allow minimal overlap
        for ind,item in enumerate(iarr):
            jind = jarr[ind]
            cells.append(self.grid[item][jind])
            radius.append(ll[item][jind])
            truthpoints.append(trutharray[item][jind])
        return cells,radius,truthpoints



    def trajectory_length(self, leu):
        vup = euclid.Vector2(0.0,1.0)
        pleft = euclid.Point2(0.0,0.0)
        pright = euclid.Point2(self.columns*self.d,0.0)
        grid_left = euclid.Line2(pleft,vup)
        grid_right = euclid.Line2(pright,vup)

        pl = leu.intersect(grid_left)
        pr = leu.intersect(grid_right)
        lseg = pl.connect(pr)
        return abs(lseg), lseg



    def multiple_scatter(self, nscatter=1, scatter_angle=3.0):
        self.slope = [] 
        self.intercept = []
        self.angle = pi/2.0
        while (self.angle>(0.5*pi-0.08) or self.angle<(0.2)):# cut vertical out
            self.angle = 0.5*random.vonmisesvariate(0.0,0) #uniform (0,pi/2)
        sl = tan(self.angle)
        ic = self.centre[1] - sl*self.centre[0]
        l = line(sl,ic)
        # book original slope, intercept
        self.slope.append(sl)
        self.intercept.append(ic)
        
        original = l.get_euclid_line()
            
        # break points
        lines = []
        breakpoints = []
        startpoint = original.p1
        frac = oldfrac = 0.0
        for i in range(nscatter+1): # number of lines needed
            length, lsegment = self.trajectory_length(original)
            while (oldfrac>=frac):
                frac = random.uniform(0.1,0.8) # breakpoint fraction on line
            oldfrac = frac # increasing fraction break points

            x = lsegment.p.x+lsegment.v.x*frac
            y = lsegment.p.y+lsegment.v.y*frac
            pbreak = euclid.Point2(x,y)
            ll = euclid.Line2(original.p1,pbreak) # shortened line
            lines.append(ll)
            breakpoints.append(pbreak)

            scat_angle = random.gauss(0.0,scatter_angle*pi/180)
            self.angle += scat_angle
            sl = tan(self.angle)
            v2 = euclid.Vector2(1.0,sl)
            original = euclid.Line2(pbreak,v2) # new line, scattered

        breakpoints.pop() # remove last
        # remove for more than one breakpoint!!
        self.breakpoint = breakpoints
        centres = []
        radii = []
        truthlist = []
        for entry in lines:
            c,r,t = self.hits(entry)
            centres.append(c)
            radii.append(r)
            truthlist.append(t)
        c,r,t = self.disentangle(centres,radii,truthlist,breakpoints)
        return c,r,t


    def disentangle(self, cs, rs, ts, bp):
        centres = []
        radii = []
        truth = []
        previous = euclid.Point2(-1.0,0.0)

        for i,entry in enumerate(bp):
            print 'break point: ',entry
            cl = cs[i]
            rl = rs[i]
            tl = ts[i]
            for c,r,t in zip(cl,rl,tl):
                # collect according to column value interval
                boundx = (entry.x+0.5*self.d,previous.x-0.5*self.d)
                if (c[0]<boundx[0] and c[0]>boundx[1]):
                    #print 'y-test'
                    #print 'wire value: ',c[1]
                    #print 'test on <entry.y+d: ',entry.y+1.5*self.d
                    #print 'test on >previous.y-d: ',previous.y-1.5*self.d
                    boundy = (entry.y+1.5*self.d,previous.y-1.5*self.d)
                    if (c[1]<=boundy[0] and c[1]>=boundy[1]):
                        if c in centres: # double wire
                            rind = centres.index(c) # where in list
                            print 'found double wire at:',rind
                            if r<radii[rind]:
                                radii[rind] = r # pick smaller radius
                                truth[rind] = t
                                # leave rest unchanged
                        else:
                            centres.append(c)
                            radii.append(r)
                            truth.append(t)
                            # print 'wires passed cut: ',c

            previous = entry
        cl = cs[i+1]
        rl = rs[i+1]
        tl = ts[i+1]
        for c,r,t in zip(cl,rl,tl):
            # collect according to column value interval
            #print 'final y-test'
            #print 'wire value: ',c[1]
            #print 'test on previous.y-d: ',previous.y-1.5*self.d
            bound = (previous.x-0.5*self.d,previous.y-1.5*self.d)
            if (c[0]>bound[0] and c[1]>=bound[1]):
                if c in centres: # double wire
                    rind = centres.index(c) # where in list
                    print 'final, found double wire at:',rind
                    if r<radii[rind]:
                        radii[rind] = r # pick smaller radius
                        truth[rind] = t
                        # leave rest unchanged
                else:
                    centres.append(c)
                    radii.append(r)
                    truth.append(t)
                    # print 'wires passed cut: ',c
                    
        return centres, radii, truth


    def manual_feed(self, slope, intercept, nscatter=1, scatter_angle=3.0):
        self.slope = [] 
        self.intercept = []
        sl = slope
        ic = intercept
        self.angle = atan(sl)
        l = line(sl,ic)
        # book original slope, intercept
        self.slope.append(sl)
        self.intercept.append(ic)
        
        original = l.get_euclid_line()
            
        # break points
        lines = []
        breakpoints = []
        startpoint = original.p1
        frac = oldfrac = 0.0
        for i in range(nscatter+1): # number of lines needed
            length, lsegment = self.trajectory_length(original)
            while (oldfrac>=frac):
                frac = random.uniform(0.1,0.8) # breakpoint fraction on line
            oldfrac = frac # increasing fraction break points

            x = lsegment.p.x+lsegment.v.x*frac
            y = lsegment.p.y+lsegment.v.y*frac
            pbreak = euclid.Point2(x,y)
            ll = euclid.Line2(original.p1,pbreak) # shortened line
            lines.append(ll)
            breakpoints.append(pbreak)

            scat_angle = random.gauss(0.0,scatter_angle*pi/180)
            self.angle += scat_angle
            sl = tan(self.angle)
            v2 = euclid.Vector2(1.0,sl)
            original = euclid.Line2(pbreak,v2) # new line, scattered

        breakpoints.pop() # remove last
        centres = []
        radii = []
        truthlist = []
        for entry in lines:
            c,r,t = self.hits(entry)
            centres.append(c)
            radii.append(r)
            truthlist.append(t)
        c,r,t = self.disentangle(centres,radii,truthlist,breakpoints)
        return c,r,t



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

