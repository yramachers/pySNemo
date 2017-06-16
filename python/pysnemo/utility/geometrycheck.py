import pysnemo.utility.euclid as EU
from pysnemo.utility.interval import Interval

class demonstratorcalo(object):
    '''
    Demonstrator calorimeter object
    Collect some useful geometry data and functionality for
    simple requests.
    All information rests in the dictionary called store
    with format key: (type, side, wall)
    type: 0 calo wall; 1 xwall; 2 gamma veto
    side 0 negative x coordinates; 1 positive
    wall 0 negative y or z coordinates, 1 positive
    values: list of [euclid plane object, list of tuples of
            [(row, column, interval object, interval object),...]]
    '''
    def __init__(self):
        self.point = []
        self.store = { }
        self.build_store()

    def __str__(self):
        s = "Geometry service store for demonstrator data."
        return s


    def build_store(self):
        '''
        type = info[0]
        side = info[1]
        wall = info[2]
        un-used integer in info as -1

        calo intervals storage as list of tuples with format: 
        (row, column, interval object, interval object)
        with interval object ordered according to coordinate axis order x, y, z
        '''
        main_wall_left = (0, 0, -1)
        main_wall_right = (0, 1, -1)
        x_wall_backl = (1, 0, 0)
        x_wall_backr = (1, 1, 0)
        x_wall_frontl = (1, 0, 1)
        x_wall_frontr = (1, 1, 1)
        gamma_veto_topl = (2, 0, 1)
        gamma_veto_topr = (2, 1, 1)
        gamma_veto_bottoml = (2, 0, 0)
        gamma_veto_bottomr = (2, 1, 0)
        
        # store the boundary planes in list at first position
        self.store[main_wall_left] = [self._calorimeter_plane(main_wall_left),self._build_calo_intervals(main_wall_left)]
        self.store[main_wall_right] = [self._calorimeter_plane(main_wall_right),self._build_calo_intervals(main_wall_right)]
        self.store[x_wall_backl] = [self._calorimeter_plane(x_wall_backl),self._build_calo_intervals(x_wall_backl)]
        self.store[x_wall_frontl] = [self._calorimeter_plane(x_wall_frontl),self._build_calo_intervals(x_wall_frontl)]
        self.store[gamma_veto_topl] = [self._calorimeter_plane(gamma_veto_topl),self._build_calo_intervals(gamma_veto_topl)]
        self.store[gamma_veto_bottoml] = [self._calorimeter_plane(gamma_veto_bottoml),self._build_calo_intervals(gamma_veto_bottoml)]
        self.store[x_wall_backr] = [self._calorimeter_plane(x_wall_backr),self._build_calo_intervals(x_wall_backr)]
        self.store[x_wall_frontr] = [self._calorimeter_plane(x_wall_frontr),self._build_calo_intervals(x_wall_frontr)]
        self.store[gamma_veto_topr] = [self._calorimeter_plane(gamma_veto_topr),self._build_calo_intervals(gamma_veto_topr)]
        self.store[gamma_veto_bottomr] = [self._calorimeter_plane(gamma_veto_bottomr),self._build_calo_intervals(gamma_veto_bottomr)]


    def _calorimeter_plane(self, info):
        type = info[0] # calo type
        side = info[1] # the tracker side
        wall = info[2] # the wall
        
        # Main Wall
        if type==0:
            # main calo sitting at x = +- 43.5cm
            if side == 0:
                calo = EU.Plane(EU.Point3(-435,0.0,0.0),EU.Point3(-435,1.0,0.0),EU.Point3(-435,0.0,1.0)) 
            else:
                calo = EU.Plane(EU.Point3(435,0.0,0.0),EU.Point3(435,1.0,0.0),EU.Point3(435,0.0,1.0)) 

        # X Wall
        elif type==1:
            # xwall calo sitting at y = +- 250.55cm
            if wall == 0:
                calo = EU.Plane(EU.Point3(0.0, -2505.5,0.0),EU.Point3(1.0, -2505.5,0.0),EU.Point3(0.0,-2505.5,1.0)) 
            else:
                calo = EU.Plane(EU.Point3(0.0, 2505.5,0.0),EU.Point3(1.0, 2505.5,0.0),EU.Point3(0.0,2505.5,1.0)) 

        # Gamma Veto
        elif type==2:
            # gveto calo sitting at z = +- 155.0cm
            if wall == 0:
                calo = EU.Plane(EU.Point3(1.0,0.0, -1550.0),EU.Point3(0.0,1.0, -1550.0),EU.Point3(0.0,0.0, -1550.0)) 
            else:
                calo = EU.Plane(EU.Point3(1.0,0.0, 1550.0),EU.Point3(0.0,1.0, 1550.0),EU.Point3(0.0,0.0, 1550.0)) 
                
        return calo # a euclid plane object in 3D


    def _build_calo_intervals(self, info):
        type = info[0] # calo type
        side = info[1] # the tracker side
        wall = info[2] # the tracker walls
        
        storagelist = []
        # Main Wall
        if type==0:
             for column in range(20):
                for row in range(13):
                    yinit = column * 259.0 - 2590.0 # 259 for modules + offset
                    zinit = row * 259.0 - 1683.5 # 259 for modules + offset
            
                    if side == 0:
                        dx = Interval(-436.0, -434.0)
                    else:
                        dx = Interval(434.0, 436.0) # bracket the plane
                    dy = Interval(yinit + 1.49, yinit + 1.5 + 256.01)
                    dz = Interval(zinit + 1.49, zinit + 1.5 + 256.01)
                    storagelist.append((type,row,column,dx,dy,dz))
        # X Wall
        elif type==1:
             for column in range(2):
                for row in range(16):
                    if side == 0:
                        xinit = -(column * 203.0 + 29.) # 202 for module + 1mm gap + offset
                        dx = Interval(xinit - 200.0 - 1.01, xinit - 0.99)
                    else:
                        xinit = column * 203.0 + 29. # 202 for module + 1mm gap + offset
                        dx = Interval(xinit + 0.99, xinit + 1.01 + 200.0)
                
                    zinit = row * 212.0 - 1696.0 # 212 for modules + offset
            
                    if wall ==0:
                        dy = Interval(-2506.0, -2505.0)
                    else:
                        dy = Interval(2505.0, 2506.0)
                    dz = Interval(zinit + 1.749, zinit + 1.75 + 208.51)
                    storagelist.append((type,row,column,dx,dy,dz))

        # Gamma Veto
        else:
            leftgv = Interval(0,8) # 8 not included!
            if side == 0:
                xinit = - 4.9 # 4.9mm offset to source
                dx = Interval(xinit - 1.01 - 290.0, xinit - 0.99)
            else:
                xinit = 4.9 # 4.9mm offset to source
                dx = Interval(xinit + 0.99, xinit + 1.01 + 290.0)
            
            for column in range(16):
                if column in leftgv: # in first block of 8
                    yinit = column * 311.5 - 2497.25 # 311.5 for modules + offset
                else:
                    yinit = column * 311.5 - 2497.25 + 10.5# 311.5 for modules + offset
            
                dy = Interval(yinit + 1.749, yinit + 1.75 + 308.01)
                if wall ==0:
                    dz = Interval(-1551.0, -1549.0)
                else:
                    dz = Interval(1549.0, 1551.0)
                storagelist.append((type,1,column,dx,dy,dz))
		
        return storagelist



    def calo_id(self, point):
        '''
        takes euclid point object as input to return info ID tuple
        of a calo block or empty tuple if not hit.
        '''
        for k in self.store.keys():
            for item in self.store[k][1]: # go through calo intervals list
                type = item[0]
                row = item[1]
                column = item[2]
                dA = item[3] # first interval object
                dB = item[4] # second
                dC = item[5] # third
                if point.x in dA and point.y in dB and point.z in dC:
                    return (-1, k[0], 0, k[1], column, row, k[2]) # full calo meta info return
        return ()



    def get_point(self, idx=0):
        if self.point[idx] is not None and len(self.point[idx])>0:
            return self.point[idx] # euclid Point3 object


    def has_overlap(self, intA, intB, info):
        '''
        takes 2 interval objects as input for main axes range of ellipse.
        intA and int B must be in order of coordinate axes order x, y, z
        for the logic test to pass.
        '''
        type = info[0] # of calorimeter
    
        for k in self.store.keys():
            if type == k[0]: # found type in key
                for item in self.store[k][1]: # go through calo intervals list
                    dA = item[2] # first interval object
                    dB = item[3] # second
                    if (intA.overlap(dA)) and (intB.overlap(dB)):
                        return True
                    else:
                        return False
                    
    def calohits(self, structure, side=0):
        '''
        Takes input line and tracker side from artificial multilines generator
        to create corresponding artificial calorimeter hits
        line is a euclid line3 object.
        '''
        self.point = []
        ci = []
        for t,s,w in self.store.keys(): # check each calo type
            thispoint = ()
            #print 'calohits side = %d'%side
            if isinstance(structure,EU.Line3):
                plane = self.store[(t,s,w)][0]
                thispoint = structure.intersect(plane)
                if isinstance(thispoint,EU.Point3):
                    caloinfo = self.calo_id(thispoint)
                    if len(caloinfo) and caloinfo[3]==side: # not an empty tuple
                        self.point.append(thispoint) # now the right units
                        ci.append(caloinfo) # found the calo hit
                        #print 'calo info found: ',caloinfo
            else: # helix object, internal length units [m]
                plane = self.store[(t,s,w)][0]
                if t==0: # main wall
                    if s==side:
                        if side==0:
                            planetup = (-435.0*1.0e-3,0.0,1.0,0.0)# input to helix method
                            thispoint = structure.intersectionXY(planetup)
                        else:
                            planetup = (435.0*1.0e-3,0.0,-1.0,0.0)# input to helix method
                            thispoint = structure.intersectionXY(planetup)
                elif t==1: # xwall
                    if s==side:
                        if w==0:
                            planetup = (0.0,-2505.5*1.0e-3,0.0,1.0)# input to helix method
                            thispoint = structure.intersectionXY(planetup)
                        else:
                            planetup = (0.0,2505.0*1.0e-3,0.0,-1.0)# input to helix method
                            thispoint = structure.intersectionXY(planetup)
                elif t==2: # gveto
                    if s==side:
                        if w==0:
                            zplane = -1550.0 # input to helix method
                            thispoint = structure.intersectionZ(zplane*1.0e-3)
                        else:
                            zplane =  1550.0# input to helix method
                            thispoint = structure.intersectionZ(zplane*1.0e-3)

                if thispoint is not None and len(thispoint)>0:
                    p = EU.Point3(thispoint[0],thispoint[1],thispoint[2])
                    p *= 1.0e3 # [m] to [mm]
                    #print 'test point: ',p
                    caloinfo = self.calo_id(p)
                    if len(caloinfo): # not an empty tuple
                        self.point.append(p) # now the right units
                        ci.append(caloinfo) # found the calo hit
                    #return caloinfo # found the calo hit
        return ci # not found


    def multi_calohits(self, lines):
        '''
        Takes input lines with tracker side from artificial multilines generator
        to create corresponding artificial calorimeter hits
        lines is a list of euclid line3 objects with side integers
        '''
        cluster = { } # each line is a cluster

        for counter, (l, s) in enumerate(lines):
            #print l
            #print 'on side = %d'%s
            ci = self.calohits(l, s)
            #print 'calo info: ',ci
            if len(ci)==0:
                return { } # try other lines
            cluster[counter+1] = (ci,self.point)
            #print 'cluster made from: ',ci,self.point
        return cluster
