from math import sqrt, atan2, sin, cos, pi, isnan
from numpy import array,linalg

class tangent_lists(object):
    """
    This object holds the full list of tangent points in tuples
    given two arbitrary cylinder objects. The tangent points result 
    from connecting straight lines to the ring-projections of the 
    cylinders and calculating the maximum four points each where
    such tangent lines can touch both rings simultaneously.
    """
    def __init__(self,cyl1=None,cyl2=None):
        """
        Parameter: takes two geiger cylinder, cyl1, cyl2
        Return: coordinatelist1, coordinatelist2 = 
                    list of dublets (tuple type);
                same for errorlist1, errorlist2
                each dublet contains coordinates (errors of) for one 
                tangent point
                in x, y with z fixed at cylinder value by default.
                These tangent points make up potential hit coordinates
                for each wire with respect to one(!) neighbour wire.
                They will depend on which two wires are compared.
                The resulting lists from this class need to be managed
                from where they have been requested. List 1 is 
                for the cylinder 1, list 2 for cylinder 2
        """
        self._cyl1 = cyl1
        self._cyl2 = cyl2
        self.coordinatelist1 = []     # coordinates on cylinder 1
        self.coordinatelist2 = []     # coordinates on cylinder 2
        self.errorlist1 = []          # errors of coordinates on cylinder 1
        self.errorlist2 = []          # errors of coordinates on cylinder 2
        self.case_manager(cyl1,cyl2)  # send the cylinders off 
                                      # to the correct case for calculation

    def __str__(self):
        s = "List of tangent points and their errors\n"
        s += "For first geiger cylinder:\n"
        for entry in self.coordinatelist1:
            s += "x=%f, y=%f\n" % entry
        s += "\nErrors:\n"
        for entry in self.errorlist1:
            s += "dr=%f, dphi=%f\n" % entry
        s += "\nSecond geiger cylinder:\n"
        for entry in self.coordinatelist2:
            s += "x=%f, y=%f\n" % entry
        s += "\nErrors:\n"
        for entry in self.errorlist2:
            s += "dr=%f, dphi=%f\n" % entry
        return s

    def case_manager(self,cyl1,cyl2):
        """
        input: cyl1, cyl2, two non-sorted tracker_hit objects
        output: none
        process: decide which case is true in relation to the two circles
                 in order to be able to calculate the tangent points of
                 connecting tangents to both circles in x, y plane
                 The calculators finally fill the answer-list.
        """
        wire1 = array([cyl1.x,cyl1.y])
        wire2 = array([cyl2.x,cyl2.y])
        distance = linalg.norm(wire2 - wire1)
        ring_overlap = (cyl1.r + cyl2.r) - distance 

        if (ring_overlap >= 0):
            #print "Found overlapping rings: %f\n" % ring_overlap
            self.coordinatelist1, self.coordinatelist2 = self.calculate_overlapping_rings(cyl1,cyl2,distance)
#            for coor in self.coordinatelist1:
#                if isnan(coor[0]):
#                    print 'found nan in TP'
#                    print cyl1
#                    print cyl2
            self.errorlist1, self.errorlist2 = self.calculate_overlapping_ringerrors(cyl1,cyl2,distance)

        else:
            self.coordinatelist1, self.coordinatelist2 = self.calculate_four_points(cyl1,cyl2,distance)            
#            for coor in self.coordinatelist1:
#                if isnan(coor[0]):
#                    print 'found nan in TP, non overlap'
#                    print cyl1
#                    print cyl2
            self.errorlist1, self.errorlist2 = self.calculate_four_pointerrors(cyl1,cyl2,distance)
        
    def checkforpi(self,d):
        if (d < -pi):
            d += 2.0*pi
        elif (d > pi):
            d -= 2.0*pi
        return d


    def check_subtraction(self,r,dr):
        if ((r-dr) <= 0.0):
            return r
        else:
            return r - dr


    def size_sorting(self,c1,c2):
        if (c1.r <= c2.r):  # size sorting
            rsmall = c1.r
            ds = c1.sigmar
            rlarge = c2.r
            dl = c2.sigmar
            xa = c1.x
            ya = c1.y
            xb = c2.x
            yb = c2.y
            order = 0
        else:
            rsmall = c2.r
            ds = c2.sigmar
            rlarge = c1.r
            dl = c1.sigmar
            xa = c2.x
            ya = c2.y
            xb = c1.x
            yb = c1.y
            order = 1

        return xa,ya,xb,yb,rsmall,ds,rlarge,dl,order


    def deltaxy_to_deltaphi(self,xyfix,dxdyarr,centrelist):
        """
        Input:
        Fixed tangent point coordinates as tuple - 
        have to be solved first before getting here.
        Their errors in x, y, as +- dxy, i.e. both sided
        and finally the centres of the ring to get the ring 
        segment as the error on the ring
        output: 2 numbers, delta phi for each ring point
        """
        dphi = []
        counter = 0
        for dublet in xyfix:
            p = atan2(dublet[1]-centrelist[counter][1],dublet[0]-centrelist[counter][0])
            xy = dublet + dxdyarr[counter] # increase xy tuple by 
            # dx dy tuple one-sided to get the correct dphi angle
            # and not the double sided interval width
            pplus1 = atan2(xy[1]-centrelist[counter][1],xy[0]-centrelist[counter][0])
            xy = dublet - dxdyarr[counter] # decrease
            pplus2 = atan2(xy[1]-centrelist[counter][1],xy[0]-centrelist[counter][0])

            pdiff1 = self.checkforpi(pplus1-p)
            pdiff2 = self.checkforpi(pplus2-p)
            pdiff = pdiff1+pdiff2
            # meanvalue of angle left and right to p as dx dy is added
            # and subtracted
            dphi.append(abs(pdiff))
            counter += 1

        return dphi


    def calc_onhorizontal(self,xwire,ywire,radius,h,xo,yo):
        xA = xwire - radius*radius/(h*h)*(xwire-xo)
        p2 = (xwire-xo)*radius/(h*h)*sqrt(h*h-radius*radius)
        yA_1 = ywire - radius*radius/(h*h)*(ywire-yo) - p2
        yA_2 = ywire - radius*radius/(h*h)*(ywire-yo) + p2

        pair1 = (xA,yA_1)
        pair2 = (xA,yA_2)
        return pair1, pair2

    def calc_onvertical(self,xwire,ywire,radius,h,xo,yo):
        yA = ywire - radius*radius/(h*h)*(ywire-yo)
        p2 = (ywire-yo)*radius/(h*h)*sqrt(h*h-radius*radius)
        xA_1 = xwire - radius*radius/(h*h)*(xwire-xo) - p2
        xA_2 = xwire - radius*radius/(h*h)*(xwire-xo) + p2

        pair1 = (xA_1,yA)
        pair2 = (xA_2,yA)
        return pair1, pair2


    def pointcalc2(self,xwire,ywire,radius,h,xo,yo):
        """
        Formulae to calculate the tangent points, all four from
        this little function. Just needs different input:
        which wire and ring radius is targeted, the 
        distance: centre to 'line connection crossing the 
        centre line between wire centres' = h, and the 
        origin or point of the tangent lines crossing: xo, yo.
        This routine should cope with singularities like 
        vertical and horizontal lines.
        """
        EPS = 1.0e-3 # define as zero
        # Special cases first:
        if (abs(ywire-yo)<=EPS): #horizontal line
            pair1, pair2 = self.calc_onhorizontal(xwire,ywire,radius,h,xo,yo)

        elif (abs(xwire-xo)<=EPS): #vertical line
            pair1, pair2 = self.calc_onvertical(xwire,ywire,radius,h,xo,yo)

        else: # normal calculation
            p2 = (ywire-yo)*radius/(h*h)*sqrt(h*h - radius*radius)
            xA_1 = xwire - radius*radius/(h*h)*(xwire-xo) + p2
            xA_2 = xwire - radius*radius/(h*h)*(xwire-xo) - p2
            
            denom1 = radius*radius+(xA_1-xwire)*(xwire-xo)
            denom2 = radius*radius+(xA_2-xwire)*(xwire-xo)
            yA_1 = (xwire*ywire*xo - xA_1*(xo*ywire + xwire*(ywire-2*yo)) + xA_1*xA_1*(ywire-yo) + radius*radius*yo - xwire*xwire*yo) / denom1
            yA_2 = (xwire*ywire*xo - xA_2*(xo*ywire + xwire*(ywire-2*yo)) + xA_2*xA_2*(ywire-yo) + radius*radius*yo - xwire*xwire*yo) / denom2
            
            pair1 = (xA_1,yA_1)
            pair2 = (xA_2,yA_2)

        return pair1,pair2


    def intersectioncalc(self,distance,rsmall,rlarge,xa,ya,xb,yb):
        """
        Formulae to calculate the intersection points of two rings
        in case they overlap.
        Input: all ring characteristics
        Output: coordinates of both intersection points in x, y as tuple
        """
        d = sqrt(abs(distance*distance*(2*rlarge*rlarge + 2*rsmall*rsmall - distance*distance) - (rlarge*rlarge - rsmall*rsmall)*(rlarge*rlarge - rsmall*rsmall)))
        p1 = (yb-ya)*d/(2*distance*distance)
        p2 = (xb-xa)*d/(2*distance*distance)

        x1 = 0.5*(xa+xb) - ((xb-xa)*(rlarge*rlarge-rsmall*rsmall))/(2*distance*distance) + p1
        x2 = 0.5*(xa+xb) - ((xb-xa)*(rlarge*rlarge-rsmall*rsmall))/(2*distance*distance) - p1
        y1 = 0.5*(ya+yb) - ((yb-ya)*(rlarge*rlarge-rsmall*rsmall))/(2*distance*distance) - p2
        y2 = 0.5*(ya+yb) - ((yb-ya)*(rlarge*rlarge-rsmall*rsmall))/(2*distance*distance) + p2
        
        pair1 = (x1,y1)
        pair2 = (x2,y2)
        return pair1,pair2


    def calculate_four_points(self,cyl1,cyl2,distance):
        """
        input: cyl1, cyl2, two non-sorted tracker_hit objects
               distance between wire centres for convenience
        output: list of tuples with tangent points in x,y plane
        """
        # non-overlapping rings still can have equal radii - capture this
        # and get the points from different functions than overlapping
        # rings.
        if (abs(cyl1.r - cyl2.r) < 0.01):
#            print "Found equal radius rings without overlap"
            overlapflag = 0
            resultA, resultB = self.equal_rings(cyl1,cyl2,distance,overlapflag)
            return resultA, resultB


        xa,ya,xb,yb,rsmall,ds,rlarge,dl,order = self.size_sorting(cyl1,cyl2)
        # getting to the first, smaller  ring
        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        pair1, pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        resultA = [pair1,pair2]

        # second, larger ring, same origin
        h = rlarge * distance / (rlarge - rsmall)
        pair1, pair2 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        resultB = [pair1,pair2]

        # origin between rings
        ha = rsmall * distance / (rlarge + rsmall)
        xorigin = (1-ha/distance)*xa + ha/distance*xb
        yorigin = (1-ha/distance)*ya + ha/distance*yb
        pair1, pair2 = self.pointcalc2(xa,ya,rsmall,ha,xorigin,yorigin)
        resultA.append(pair1)
        resultA.append(pair2)

        # same origin between rings, for larger ring
        hb = rlarge * distance / (rlarge + rsmall)
        pair1, pair2 = self.pointcalc2(xb,yb,rlarge,hb,xorigin,yorigin)
        resultB.append(pair1)
        resultB.append(pair2)

        if (order<1):
            return resultA, resultB
        else:
            return resultB, resultA


    def calculate_four_pointerrors(self,cyl1,cyl2,distance):
        """
        input: cyl1, cyl2, two non-sorted tracker_hit objects
               distance between wire centres for convenience
        output: list of tuples with errors on tangent points 
                in the polar r, phi plane, i.e. around the ring
        """
        xa,ya,xb,yb,rstemp,dr_small,rltemp,dr_large,order = self.size_sorting(cyl1,cyl2)
        
        swapped = False
        # make sure that the origin lands on the correct side
        if ((rstemp+dr_small) >= (rltemp-dr_large)):
            rlarge = rstemp + dr_small
            rsmall = self.check_subtraction(rltemp,dr_large)
            # swap centres
	    xd = xa
	    xa = xb
	    xb = xd
	    yd = ya
	    ya = yb
	    yb = yd
#            print "SWAPPED in point errors"
            swapped = True
        else:
            rsmall = rstemp + dr_small
            rlarge = self.check_subtraction(rltemp,dr_large)
            
#        order = order is not swapped
#        print "order is: %d" % order
        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        l_pair1, l_pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        l1 = array(l_pair1)
        l2 = array(l_pair2)
        h = rlarge * distance / (rlarge - rsmall)
        l_pair3, l_pair4 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        l3 = array(l_pair3)
        l4 = array(l_pair4)
#        print l1,l2,l3,l4

        rsmall = self.check_subtraction(rstemp,dr_small)
        rlarge = rltemp + dr_large
        if (swapped):
            # swap centres back to initial
            xd = xa
            xa = xb
            xb = xd
            yd = ya
            ya = yb
            yb = yd

        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        s_pair1, s_pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        s1 = array(s_pair1)
        s2 = array(s_pair2)
        h = rlarge * distance / (rlarge - rsmall)
        s_pair3, s_pair4 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        s3 = array(s_pair3)
        s4 = array(s_pair4)
#        print s1,s2,s3,s4
        
        if (order):
            fix1 = self.coordinatelist2[0] # the tangent point coordinates
            fix2 = self.coordinatelist2[1] # the tangent point coordinates
            fix3 = self.coordinatelist1[0] # the tangent point coordinates
            fix4 = self.coordinatelist1[1] # the tangent point coordinates
            if (swapped):
                pair1 = l1 - s4 # first point on ring, delta x, delta y here
                pair2 = l2 - s3 # second point on ring, delta x, delta y here
                pair3 = l3 - s2 # first point on ring, delta x, delta y here
                pair4 = l4 - s1 # second point on ring, delta x, delta y here
            else:
                pair1 = l1 - s1 # first point on ring, delta x, delta y here
                pair2 = l2 - s2 # second point on ring, delta x, delta y here
                pair3 = l3 - s3 # first point on ring, delta x, delta y here
                pair4 = l4 - s4 # second point on ring, delta x, delta y here
        else:
            fix1 = self.coordinatelist1[0] # the tangent point coordinates
            fix2 = self.coordinatelist1[1] # the tangent point coordinates
            fix3 = self.coordinatelist2[0] # the tangent point coordinates
            fix4 = self.coordinatelist2[1] # the tangent point coordinates
            if (swapped):
                pair1 = l1 - s4 # first point on ring, delta x, delta y here
                pair2 = l2 - s3 # second point on ring, delta x, delta y here
                pair3 = l3 - s2 # first point on ring, delta x, delta y here
                pair4 = l4 - s1 # second point on ring, delta x, delta y here
            else:
                pair1 = l1 - s1 # first point on ring, delta x, delta y here
                pair2 = l2 - s2 # second point on ring, delta x, delta y here
                pair3 = l3 - s3 # first point on ring, delta x, delta y here
                pair4 = l4 - s4 # second point on ring, delta x, delta y here

        centrelist = [(xa,ya),(xa,ya),(xb,yb),(xb,yb)]
        xy = [fix1,fix2,fix3,fix4]
        xyfix = array(xy)
        dxdy = [pair1,pair2,pair3,pair4]
        dxdyarr = array(dxdy)
        dphi = self.deltaxy_to_deltaphi(xyfix,dxdyarr,centrelist)

        resultA = [(dr_small,dphi[0]),(dr_small,dphi[1])]
        resultB = [(dr_large,dphi[2]),(dr_large,dphi[3])]

#        print "Between rings now in point errors"
        # origin between rings
        # here vary ring radii in the same direction to get 
        # the error variation on the rings - change of symmetry 
        # compared to outer lines
        rsmall = rstemp + dr_small      # larger ring radius
        rlarge = rltemp + dr_large      # larger ring radius
        ring_overlap = (rlarge + rsmall) - distance 
        if (ring_overlap >= 0):
#            print "Overlapping on larger errors"
            # don't get between rings when they overlap
            # the two points are taken as the intersection
            # points of the two rings
            l_pair1, l_pair2 = self.intersectioncalc(distance,rsmall,rlarge,xa,ya,xb,yb)
            l1 = array(l_pair1)
            l2 = array(l_pair2)
            l3 = array(l_pair2)
            l4 = array(l_pair1)
#            print l1,l2,l3,l4
        else:        
#            print "Not overlapping on larger errors"
            ha = rsmall * distance / (rlarge + rsmall)
            xorigin = (1-ha/distance)*xa + ha/distance*xb
            yorigin = (1-ha/distance)*ya + ha/distance*yb
            l_pair1, l_pair2 = self.pointcalc2(xa,ya,rsmall,ha,xorigin,yorigin)
            l1 = array(l_pair1)
            l2 = array(l_pair2)
            hb = rlarge * distance / (rlarge + rsmall)
            l_pair3, l_pair4 = self.pointcalc2(xb,yb,rlarge,hb,xorigin,yorigin)
            l3 = array(l_pair3)
            l4 = array(l_pair4)
#            print l1,l2,l3,l4

        # no need to check here for overlap
        rsmall = self.check_subtraction(rstemp,dr_small) # smaller ring radius
        rlarge = self.check_subtraction(rltemp,dr_large) # smaller ring radius
        ha = rsmall * distance / (rlarge + rsmall)
        xorigin = (1-ha/distance)*xa + ha/distance*xb
        yorigin = (1-ha/distance)*ya + ha/distance*yb
        s_pair1, s_pair2 = self.pointcalc2(xa,ya,rsmall,ha,xorigin,yorigin)
        s1 = array(s_pair1)
        s2 = array(s_pair2)
        hb = rlarge * distance / (rlarge + rsmall)
        s_pair3, s_pair4 = self.pointcalc2(xb,yb,rlarge,hb,xorigin,yorigin)
        s3 = array(s_pair3)
        s4 = array(s_pair4)
#        print s1,s2,s3,s4

        centrelist = [(xa,ya),(xa,ya),(xb,yb),(xb,yb)]
        if (order):
            fix1 = self.coordinatelist2[2] # the tangent point coordinates
            pair1 = l1 - s1 # first point on ring, delta x, delta y here
            fix2 = self.coordinatelist2[3] # the tangent point coordinates
            pair2 = l2 - s2 # second point on ring, delta x, delta y here
            fix3 = self.coordinatelist1[2] # the tangent point coordinates
            pair3 = l3 - s3 # first point on ring, delta x, delta y here
            fix4 = self.coordinatelist1[3] # the tangent point coordinates
            pair4 = l4 - s4 # second point on ring, delta x, delta y here
        else:
            fix1 = self.coordinatelist1[2] # the tangent point coordinates
            pair1 = l1 - s1 # first point on ring, delta x, delta y here
            fix2 = self.coordinatelist1[3] # the tangent point coordinates
            pair2 = l2 - s2 # second point on ring, delta x, delta y here
            fix3 = self.coordinatelist2[2] # the tangent point coordinates
            pair3 = l3 - s3 # first point on ring, delta x, delta y here
            fix4 = self.coordinatelist2[3] # the tangent point coordinates
            pair4 = l4 - s4 # second point on ring, delta x, delta y here

        xy = [fix1,fix2,fix3,fix4]
        xyfix = array(xy)
        dxdy = [pair1,pair2,pair3,pair4]
        dxdyarr = array(dxdy)
        dphi = self.deltaxy_to_deltaphi(xyfix,dxdyarr,centrelist)

        resultA.append((dr_small,dphi[0]))
        resultA.append((dr_small,dphi[1]))
        resultB.append((dr_large,dphi[2]))
        resultB.append((dr_large,dphi[3]))

        if (order<1):
            return resultA, resultB
        else:
            return resultB, resultA


    def calculate_overlapping_rings(self,cyl1,cyl2,distance):
        """
        input: cyl1, cyl2, two non-sorted tracker_hit objects
               distance between wire centres for convenience
        output: list of tuples with tangent points in x,y plane
        """
        # overlapping rings can have equal radii - capture this
        # and get the points from different functions than overlapping
        # rings.
        if (abs(cyl1.r - cyl2.r) < 0.01):
#            print "Found equal radius rings with overlap"
            overlapflag = 1
            resultA, resultB = self.equal_rings(cyl1,cyl2,distance,overlapflag)
            return resultA, resultB


        xa,ya,xb,yb,rsmall,ds,rlarge,dl,order = self.size_sorting(cyl1,cyl2)
        
        # getting to the first, smaller  ring
        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        pair1, pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        resultA = [pair1,pair2]

        # second, larger ring, same origin
        h = rlarge * distance / (rlarge - rsmall)
        pair1, pair2 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        resultB = [pair1,pair2]

        # the final two points are taken as the intersection
        # points of the two rings
        pair1, pair2 = self.intersectioncalc(distance,rsmall,rlarge,xa,ya,xb,yb)

        resultA.append(pair1)
        resultA.append(pair2)
        resultB.append(pair1)
        resultB.append(pair2)

        if (order<1):
            return resultA, resultB
        else:
            return resultB, resultA


    def calculate_overlapping_ringerrors(self,cyl1,cyl2,distance):
        """
        input: cyl1, cyl2, two non-sorted tracker_hit objects
               distance between wire centres for convenience
        output: list of tuples with errors on tangent points 
                in the polar r, phi plane, i.e. around the ring
        """
        xa,ya,xb,yb,rstemp,dr_small,rltemp,dr_large,order = self.size_sorting(cyl1,cyl2)

        swapped = False
        # make sure that the origin lands on the correct side
        if ((rstemp+dr_small) > (rltemp-dr_large)):
            rlarge = rstemp + dr_small
            rsmall = self.check_subtraction(rltemp,dr_large)
            # swap centres
	    xd = xa
	    xa = xb
	    xb = xd
	    yd = ya
	    ya = yb
	    yb = yd
            #print "SWAPPED in overlapping ring errors"
            swapped = True
        else:
            rsmall = rstemp + dr_small
            rlarge = self.check_subtraction(rltemp,dr_large)
            

#        order = order is not swapped
        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        l_pair1, l_pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        l1 = array(l_pair1)
        l2 = array(l_pair2)
        h = rlarge * distance / (rlarge - rsmall)
        l_pair3, l_pair4 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        l3 = array(l_pair3)
        l4 = array(l_pair4)
        #print 'into intersection calc: ',distance,rsmall,rlarge,xa,ya,xb,yb
        l_pair5, l_pair6 = self.intersectioncalc(distance,rsmall,rlarge,xa,ya,xb,yb)
        l5 = array(l_pair5)
        l6 = array(l_pair6)

        rsmall = self.check_subtraction(rstemp,dr_small)
        rlarge = rltemp + dr_large
        if (swapped):
            # swap centres back to initial
            xd = xa
            xa = xb
            xb = xd
            yd = ya
            ya = yb
            yb = yd
            
        h = rsmall * distance / (rlarge - rsmall)
        xorigin = (1+h/distance)*xa - h/distance*xb
        yorigin = (1+h/distance)*ya - h/distance*yb
        s_pair1, s_pair2 = self.pointcalc2(xa,ya,rsmall,h,xorigin,yorigin)
        s1 = array(s_pair1)
        s2 = array(s_pair2)
        h = rlarge * distance / (rlarge - rsmall)
        s_pair3, s_pair4 = self.pointcalc2(xb,yb,rlarge,h,xorigin,yorigin)
        s3 = array(s_pair3)
        s4 = array(s_pair4)
        s_pair5, s_pair6 = self.intersectioncalc(distance,rsmall,rlarge,xa,ya,xb,yb)
        s5 = array(s_pair5)
        s6 = array(s_pair6)


        # Caution: swapping origins swaps list order
        # not to compare the wrong points to each other!
        if (swapped):
            fix1 = self.coordinatelist2[0] # the tangent point coordinates
            pair1 = l1 - s4 # first point on ring, delta x, delta y here
            fix2 = self.coordinatelist2[1] # the tangent point coordinates
            pair2 = l2 - s3 # second point on ring, delta x, delta y here
            fix3 = self.coordinatelist1[0] # the tangent point coordinates
            pair3 = l3 - s2 # first point on ring, delta x, delta y here
            fix4 = self.coordinatelist1[1] # the tangent point coordinates
            pair4 = l4 - s1 # second point on ring, delta x, delta y here
            fix5 = self.coordinatelist1[3] # the tangent point coordinates
            pair5 = l5 - s6 # first intersection point on ring
            fix6 = self.coordinatelist1[2] # the tangent point coordinates
            pair6 = l6 - s5 # second intersection point on ring
            centrelist = [(xb,yb),(xb,yb),(xa,ya),(xa,ya),(xa,ya),(xa,ya)]
        else:
            fix1 = self.coordinatelist1[0] # the tangent point coordinates
            pair1 = l1 - s1 # first point on ring, delta x, delta y here
            fix2 = self.coordinatelist1[1] # the tangent point coordinates
            pair2 = l2 - s2 # second point on ring, delta x, delta y here
            fix3 = self.coordinatelist2[0] # the tangent point coordinates
            pair3 = l3 - s3 # first point on ring, delta x, delta y here
            fix4 = self.coordinatelist2[1] # the tangent point coordinates
            pair4 = l4 - s4 # second point on ring, delta x, delta y here
            fix5 = self.coordinatelist1[2] # the tangent point coordinates
            pair5 = l5 - s5 # first intersection point on ring
            fix6 = self.coordinatelist1[3] # the tangent point coordinates
            pair6 = l6 - s6 # second intersection point on ring
            centrelist = [(xa,ya),(xa,ya),(xb,yb),(xb,yb),(xa,ya),(xa,ya)]

        # intersection points are common to both rings by definition
        xy = [fix1,fix2,fix3,fix4,fix5,fix6]
        xyfix = array(xy)
        dxdy = [pair1,pair2,pair3,pair4,pair5,pair6]
        dxdyarr = array(dxdy)
        dphi = self.deltaxy_to_deltaphi(xyfix,dxdyarr,centrelist)

        resultA = [(dr_small,dphi[0]),(dr_small,dphi[1]),(dr_small,dphi[4]),(dr_small,dphi[5])]
        resultB = [(dr_large,dphi[2]),(dr_large,dphi[3]),(dr_large,dphi[4]),(dr_large,dphi[5])]

        if (order<1):
            return resultA, resultB
        else:
            return resultB, resultA


    def equal_rings(self,cyl1,cyl2,distance,overlapflag):
        """
        input: cyl1, cyl2, two similar to equal radius tracker_hit objects
               distance between wire centres for convenience
        output: list of tuples with tangent points in x,y plane
        """
        xa,ya,xb,yb,rsmall,ds,rlarge,dl,order = self.size_sorting(cyl1,cyl2)
        
        # for left and right tangents, shift centre connecting 
        # line to ring radius - tangents parallel to centre
        # line by definition
        alpha = atan2((yb-ya),(xb-xa)) #can handle singularity - useful here

        x1 = xa + rsmall*sin(alpha)
        x2 = xa - rsmall*sin(alpha)
        y1 = ya - rsmall*cos(alpha)
        y2 = ya + rsmall*cos(alpha)
        resultA = [(x1,y1),(x2,y2)]

        x1 = xb + rlarge*sin(alpha)
        x2 = xb - rlarge*sin(alpha)
        y1 = yb - rlarge*cos(alpha)
        y2 = yb + rlarge*cos(alpha)
        resultB = [(x1,y1),(x2,y2)]

        if (overlapflag < 1):
            # origin between rings
            ha = rsmall * distance / (rlarge + rsmall)
            xorigin = (1-ha/distance)*xa + ha/distance*xb
            yorigin = (1-ha/distance)*ya + ha/distance*yb
            pair1, pair2 = self.pointcalc2(xa,ya,rsmall,ha,xorigin,yorigin)
            resultA.append(pair1)
            resultA.append(pair2)
            
            # same origin between rings, for larger ring
            hb = rlarge * distance / (rlarge + rsmall)
            pair1, pair2 = self.pointcalc2(xb,yb,rlarge,hb,xorigin,yorigin)
            resultB.append(pair1)
            resultB.append(pair2)

            if (order<1):
                return resultA, resultB
            else:
                return resultB, resultA

        else: # don't get between rings when they overlap
            # the final two points are taken as the intersection
            # points of the two rings
            pair1, pair2 = self.intersectioncalc(distance,rsmall,rlarge,xa,ya,xb,yb)
            resultA.append(pair1)
            resultA.append(pair2)
            resultB.append(pair1)
            resultB.append(pair2)

            if (order<1):
                return resultA, resultB
            else:
                return resultB, resultA
