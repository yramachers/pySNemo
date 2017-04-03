from math import atan2, cos, sin, pi, isnan
import pysnemo.utility.interval as I

class tp_controller(object):
    """
    This object decides how many tangent points are actually distinct
    points or identical to the precision of their respective errors.
    Input: tangent_lists object and sigma factor for error interval width
    Process: filtered coordinate and error lists of tuples
    """
    def __init__(self,tp_lists = None, sigmafactor = 1):
        self.filtered_coordinates1 = []   # distinct coordinates on 1
        self.filtered_coordinates2 = []   # distinct coordinates on 2
        self.filtered_errors1 = []    # errors of coordinates on cylinder 1
        self.filtered_errors2 = []    # errors of coordinates on cylinder 2
        self.pairs = []               # point to point pairs for connections
        self.pair_errors = []         # point to point errors
        self.idlist = []              # point to point id (counter,cyl_id)
        self.process_lists(tp_lists,sigmafactor)   # send the cylinders off 


    def __str__(self):
        s = "Filtered list of tangent points and their errors\n"
        s += "For first geiger cylinder:\n"
        for entry in self.filtered_coordinates1:
            s += "x=%f, y=%f\n" % entry
        s += "\nErrors:\n"
        for entry in self.filtered_errors1:
            s += "dr=%f, dphi=%f\n" % entry
        s += "\nSecond geiger cylinder:\n"
        for entry in self.filtered_coordinates2:
            s += "x=%f, y=%f\n" % entry
        s += "\nErrors:\n"
        for entry in self.filtered_errors2:
            s += "dr=%f, dphi=%f\n" % entry
        return s


    def checkforpi(self,d):
        if (d < -pi):
            d += 2.0*pi
        elif (d > pi):
            d -= 2.0*pi
        return d

    def checkforpihalf(self,d):
        sign = 1
        if (d < -pi/2.0):
            d += pi
            sign = -1
        elif (d > pi/2.0):
            d -= pi
            sign = -1
        return sign,d

    def check_overlaps(self,tp,inval1,inval2):
        """
        return list of index tuples of overlapping intervals in 
        lists inval1 and inval2 corresponding to intervals on
        ring1 and ring2. Knows the non-overlapping points, i.e.
        copies those to result lists in here.
        """
        c1 = tp.coordinatelist1
        c2 = tp.coordinatelist2
        err1 = tp.errorlist1
        err2 = tp.errorlist2
        id1 = tp._cyl1.meta_info[0] # geiger id
        id2 = tp._cyl2.meta_info[0] # geiger id

        ovlist1 = []
        ovlist2 = []
        # Ring1
        if (inval1[0].angle_overlap(inval1[1])):
            ovlist1.append((0,1))
        else:
            self.filtered_coordinates1.append(c1[0])
            self.filtered_errors1.append(err1[0])
            self.filtered_coordinates1.append(c1[1])
            self.filtered_errors1.append(err1[1])
        if (inval1[2].angle_overlap(inval1[3])):
            ovlist1.append((2,3))
        else:
            self.filtered_coordinates1.append(c1[2])
            self.filtered_errors1.append(err1[2])
            self.filtered_coordinates1.append(c1[3])
            self.filtered_errors1.append(err1[3])

        # Ring2
        if (inval2[0].angle_overlap(inval2[1])):
            ovlist2.append((0,1))
        else:
            self.filtered_coordinates2.append(c2[0])
            self.filtered_errors2.append(err2[0])
            self.filtered_coordinates2.append(c2[1])
            self.filtered_errors2.append(err2[1])
        if (inval2[2].angle_overlap(inval2[3])):
            ovlist2.append((2,3))
        else:
            self.filtered_coordinates2.append(c2[2])
            self.filtered_errors2.append(err2[2])
            self.filtered_coordinates2.append(c2[3])
            self.filtered_errors2.append(err2[3])
        
        if (len(ovlist1)<1 and len(ovlist2)<1):
            for (t1,t2) in zip(self.filtered_coordinates1,self.filtered_coordinates2):
                self.pairs.append((t1,t2))
                self.idlist.append(id1)
                self.idlist.append(id2)
                #print 'in checkoverlap: pair/ids: ',(t1,t2),id1,id2
            for (t1,t2) in zip(self.filtered_errors1,self.filtered_errors2):
                self.pair_errors.append((t1,t2))

        return ovlist1, ovlist2


    def process_lists(self,tp_lists,sigma):
        """
        Loop through lists and check whether points with errors 
        overlap to be combined into single points at mean positions
        and combined error intervall.
        - The result forms the final list of tuples of valid hits in 
        the tracker with associated errors.
        - Only call for FOUR tangent points in lists - makes no sense
        for fewer tangent points identified earlier, see case manager
        method in tangent_lists object.
        - Fills the attribute lists from constructor: filtered_...
        """
        c1 = tp_lists.coordinatelist1
        radius1 = tp_lists._cyl1.r
        id1 = tp_lists._cyl1.meta_info[0] # geiger id
        c2 = tp_lists.coordinatelist2
        radius2 = tp_lists._cyl2.r
        id2 = tp_lists._cyl2.meta_info[0] # geiger id
        err1 = tp_lists.errorlist1
        err2 = tp_lists.errorlist2
        phi1 = []
        phi2 = []
        for dublet in c1:
#            if isnan(dublet[0]) or isnan(dublet[1]):
#                print 'found nan'
            phi1.append(atan2(dublet[1]-tp_lists._cyl1.y,dublet[0]-tp_lists._cyl1.x)) # can handle an x=0 case 
        for dublet in c2:
            phi2.append(atan2(dublet[1]-tp_lists._cyl2.y,dublet[0]-tp_lists._cyl2.x)) # can handle an x=0 case 

        counter = 0
        # list of intervals, at sigma error widths
        inval1 = []
        inval2 = []
        for dublet in err1:
            diff = phi1[counter]-sigma*dublet[1]
            diffstart = self.checkforpi(diff)
            diff = phi1[counter]+sigma*dublet[1]
            diffend = self.checkforpi(diff)
            if (diffend<diffstart):
                inval1.append(I.Interval(diffend,diffstart))
            else:
                inval1.append(I.Interval(diffstart,diffend))
            counter += 1
        counter = 0
        for dublet in err2:
            diff = phi2[counter]-sigma*dublet[1]
            diffstart = self.checkforpi(diff)
            diff = phi2[counter]+sigma*dublet[1]
            diffend = self.checkforpi(diff)
            if (diffend<diffstart):
                inval2.append(I.Interval(diffend,diffstart))
            else: 
                inval2.append(I.Interval(diffstart,diffend))
            counter += 1

        # Now check for the overlaps of error intervals
        ovlist1,ovlist2 = self.check_overlaps(tp_lists,inval1,inval2)

        for entry in ovlist1:
            i = entry[0]
            j = entry[1]
            #print "found overlap"
            h = inval1[i].hull(inval1[j]) # combined interval when overlap
            sign1,dstart = self.checkforpihalf(h.start)
            sign2,dend = self.checkforpihalf(h.end)
            if (sign1<0 and sign2<0): 
                if (dend<=dstart):
                    h = I.Interval(dend,dstart)
                else:
                    h = I.Interval(dstart,dend)
                signx = -1
            else: 
                signx = 1
            centre = h.midinterval()
            xc = signx*radius1*cos(centre)
            yc = radius1*sin(centre)

            # set new centre point as hit coordinate
            # and errors, dr as before, and dphi as in phi+-dphi
            self.filtered_coordinates1.insert(i,(xc + tp_lists._cyl1.x,yc + tp_lists._cyl1.y))
            self.filtered_errors1.insert(i,(err1[0][0],h.end-centre))

        for entry in ovlist2:
            i = entry[0]
            j = entry[1]
            #print "found overlap 2"
            h = inval2[i].hull(inval2[j]) # combined interval when overlap
            sign1,dstart = self.checkforpihalf(h.start)
            sign2,dend = self.checkforpihalf(h.end)
            if (sign1<0 and sign2<0): 
                if (dend<=dstart):
                    h = I.Interval(dend,dstart)
                else:
                    h = I.Interval(dstart,dend)
                signx = -1
            else: 
                signx = 1
            centre = h.midinterval()
            xc = signx*radius2*cos(centre)
            yc = radius2*sin(centre)

            # set new centre point as hit coordinate
            # and errors, dr as before, and dphi as in phi+-dphi
            self.filtered_coordinates2.insert(i,(xc + tp_lists._cyl2.x,yc + tp_lists._cyl2.y))
            self.filtered_errors2.insert(i,(err2[0][0],h.end-centre))

        # Now build pairs in case it hasn't happened already
        if (len(self.pairs)<1):
            # fill completely first
            for (entry1,entry2) in zip(self.filtered_coordinates1,self.filtered_coordinates2):
                self.pairs.append((entry1,entry2))
                self.idlist.append(id1)
                self.idlist.append(id2)
                #print 'no pairs yet: pair/ids: ',(entry1,entry2),id1,id2
            for entry1,entry2 in zip(self.filtered_errors1,self.filtered_errors2):
                self.pair_errors.append((entry1,entry2))

            # final permutations; both empty and equal are done
            if (len(ovlist1)<1 and len(ovlist2)==1): # ovlist1 is empty
                for entry in ovlist2:
                    self.pairs.insert(entry[1],(self.filtered_coordinates1[entry[1]],self.filtered_coordinates2[entry[0]]))
                    self.pair_errors.insert(entry[1],(self.filtered_errors1[entry[1]],self.filtered_errors2[entry[0]]))
                    self.idlist.insert(entry[1],id1)
                    self.idlist.insert(entry[1]+1,id2)

            # vice versa
            if (len(ovlist2)<1 and len(ovlist1)==1): # ovlist2 is empty
                for entry in ovlist1:
                    self.pairs.insert(entry[1],(self.filtered_coordinates1[entry[0]],self.filtered_coordinates2[entry[1]]))
                    self.pair_errors.insert(entry[1],(self.filtered_errors1[entry[0]],self.filtered_errors2[entry[1]]))
                    self.idlist.insert(entry[1],id1)
                    self.idlist.insert(entry[1]+1,id2)

            # double overlap on one ring = one p to four p on second ring
            if (len(ovlist1)<1 and len(ovlist2)==2): # ovlist1 is empty
                for entry in ovlist2:
                    if entry[1]>1:
                        self.pairs.insert(entry[0],(self.filtered_coordinates1[entry[0]],self.filtered_coordinates2[0]))
                        self.pairs.insert(entry[1],(self.filtered_coordinates1[entry[1]],self.filtered_coordinates2[1]))
                        self.pair_errors.insert(entry[0],(self.filtered_errors1[entry[0]],self.filtered_errors2[0]))
                        self.pair_errors.insert(entry[1],(self.filtered_errors1[entry[1]],self.filtered_errors2[1]))
                    self.idlist.insert(entry[0],id1)
                    self.idlist.insert(entry[1],id2)

            # vice versa
            if (len(ovlist2)<1 and len(ovlist1)==2): # ovlist2 is empty
                for entry in ovlist1:
                    if entry[1]>1:
                        self.pairs.insert(entry[0],(self.filtered_coordinates1[0],self.filtered_coordinates2[entry[0]]))
                        self.pairs.insert(entry[1],(self.filtered_coordinates1[1],self.filtered_coordinates2[entry[1]]))
                        self.pair_errors.insert(entry[0],(self.filtered_errors1[0],self.filtered_errors2[entry[0]]))
                        self.pair_errors.insert(entry[1],(self.filtered_errors1[1],self.filtered_errors2[entry[1]]))
                        self.idlist.insert(entry[1],id1)
                        self.idlist.insert(entry[0],id2)
                    
            # Both are non-empty but not equal length
            if (len(ovlist1)==1 and len(ovlist2)>len(ovlist1)):
                tt = ovlist1[0]
                for entry in ovlist2:
                    if (entry[0] != tt[0]):
                        self.pairs.insert(entry[0],(self.filtered_coordinates1[entry[1]-1],self.filtered_coordinates2[entry[0]-1]))
                        self.pair_errors.insert(entry[0],(self.filtered_errors1[entry[1]-1],self.filtered_errors2[entry[0]-1]))
                        self.idlist.insert(entry[0],id1)
                        self.idlist.insert(entry[0]+1,id2)
                        
            if (len(ovlist2)==1 and len(ovlist1)>len(ovlist2)):
                tt = ovlist2[0]
                for entry in ovlist1:
                    if (entry[0] != tt[0]):
                        self.pairs.insert(entry[0],(self.filtered_coordinates1[entry[0]-1],self.filtered_coordinates2[entry[1]-1]))
                        self.pair_errors.insert(entry[0],(self.filtered_errors1[entry[0]-1],self.filtered_errors2[entry[1]-1]))
                        self.idlist.insert(entry[0],id1)
                        self.idlist.insert(entry[0]+1,id2)

        # Special case of identical pair points on overlapping rings
        eps = 1.0e-6 # small number
        for counter in xrange(0,len(self.pairs)-1):
            (p1,p2),(p3,p4) = self.pairs[counter]
            (p1e,p2e),(p3e,p4e) = self.pair_errors[counter]
            (id1,id2) = (self.idlist[counter],self.idlist[counter+1])
            length = abs((p1-p3)*(p1-p3)+(p2-p4)*(p2-p4))
            if (length <= eps):
                # swap pair partners
                (p1sw,p2sw),(p3sw,p4sw) = self.pairs[counter+1]
                (p1swe,p2swe),(p3swe,p4swe) = self.pair_errors[counter+1]
                (id1sw,id2sw) = (self.idlist[counter+2],self.idlist[counter+3])
                newp = ((p1,p2),(p3sw,p4sw))
                newid = (id1,id2sw)
                self.pairs[counter] = newp
                (self.idlist[counter],self.idlist[counter+1]) = newid
                newp = ((p1sw,p2sw),(p3,p4))
                newid = (id1sw,id2)
                self.pairs[counter+1] = newp
                (self.idlist[counter+2],self.idlist[counter+3]) = newid

                newp = ((p1e,p2e),(p3swe,p4swe))
                self.pair_errors[counter] = newp
                newp = ((p1swe,p2swe),(p3e,p4e))
                self.pair_errors[counter+1] = newp

