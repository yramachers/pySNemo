from math import sqrt
from numpy import array, where, abs
from pysnemo.io.edm import tracker_hit
import pysnemo.utility.tangentpoints as TP
import pysnemo.utility.tpcontroller as TC
import numpy as np

class raw_hits(object):
    """
    Raw Geiger Hits object
    Input: get the grid constant, i.e. the distance between wires
           and secondly the list of geiger cylinders

    Provides: make_hitdictionary() method with

    Output:   dict of actual hits for each ring: format is
        key = id number taken from input list sequences
        value is a list of triplets: 2 cylinder objects to make a pair
        and one resulting tp_controller object holding the 
        actual filtered hits on both rings in that pair (and errors)
    """
    def __init__(self, sigma, lattice=4.4, data = []):
        self.data = data
        self.lattice = lattice
        self.sigma = sigma
        
        #holds the dictionary on who is neighbour to whom
        self.neighbours = self.fill_nearest()
        #print self.neighbours
        self.hd = {} # full hits dictionary
        
    def __str__(self):
        s = "Each wire from the list of wire coordinates:\n"
        s += str(self.centres)
        s += "\nhas neighbours in the lattice_constant vicinity,\n"
        s += "in dictionary format\n"
        s += "{tuple(x,y,z,r,id):list of neighbours (x,y,z,r)}:\n"
        for k,v in self.neighbours.iteritems():
            s += "Key:%s\n"% str(k)
            s += "List of neighbour circles: %s\n" % str(v)
        return s

    def norm(self,t1,t2):
        # compute distance between x,y tuple objects
        temp = (t1[0]-t2[0])*(t1[0]-t2[0]) + (t1[1]-t2[1])*(t1[1]-t2[1])
        return sqrt(temp)
        
    def fill_nearest(self):
        """
        Get the centre points one by one and find the nearest wires to each.
        Then put those found into the dictionary
        """
        # split the tracker_hit objects for key construction
        self.wires = []
        self.radii = []
        self.dr = []
        self.dz = []
        self.mi = []
        for gcyl in self.data:
            self.wires.append((gcyl.x,gcyl.y,gcyl.z))
            self.radii.append(gcyl.r)
            self.dr.append(gcyl.sigmar)
            self.dz.append(gcyl.sigmaz)
            self.mi.append(gcyl.meta_info)

        nn = {}
        mindist = sqrt(2.0)*self.lattice + 0.1*self.lattice # diagonal+10%
        rel_distance = []
        templist = []
        indexlist = []
        for counter,entry in enumerate(self.wires):
            # want a tuple, the counter is the id, identical to input list
            # of wire centres and radii, carries over to results dict
            key = (entry[0],entry[1],entry[2],self.dz[counter],self.radii[counter],self.dr[counter],self.mi[counter],counter)
            for t in self.wires:
                distance = self.norm(t,entry)
                rel_distance.append(distance)
            for i in xrange(len(rel_distance)):
                if (rel_distance[i]<=mindist and rel_distance[i]>self.dr[i]):
                    indexlist.append(i)
            # book results in format for geigerzylinder
            # that is (xwire, ywire, zwire, radius) triplet+
            for i in indexlist: 
                #tup = (self.wires[i][0],self.wires[i][1],self.wires[i][2],self.dz[i],self.radii[i],self.dr[i])
                templist.append(self.data[i]) # book full tracker_hit
            nn[key] = templist
            templist = []
            indexlist = []
            rel_distance = []
        return nn

    
    def make_hitdictionary(self):
        """
        Results in pairs of geiger cell objects for analysis of tangent points
        Radial error not taken into account directly yet -> can have 
        small radial discrepancies between tangent points.
        """
        hitsdict = {}
        resultlist = []
        for k, v in self.neighbours.iteritems():
            # return key to tracker_hit
            c1 = tracker_hit(k[0],k[1],k[2],k[3],k[4],k[5])
            minfo = k[6] # a tuple
            c1.set_info(minfo[0],minfo[1],minfo[2],minfo[3],minfo[4],minfo[5])
            idnumber = k[7]
            for entry in v:
                #c2 = tracker_hit(entry[0],entry[1],entry[2],entry[3],entry[4],entry[5])
                c2 = entry # is a tracker_hit already
                tp = TP.tangent_lists(c1,c2)
                tpc = TC.tp_controller(tp,self.sigma) # full list of hits
                resultlist.append((c1,c2,tpc))
            if (len(v)):
                hitsdict[idnumber] = resultlist
            resultlist = []
        self.hd = hitsdict # keep copy for other functions
 #       print 'hitdict: '
 #       for val in self.hd.values():
 #           for tripl in val:
 #               print 'c1: ',tripl[0]
 #               print 'c2: ',tripl[1]
 #               print 'TP: ',tripl[2]
 #               print ''
        return hitsdict
    
    def append_points(self,ditem,plist):
        '''
        Repeated filling of start and end lists gets its own helper function
        '''
        for v in ditem:
            startlist = v[2].filtered_coordinates1
            for entry in startlist: # construct triplets of coordinates
                plist.append((entry[0],entry[1],v[0].z))
        return plist

                
    def start_end_points(self):
        '''
        return all potential start and end wires (could be more than one each)
        including all the calculated points on them. Order is from 
        left to right, i.e. start to end, simply to have the extreme 
        ends of a track. Doesn't need to be real start and end points, order
        could be swapped upstream.
        Output: 2 lists - list of start coordinates,
                list of end coordinates
        '''
        if len(self.hd)<1: 
            print "Warning: Fill hits dictionary first!"
            return [],[] # empty start and end lists return

        klist = []
        xlist = []
        ylist = []
        for k,val in self.hd.iteritems(): # cycle through all wires
            for tup in val:
                c1 = tup[0] # the cylinder
                klist.append(k)
                xlist.append(c1.x)
                ylist.append(c1.y)

        xarr = abs(xlist) # repeats from number of neighbours
        yarr = array(ylist)
        #print 'Arrays: x and y'
        #print xarr
        #print yarr
        xposmin = where(xarr==min(abs(xlist))) # could be many
        xposmax = where(xarr==max(abs(xlist))) # could be many
        yposmin = where(yarr==min(ylist)) # could be many
        yposmax = where(yarr==max(ylist)) # could be many
        ixmin = xposmin[0][0] # take the first
        iymin = yposmin[0][0] # take the first
        ixmax = xposmax[0][-1] # take the last
        iymax = yposmax[0][-1] # take the last

        #print 'xposmin,xposmax,yposmin,yposmax:'
        #print xposmin,xposmax,yposmin,yposmax

        sxmin = set(xposmin[0])
        sxmax = set(xposmax[0])
        symin = set(yposmin[0])
        symax = set(yposmax[0])

        if (len(sxmin & symin)>0): # pos. slope flag
            pslopeid = sxmin & symin
            #print "pos. slope id"
            #print pslopeid
            ixmin = list(pslopeid)[0] # take the first
            iymin = ixmin
            #pslopeid = sxmax & symax
            #print pslopeid
            #ixmax = list(pslopeid)[0] # take the first
            #iymax = ixmax
            ixmax, iymax = self.posslopeid(xarr,yarr)
            flag = True
        elif (len(sxmin & symax)>0):
            nslopeid = sxmin & symax
            ixmin = list(nslopeid)[0] # take the first
            iymax = ixmin
            #nslopeid = sxmax & symin
            #ixmax = list(nslopeid)[0] # take the first
            #iymin = ixmax
            ixmax, iymin = self.negslopeid(xarr,yarr)
            flag = False
        else:
            flag, ixmin, ixmax, iymin, iymax = self.bentid(xarr,yarr)
            

        #print "Found indices: %d, %d, %d, %d" % (ixmin,ixmax,iymin,iymax)
        #print "Min x, y: (%f, %f)" %(xlist[ixmin],ylist[iymin])
        #print "Max x, y: (%f, %f)" %(xlist[ixmax],ylist[iymax])
        
        start = []
        end = []
        # check for the possible extreme cases
        # vertical
        if (xlist[ixmin]==xlist[ixmax]):
            ikey = klist[iymin]
            startpoint = self.hd[ikey]
            start = self.append_points(startpoint,start)
            ikey = klist[iymax]
            endpoint = self.hd[ikey]
            end = self.append_points(endpoint,end)
            return start, end

        elif (ylist[iymin]==ylist[iymax]): # horizontal
#            print "found horizontal case"
            ikey = klist[ixmin]
            startpoint = self.hd[ikey]
            start = self.append_points(startpoint,start)
            ikey = klist[ixmax]
            endpoint = self.hd[ikey]
            end = self.append_points(endpoint,end)
            return start, end

        if (not flag): # neg. slope, found top left wire
            ikey = klist[ixmin]
            startpoint = self.hd[ikey]
            start = self.append_points(startpoint,start)
            # check bottom right for end
            nbotright = sxmax & symin
            if (len(nbotright)>0): # there is a corner wire
                ixmax = list(nbotright)[0]
                ikey = klist[ixmax]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
            else: # take nearest two to corner wire
                ikey = klist[iymin]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
                ikey = klist[ixmax]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
            return start, end

        else: # pos. slope, found bottom left wire
            #print 'positive slope found'
            ikey = klist[ixmin]
            startpoint = self.hd[ikey]
            start = self.append_points(startpoint,start)
            #print 'starts: '
            #print start
            # check top right for end
            if (ixmax in yposmax[0]):
                ikey = klist[ixmax]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
                #print 'ixmax in yposmax: adding endpoint: ',end
            else: # take nearest two to corner wire
                ikey = klist[iymax]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
                #print 'else: adding endpoint: ',end
                ikey = klist[ixmax]
                endpoint = self.hd[ikey]
                end = self.append_points(endpoint,end)
                #print 'else: adding endpoint: ',end
            return start, end



    def posslopeid(self,xarr,yarr):
        '''
        return indices of start/end points for positive slope structures.
        '''
        xposmax = where(xarr==np.max(xarr)) # could be many
        ymaxarray = yarr[xposmax[0]] # force ymax to be at xmax layer
        yposmax = where(yarr==np.max(ymaxarray)) # should be one
        sxmax = set(xposmax[0])
        symax = set(yposmax[0])
        pslopeid = sxmax & symax # must have an entry
        ixmax = iymax = list(pslopeid)[0]
        return ixmax, iymax

        
    def negslopeid(self,xarr,yarr):
        '''
        return indices of start/end points for negative slope structures.
        '''
        xposmax = where(xarr==np.max(xarr)) # could be many
        yminarray = yarr[xposmax[0]] # force ymin to be at xmax layer
        yposmin = where(yarr==np.min(yminarray)) # should be one
        sxmax = set(xposmax[0])
        symin = set(yposmin[0])
        nslopeid = sxmax & symin # must have an entry
        ixmax = iymin = list(nslopeid)[0]
        return ixmax, iymin

        
    def bentid(self,xarr,yarr):
        '''
        return flag+indices of start/end points for strongly bent structures.
        '''
        xposmin = where(xarr==np.min(xarr)) # could be many
        xposmax = where(xarr==np.max(xarr)) # could be many
        sxmin = set(xposmin[0])
        sxmax = set(xposmax[0])

        # pos slope: ymin at xmin < ymax at xmax
        yminarr = yarr[xposmin[0]] # force ymin to be at xmin layer
        ymaxarr = yarr[xposmax[0]] # force ymax to be at xmax layer
        yposmin = where(yarr==np.min(yminarr)) # should be one
        yposmax = where(yarr==np.max(ymaxarr)) # should be one
        symin = set(yposmin[0])
        symax = set(yposmax[0])
        #print 'pslope bentid sxmin, sxmax',sxmin,sxmax
        #print 'pslope bentid symin, symax',symin,symax
        pslopeid = sxmax & symax # must have an entry
        ixmax = iymax = list(pslopeid)[0]
        pslopeid = sxmin & symin # must have an entry
        ixmin = iymin = list(pslopeid)[0]
        if yarr[iymin]<=yarr[iymax]:
            #print 'min y, max y in pslope: ',yarr[iymin],yarr[iymax]
            flag = True
            return flag,ixmin,ixmax,iymin,iymax
        else:
            flag = False
            dummy = iymin # swap min/max
            iymin = iymax
            iymax = dummy
            #print 'min y, max y in nslope: ',yarr[iymin],yarr[iymax]
            return flag,ixmin,ixmax,iymin,iymax
        
        # neg slope: ymax at xmin > ymin at xmax
#        yposmin = where(yarr==np.max(yminarr)) # should be one
#        yposmax = where(yarr==np.min(ymaxarr)) # should be one
#        symin = set(yposmin[0])
#        symax = set(yposmax[0])
#        print 'nslope bentid symin, symax',symin,symax
#        nslopeid = sxmax & symin # must have an entry
#        ixmax = iymin = list(nslopeid)[0]
#        nslopeid = sxmin & symax # must have an entry
#        ixmin = iymax = list(nslopeid)[0]
#        flag = False
#        return flag,ixmin,ixmax,iymin,iymax


    def condition(self,t1=None,t2=None,err=0.1):
        if ((abs(t1[0]-t2[0])<=self.sigma*err) and (abs(t1[1]-t2[1])<=self.sigma*err)):
            return True
        else:
            return False
        

    def adjust_edges(self,hitsdict=None, err=0.0):
        '''
        allow for radial error by 3 sigma to identify tangent points
        as equal to remove breaks in connections on the rings due to 
        tiny differences in point positions
        '''
        pairlist = []
        for i in hitsdict.keys():
            for val in hitsdict[i]:
                pair = []
                p = val[2].pairs # pick the tangent point pairs
                for entry in p:
                    stpoint=(entry[0][0],entry[0][1],val[0].z)
                    endpoint=(entry[1][0],entry[1][1],val[1].z)
                    pair.append((stpoint,endpoint))

                pairlist.append(pair)

        ll=[]
        for pair in pairlist:
            lltemp = []
            for i,entry in enumerate(pair):
                for j,t1 in enumerate(entry):
                    for m in xrange(i+1,len(pair)):
                        for n,t2 in enumerate(pair[m]):
                            if (self.condition(t1,t2,err)):
                                lltemp.append((i,j,m,n))
            ll.append(lltemp)

        newpairs = []
        for it in xrange(len(pairlist)):
            pair = pairlist[it]
            lltemp = ll[it]
            for entry in lltemp: #entry is a tuple
                (i,j,m,n)=entry
                if (j and n): # position: end, end
                    t=(pair[m][0],pair[i][1])
                elif (not j and not n): # position: start, start
                    t=(pair[i][0],pair[m][1])
                elif (not j and n): # end replaced by start
                    t=(pair[m][0],pair[i][0])
                elif (j and not n): # start replaced by end
                    t=(pair[i][1],pair[m][1])
                pair[m] = t 
            newpairs.append(pair)
        
        return newpairs


    def cleanup(self, datalist):
        tolerance = 0.01
        clean = []
        for entry in datalist:
            for pair in entry:
                for hit in pair:
                    clean.append(hit)

        # sort along 0 coordinate
        clean = sorted(clean,key=lambda t: t[0]) # in-place sort 

        # compare for identicals to clean
        new =  []
        temp = clean[0][0]
        for it in xrange(1,len(clean)):
            if (abs(clean[it][0]-temp)>tolerance): # take hit
                new.append(tuple(clean[it-1]))
                final_flag = True
            else:
                final_flag = False
            temp = clean[it][0]
        if final_flag: # take also the final data hit
            new.append(tuple(clean[len(clean)]-1))

        return new
