import math
import numpy as np
import ROOT as root
import multilines as ML
from scipy.ndimage import label
from pysnemo.utility import geometrycheck as gcheck

root.gROOT.ProcessLine(
"struct DataStruct{\
   vector<double>* dirx;\
   vector<double>* diry;\
   vector<double>* dirz;\
   vector<double>* pointx;\
   vector<double>* pointy;\
   vector<double>* pointz;\
   vector<double>* breakpointx;\
   vector<double>* breakpointy;\
   vector<double>* bpangle;\
   vector<double>* radius;\
   vector<double>* wirex;\
   vector<double>* wirey;\
   vector<double>* wirez;\
   vector<int>*    grid_id;\
   vector<int>*    grid_side;\
   vector<int>*    grid_layer;\
   vector<int>*    grid_column;\
   vector<int>*    break_layer;\
   vector<int>*    charge;\
   vector<int>*    calo_id;\
   vector<int>*    calo_type;\
   vector<int>*    calo_side;\
   vector<int>*    calo_wall;\
   vector<int>*    calo_column;\
   vector<int>*    calo_row;\
};");

def cleanpicture(clist, rlist, info, yintercept):
    '''
    search for connected structures, 
    select the one close to the reference point of the helix which
    removes the returning branch if it exists.
    '''
    pic = np.zeros((113,9)) # half tracker hardwired
    side = info[0][0] # same for all
    rowpos = int(math.floor((yintercept+2464.0) / 44.0))

    # turn wireinfo into picture
    for entry in info:
        pic[entry[2],entry[1]] = 1
    anyslope = [[1,1,1], [1,1,1], [1,1,1]] # any slope structure

    labels, nstruc = label(pic,anyslope)

    if nstruc <2: # nothing else to do
        return clist, rlist, info
    else:
        ncells = []
        nradii = []
        ninfo  = []
        clusteridx = { }
        collection = []
        for n in range(1,nstruc+1):
            clusteridx[n] = np.where(labels==n) # all indices in labels
        for k in clusteridx.keys():
            for row,column in zip(clusteridx[k][0],clusteridx[k][1]):
                ninfo.append((side,column,row))
                nradii.append(rlist[info.index((side,column,row))])
                ncells.append(clist[info.index((side,column,row))])
            collection.append((ncells, nradii, ninfo)) # structures fully separated in collection
            ninfo  = []
            ncells = []
            nradii = []
        pick = 0
        for n, entry in enumerate(collection):
            if (side, 0, rowpos) in entry[2]:
                pick = n
            elif (side, 0, rowpos-1) in entry[2]: # allow for nearest positions
                pick = n
            elif (side, 0, rowpos+1) in entry[2]: # +-1 around the integer row
                pick = n
        return collection[pick]


def main_wall_test(caloinfo, clusterentry):
    info = clusterentry[2]
    cells = clusterentry[0]
    # checks piercings of the main wall in case there are two
    rowpos = caloinfo[5]*259.0 - 1683.5 # pos in z
    
    for winfo, entry in zip(info, cells):
        if winfo[1]==8: # for final wire row only, check position in z
            if entry[2] >= rowpos-259.0 and entry[2] <= rowpos+259.0: # 50cm calo z interval
                return True
            else:
                return False
            

Nsims = 1000 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("/tmp/2helix_calo.tsim","recreate")
tree = root.TTree("hit_tree","Hit data")
tree.SetDirectory(file)

dataStruct = root.DataStruct()
dataStruct.dirx = root.std.vector('double')()
dataStruct.diry = root.std.vector('double')()
dataStruct.dirz = root.std.vector('double')()
dataStruct.pointx = root.std.vector('double')()
dataStruct.pointy = root.std.vector('double')()
dataStruct.pointz = root.std.vector('double')()
dataStruct.bpointx = root.std.vector('double')()
dataStruct.bpointy = root.std.vector('double')()
dataStruct.bangle = root.std.vector('double')()
dataStruct.radius = root.std.vector('double')()
dataStruct.wirex = root.std.vector('double')()
dataStruct.wirey = root.std.vector('double')()
dataStruct.wirez = root.std.vector('double')()
dataStruct.gridid = root.std.vector('int')()
dataStruct.gridside = root.std.vector('int')()
dataStruct.gridlayer = root.std.vector('int')()
dataStruct.gridcolumn = root.std.vector('int')()
dataStruct.breaklayer = root.std.vector('int')()
dataStruct.charge = root.std.vector('int')()
dataStruct.caloid = root.std.vector('int')()
dataStruct.calotype = root.std.vector('int')()
dataStruct.calowall = root.std.vector('int')()
dataStruct.caloside = root.std.vector('int')()
dataStruct.calorow = root.std.vector('int')()
dataStruct.calocolumn = root.std.vector('int')()

tree.Branch('dirx', dataStruct.dirx)
tree.Branch('diry', dataStruct.diry)
tree.Branch('dirz', dataStruct.dirz)
tree.Branch('pointx', dataStruct.pointx)
tree.Branch('pointy', dataStruct.pointy)
tree.Branch('pointz', dataStruct.pointz)
tree.Branch('breakpointx', dataStruct.bpointx)
tree.Branch('breakpointy', dataStruct.bpointy)
tree.Branch('bpangle', dataStruct.bangle)
tree.Branch('radius', dataStruct.radius)
tree.Branch('wirex',  dataStruct.wirex)
tree.Branch('wirey',  dataStruct.wirey)
tree.Branch('wirez',  dataStruct.wirez)
tree.Branch('grid_id', dataStruct.gridid)
tree.Branch('grid_side', dataStruct.gridside)
tree.Branch('grid_layer', dataStruct.gridlayer)
tree.Branch('grid_column', dataStruct.gridcolumn)
tree.Branch('break_layer', dataStruct.breaklayer)
tree.Branch('charge', dataStruct.charge)
tree.Branch('calo_id', dataStruct.caloid)
tree.Branch('calo_type', dataStruct.calotype)
tree.Branch('calo_side', dataStruct.caloside)
tree.Branch('calo_wall', dataStruct.calowall)
tree.Branch('calo_row', dataStruct.calorow)
tree.Branch('calo_column', dataStruct.calocolumn)

wgr = ML.demonstratorgrid()
tgen = ML.helix_generator()
dcalo = gcheck.demonstratorcalo()
yinterc = 0.0
nonhit = 0
for i in range(Nsims):
    struc = tgen.double_helix() # x=0 fixed for this generator
    both = tgen.getHelices()
    lines = []

    # enable one line on the left
    lines.append((both[0],0))

    # second line on the right
    lines.append((both[1],1))

    # all hits related truth data in cluster
    cluster = wgr.multi_track_hits(lines)
    cluster2= dcalo.multi_calohits(lines)

    while len(cluster2) < 2: # no calo was hit, try again
        struc = tgen.double_helix() # x=0 fixed for this generator
        both = tgen.getHelices()
        lines = []
        lines.append((both[0],0))
        lines.append((both[1],1))
        cluster = wgr.multi_track_hits(lines)
        cluster2= dcalo.multi_calohits(lines)

    # clean tracker traces
    for k,item in cluster.iteritems():
        c, r, info = item
        # print 'cells before cleaning, key,cells',k,c
        ncells, nradii, ninfo = cleanpicture(c, r, info, yinterc) # remove returning helix branch
        # print 'cells after cleaning',ncells
        # overwrite
        cluster[k] = (ncells, nradii, ninfo)

        
    file.cd()
    # Prepare data structure for this line
    dataStruct.dirx.clear()
    dataStruct.diry.clear()
    dataStruct.dirz.clear()
    dataStruct.pointx.clear()
    dataStruct.pointy.clear()
    dataStruct.pointz.clear()
    dataStruct.bpointx.clear()
    dataStruct.bpointy.clear()
    dataStruct.bangle.clear()
    dataStruct.radius.clear()
    dataStruct.wirex.clear()
    dataStruct.wirey.clear()
    dataStruct.wirez.clear()
    dataStruct.gridid.clear()
    dataStruct.gridside.clear()
    dataStruct.gridlayer.clear() 
    dataStruct.gridcolumn.clear()
    dataStruct.breaklayer.clear()
    dataStruct.charge.clear()
    dataStruct.caloid.clear()
    dataStruct.calotype.clear()
    dataStruct.caloside.clear()
    dataStruct.calowall.clear()
    dataStruct.calorow.clear() 
    dataStruct.calocolumn.clear()

    # clean multiple calo hits
    write = True
    for k, caloinfo in cluster2.iteritems():
        (cilist,pointlist) = caloinfo
        ci = cilist[0]
        point = pointlist[0]
        if ci[1]==0:  # main wall
            if ci[3]==k-1: # correct side
                if main_wall_test(ci, cluster[k]):
                    type = ci[1]
                    side = ci[3]
                    col  = ci[4]
                    row  = ci[5]
                    wall = ci[6]
                #print 'picked: ',ci
                    calo_hit_point = point # back as tuple here
                #print 'at impact point: ',calo_hit_point
                else:
                    write = False
                    continue # loose that helix
            else:
                write = False
                continue # loose that helix
        else:
            write = False
            continue # loose that helix
        if write:
            dataStruct.pointx.push_back(calo_hit_point.x) # hit point on calo
            dataStruct.pointy.push_back(calo_hit_point.y)
            dataStruct.pointz.push_back(calo_hit_point.z)
            dataStruct.caloid.push_back(k-1)
            dataStruct.calorow.push_back(row)
            dataStruct.calocolumn.push_back(col)
            dataStruct.calotype.push_back(type)
            dataStruct.caloside.push_back(side)
            dataStruct.calowall.push_back(wall)
    counter = 0
    if write:
        for k,val in cluster.iteritems():
            struc = both[k-1]  # helix object
            cells = val[0] # first list in cluster tuple 
            radii = val[1] # as list
            info = val[2]  # as list
            
            mom  = struc.momentum       # triplet, z=0 by default
            charge = int(struc.charge)       # number
            dataStruct.dirx.push_back(mom[0]) # momenta in here
            dataStruct.diry.push_back(mom[1])
            dataStruct.dirz.push_back(mom[2])
            dataStruct.charge.push_back(charge)
            
            for w,r,mi in zip(cells,radii,info):
                dataStruct.radius.push_back(r)
                dataStruct.wirex.push_back(w[0])
                dataStruct.wirey.push_back(w[1])
                dataStruct.wirez.push_back(w[2])
                dataStruct.gridid.push_back(counter)
                side = mi[0] # wire side
                row = mi[1] # wire column
                col = mi[2] # wire layer
                dataStruct.gridlayer.push_back(row)
                dataStruct.gridcolumn.push_back(col)
                dataStruct.gridside.push_back(side) # not covered yet 
                counter += 1 # count up all hits for entire event


        tree.Fill() # data structure fully filled, lines done
    else:
        nonhit += 1
print 'made %d events with %d non hit events'%(Nsims-nonhit,nonhit)
tree.Write() # write all lines to disk
file.Close()
