import math
import numpy as np
from scipy.ndimage import label
import ROOT as root
import multilines as ML

root.gROOT.ProcessLine(
"struct DataStruct{\
   vector<double>* dirx;\
   vector<double>* diry;\
   vector<double>* dirz;\
   vector<double>* pointx;\
   vector<double>* pointy;\
   vector<double>* pointz;\
   vector<double>* radius;\
   vector<double>* wirex;\
   vector<double>* wirey;\
   vector<double>* wirez;\
   vector<int>*    grid_id;\
   vector<int>*    grid_side;\
   vector<int>*    grid_layer;\
   vector<int>*    grid_column;\
};");


def cleanpicture(c, r, info, yintercept):
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
        return c, r, info
    else:
        ncells = []
        nradii = []
        ninfo  = []
        clusteridx = { }
        collection = []
        for n in range(1,nstruc+1):
            clusteridx[n] = np.where(labels==n) # all indices in labels
        for k in clusteridx.keys():
            for r,c in zip(clusteridx[k][0],clusteridx[k][1]):
                ninfo.append((side,c,r))
                nradii.append(radii[info.index((side,c,r))])
                ncells.append(cells[info.index((side,c,r))])
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


Nsims = 1 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("singlerighthelix.tsim","recreate")
tree = root.TTree("hit_tree","Hit data")
tree.SetDirectory(file)

dataStruct = root.DataStruct()
dataStruct.dirx = root.std.vector('double')()
dataStruct.diry = root.std.vector('double')()
dataStruct.dirz = root.std.vector('double')()
dataStruct.pointx = root.std.vector('double')()
dataStruct.pointy = root.std.vector('double')()
dataStruct.pointz = root.std.vector('double')()
dataStruct.radius = root.std.vector('double')()
dataStruct.wirex = root.std.vector('double')()
dataStruct.wirey = root.std.vector('double')()
dataStruct.wirez = root.std.vector('double')()
dataStruct.gridid = root.std.vector('int')()
dataStruct.gridside = root.std.vector('int')()
dataStruct.gridlayer = root.std.vector('int')()
dataStruct.gridcolumn = root.std.vector('int')()

tree.Branch('dirx', dataStruct.dirx)
tree.Branch('diry', dataStruct.diry)
tree.Branch('dirz', dataStruct.dirz)
tree.Branch('pointx', dataStruct.pointx)
tree.Branch('pointy', dataStruct.pointy)
tree.Branch('pointz', dataStruct.pointz)
tree.Branch('radius', dataStruct.radius)
tree.Branch('wirex',  dataStruct.wirex)
tree.Branch('wirey',  dataStruct.wirey)
tree.Branch('wirez',  dataStruct.wirez)
tree.Branch('grid_id', dataStruct.gridid)
tree.Branch('grid_side', dataStruct.gridside)
tree.Branch('grid_layer', dataStruct.gridlayer)
tree.Branch('grid_column', dataStruct.gridcolumn)

wgr = ML.demonstratorgrid()
tgen = ML.helix_generator() # only 2D helices!
yinterc = 0.0
struc = tgen.single_random_momentum(yinterc) # x=0 fixed for this generator

cells, radii = wgr.hits(struc, 1) # right tracker half
info = wgr.wireinfo


file.cd()
# Prepare data structure for this line
dataStruct.dirx.clear()
dataStruct.diry.clear()
dataStruct.dirz.clear()
dataStruct.pointx.clear()
dataStruct.pointy.clear()
dataStruct.pointz.clear()
dataStruct.radius.clear()
dataStruct.wirex.clear()
dataStruct.wirey.clear()
dataStruct.wirez.clear()
dataStruct.gridid.clear()
dataStruct.gridside.clear()
dataStruct.gridlayer.clear() 
dataStruct.gridcolumn.clear()

# save the geometry truth data
# save a helix as parameter set
refp = struc.referencePoint # triplet
mom  = struc.momentum       # triplet, z=0 by default
charge = struc.charge       # number
dataStruct.dirx.push_back(mom[0]) # momenta in here
dataStruct.diry.push_back(mom[1])
dataStruct.dirz.push_back(charge) # use dirz for charge
dataStruct.pointx.push_back(refp[0]) # reference point here
dataStruct.pointy.push_back(refp[1])
dataStruct.pointz.push_back(refp[2])

print 'for helix: '
print struc
print info
if len(info): # only if there is any data at all
    ncells, nradii, ninfo = cleanpicture(cells, radii, info, yinterc) # remove returning helix branch

print 'entries: %d'%len(ninfo)
for w,r,mi in zip(ncells, nradii, ninfo):
    dataStruct.radius.push_back(r)
    dataStruct.wirex.push_back(w[0])
    dataStruct.wirey.push_back(w[1])
    dataStruct.wirez.push_back(0.0)
    dataStruct.gridid.push_back(0)
    side = mi[0] # wire side
    row = mi[1] # wire column
    col = mi[2] # wire layer
    dataStruct.gridlayer.push_back(row)
    dataStruct.gridcolumn.push_back(col)
    dataStruct.gridside.push_back(side)
    print mi

tree.Fill() # data structure fully filled, lines done

tree.Write() # write all lines to disk
file.Close()
