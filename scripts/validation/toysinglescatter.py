import math
import random
import numpy as np
import ROOT as root
import multilines as ML
from pysnemo.utility import euclid

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
   vector<double>* radius;\
   vector<double>* wirex;\
   vector<double>* wirey;\
   vector<double>* wirez;\
   vector<int>*    grid_id;\
   vector<int>*    grid_side;\
   vector<int>*    grid_layer;\
   vector<int>*    grid_column;\
   vector<int>*    break_layer;\
};");


def remove_hits(fromlayer, tolayer, cells, radii, info):
    infoarr = np.array(info)
    cellarr = np.array(cells)
    radarr  = np.array(radii)
    for lay in range(fromlayer,tolayer+1):
        indx = np.where(infoarr[:,1]!=lay)[0]
        infoarr = infoarr[indx]
        cellarr = cellarr[indx]
        radarr  = radarr[indx]
    return cellarr, radarr, infoarr


Nsims = 1000 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("/tmp/singlebreakp.tsim","recreate")
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
dataStruct.radius = root.std.vector('double')()
dataStruct.wirex = root.std.vector('double')()
dataStruct.wirey = root.std.vector('double')()
dataStruct.wirez = root.std.vector('double')()
dataStruct.gridid = root.std.vector('int')()
dataStruct.gridside = root.std.vector('int')()
dataStruct.gridlayer = root.std.vector('int')()
dataStruct.gridcolumn = root.std.vector('int')()
dataStruct.breaklayer = root.std.vector('int')()

tree.Branch('dirx', dataStruct.dirx)
tree.Branch('diry', dataStruct.diry)
tree.Branch('dirz', dataStruct.dirz)
tree.Branch('pointx', dataStruct.pointx)
tree.Branch('pointy', dataStruct.pointy)
tree.Branch('pointz', dataStruct.pointz)
tree.Branch('breakpointx', dataStruct.bpointx)
tree.Branch('breakpointy', dataStruct.bpointy)
tree.Branch('radius', dataStruct.radius)
tree.Branch('wirex',  dataStruct.wirex)
tree.Branch('wirey',  dataStruct.wirey)
tree.Branch('wirez',  dataStruct.wirez)
tree.Branch('grid_id', dataStruct.gridid)
tree.Branch('grid_side', dataStruct.gridside)
tree.Branch('grid_layer', dataStruct.gridlayer)
tree.Branch('grid_column', dataStruct.gridcolumn)
tree.Branch('break_layer', dataStruct.breaklayer)

wgr = ML.demonstratorgrid()
tgen = ML.track_generator()

for nsims in range(Nsims):
    lines = []
    bpoints = []
    blayer = []
    scatter_angle = 3.0 # multiple scattering angle width [degrees]
    lrtracker = random.randint(0,1) # random left or right tracker side

    # random line slope for the first straight line
    angle = math.pi/2.0
    while (angle>(0.5*math.pi-0.08) or angle<(-0.5*math.pi+0.08)):# cut vertical out
        angle = 0.5*random.vonmisesvariate(0.0,0) #uniform angle (0,pi)
    sl = math.tan(angle)

    # make a first line with cells etc.
    dummy = tgen.single_line_manual(sl,0.0) # Line3 with vertex on foil at x=0,y=0
    cells, radii = wgr.hits(dummy, lrtracker) # left/right tracker half
    info = wgr.wireinfo
    # 2D projection not an issue since is in z=0 plane by construction
    original = euclid.Line2(euclid.Point2(dummy.p.x, dummy.p.y), euclid.Vector2(dummy.v.x, dummy.v.y))
    lines.append(dummy)

    # break first line and pick scatter angle for continuation
    bl = random.randint(1,7) # tracker layer with break
    cellvariation = random.uniform(-21.0, 21.0)
    layerline = euclid.Line2(euclid.Point2(53.0 + bl*44.0 + cellvariation,0.0), euclid.Vector2(0.0,1.0))
    breakpoint = original.intersect(layerline)
    c, r, i = remove_hits(bl, 8, cells, radii, info) # return numpy arrays
    bpoints.append((breakpoint.x, breakpoint.y))
    blayer.append(bl)
    
    scat_angle = random.gauss(0.0,scatter_angle*math.pi/180.0) # random scattering angle

    # next line continuing from breakpoint
    angle += scat_angle # altering original slope angle with new angle
    sl = math.tan(angle)
    nextdummy = euclid.Line3(euclid.Point3(breakpoint.x, breakpoint.y,0.0), euclid.Vector3(1.0, sl, 0.0))
    lines.append(nextdummy)
    ncells, nradii = wgr.hits(nextdummy, lrtracker) # left/right tracker half
    ninfo = wgr.wireinfo
    if len(ninfo)>1:
        c2, r2, i2 = remove_hits(0, bl-1, ncells, nradii, ninfo) # return numpy arrays

        allcells = np.concatenate((c, c2)) # concatenate
        allradii = np.concatenate((r, r2))
        allinfo  = np.concatenate((i, i2))
        
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
        dataStruct.radius.clear()
        dataStruct.wirex.clear()
        dataStruct.wirey.clear()
        dataStruct.wirez.clear()
        dataStruct.gridid.clear()
        dataStruct.gridside.clear()
        dataStruct.gridlayer.clear() 
        dataStruct.gridcolumn.clear()
        dataStruct.breaklayer.clear()
        
        for entry in lines: # truth lines
            dataStruct.dirx.push_back(entry.v.x)
            dataStruct.diry.push_back(entry.v.y)
            dataStruct.dirz.push_back(entry.v.z)
            dataStruct.pointx.push_back(entry.p.x)
            dataStruct.pointy.push_back(entry.p.y)
            dataStruct.pointz.push_back(entry.p.z)
            
        for bp, bl in zip(bpoints, blayer):
            dataStruct.bpointx.push_back(bp[0])
            dataStruct.bpointy.push_back(bp[1])
            dataStruct.breaklayer.push_back(bl)

    
        counter = 0
        for w,r,mi in zip(allcells.tolist(),allradii.tolist(),allinfo.tolist()):
            dataStruct.radius.push_back(r)
            dataStruct.wirex.push_back(w[0])
            dataStruct.wirey.push_back(w[1])
            dataStruct.wirez.push_back(0.0)
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
        continue
    
tree.Write() # write all lines to disk
file.Close()
