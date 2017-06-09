import math
import random
import numpy as np
import ROOT as root
import multilines as ML
from pysnemo.utility import geometrycheck as gcheck
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

def remove_hits(fromlayer, tolayer, cells, radii, info):
    infoarr = np.array(info)
    cellarr = np.array(cells)
    radarr  = np.array(radii)
    for lay in range(fromlayer,tolayer):
        indx = np.where(infoarr[:,1]!=lay)[0]
        infoarr = infoarr[indx]
        cellarr = cellarr[indx]
        radarr  = radarr[indx]
    return cellarr.tolist(), radarr.tolist(), infoarr.tolist()


def remove_doubles(cluster):
    k=1 # keys start at 1
    if len(cluster[k][2])>0 and len(cluster[k+1][2])>0: # one empty = no doubles
        for pos, mi in enumerate(cluster[k][2]): # check info entries against all of next cluster entry
            nextinfo = np.array(cluster[k+1][2])  # from modified next cluster
            if len(nextinfo)<1:
                break # nothing to compare to
            dublet = list(mi[-2:]) # comes as tuple -> list convert
            nilist = nextinfo[:,-2:].tolist() # smart slicing, back to list
            
            if dublet in nilist: # found a double entry
                indx = nilist.index(dublet) # first entry; should be only one
                nextradius = cluster[k+1][1][indx] # choose according to radius
                if cluster[k][1][pos] > nextradius: # smaller wins
                    cluster[k][1][pos] = nextradius # in-place change
                cluster[k+1][0].pop(indx) # change cluster k+1
                cluster[k+1][1].pop(indx)
                cluster[k+1][2].pop(indx)
    return cluster


Nsims = 1000 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("/tmp/multiscatter_calo.tsim","recreate")
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
tgen = ML.track_generator()
dcalo = gcheck.demonstratorcalo()

for i in range(Nsims):
    cluster = { }
    lines = []
    bpoints = []
    bangles = []
    blayer = []
    scatter_angle = 3.0 # multiple scattering angle width [degrees]
    lrtracker = random.randint(0,1) # random left or right tracker side
    sign = -1*(-1)**lrtracker

    # random line xy slope for the first straight line
    angle = random.uniform(-math.pi*0.5+0.17, math.pi*0.5-0.17) # taking vertical out
    sl = math.tan(angle)

    # make a first line with cells etc.
    dummy = tgen.single_line_manual_with_z(sl,0.0,0.0,5.0) # Line3 with vertex on foil at x=0,y=0,z=5.0
    cells, radii = wgr.hits(dummy, lrtracker) # left/right tracker half
    info = wgr.wireinfo
    if len(info)>0:
        # 2D projection not an issue since is parallel to z=0 plane by construction
        original = euclid.Line2(euclid.Point2(dummy.p.x, dummy.p.y), euclid.Vector2(dummy.v.x, dummy.v.y))
        lines.append(dummy)

        # break first line and pick scatter angle for continuation
        bl = random.randint(1,7) # tracker layer with break
        cellvariation = random.uniform(-21.0, 21.0)
        layerline = euclid.Line2(euclid.Point2(sign*53.0 + sign*bl*44.0 + cellvariation,0.0), euclid.Vector2(0.0,1.0))
        breakpoint = original.intersect(layerline)
        bpoints.append((breakpoint.x, breakpoint.y))
        blayer.append(bl)
        c, r, i = remove_hits(bl+1, 9, cells, radii, info) # return numpy arrays
        cluster[1] = (c, r, i)

        # next line continuing from breakpoint
        scat_angle = random.gauss(0.0,scatter_angle*math.pi/180.0) # random scattering angle
        newangle = angle + scat_angle # altering original slope angle with new angle
        sl = math.tan(newangle)
        nextdummy = euclid.Line3(euclid.Point3(breakpoint.x, breakpoint.y, 5.0), euclid.Vector3(1.0, sl, 0.0))
        caloinfo= dcalo.calohits(nextdummy, lrtracker)
        while len(caloinfo) < 1: # no calo was hit, try again
            scat_angle = random.gauss(0.0,scatter_angle*math.pi/180.0) # random scattering angle
            newangle = angle + scat_angle # altering original slope angle with new angle
            sl = math.tan(newangle)
            nextdummy = euclid.Line3(euclid.Point3(breakpoint.x, breakpoint.y, 5.0), euclid.Vector3(1.0, sl, 0.0))
            caloinfo= dcalo.calohits(nextdummy, lrtracker)

        bangles.append(scat_angle)
        lines.append(nextdummy)
        ncells, nradii = wgr.hits(nextdummy, lrtracker) # left/right tracker half
        ninfo = wgr.wireinfo

        if len(ninfo)>0:
            c2, r2, i2 = remove_hits(0, bl, ncells, nradii, ninfo) # return numpy arrays
            cluster[2] = (c2, r2, i2)
            cluster = remove_doubles(cluster)

            allcells = cluster[1][0] + cluster[2][0] # concatenate
            allradii = cluster[1][1] + cluster[2][1] # concatenate
            allinfo  = cluster[1][2] + cluster[2][2] # concatenate

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

            for entry in lines: # truth lines
                dataStruct.dirx.push_back(entry.v.x)
                dataStruct.diry.push_back(entry.v.y)
                dataStruct.dirz.push_back(entry.v.z)
                dataStruct.pointx.push_back(entry.p.x)
                dataStruct.pointy.push_back(entry.p.y)
                dataStruct.pointz.push_back(entry.p.z)
                dataStruct.charge.push_back(0)

            for bp, bang, bl in zip(bpoints, bangles, blayer):
                dataStruct.bpointx.push_back(bp[0])
                dataStruct.bpointy.push_back(bp[1])
                dataStruct.bangle.push_back(bang)
                dataStruct.breaklayer.push_back(bl)

            counter = 0
            for w,r,mi in zip(allcells,allradii,allinfo):
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

            type = caloinfo[0][1]
            side = caloinfo[0][3]
            col  = caloinfo[0][4]
            row  = caloinfo[0][5]
            wall = caloinfo[0][6]
            dataStruct.caloid.push_back(0)
            dataStruct.calorow.push_back(row)
            dataStruct.calocolumn.push_back(col)
            dataStruct.calotype.push_back(type)
            dataStruct.caloside.push_back(side)
            dataStruct.calowall.push_back(wall)
        

            tree.Fill() # data structure fully filled, lines done
    
tree.Write() # write all lines to disk
file.Close()
