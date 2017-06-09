import math
import random
import ROOT as root
import multilines as ML
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

Nsims = 1000 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("/tmp/illumination.tsim","recreate")
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
    # random line xy slope for the straight line
    angle = random.uniform(-math.pi*0.5+0.1, math.pi*0.5-0.1) # taking vertical out
    sl = math.tan(angle)
    # random line xz slope for the straight line
    angle2 = random.uniform(-math.pi*0.5+0.1, math.pi*0.5-0.1) # taking vertical out
    slxz = math.tan(angle2)

    dummy = tgen.single_line_manual_with_z(sl, slxz, 0.0, 0.0) # vertex on foil at x=0,y=0,z=0
    lrtracker = random.randint(0,1) # pick left or right

    cells, radii = wgr.hits(dummy, lrtracker) # left/right tracker half
    info = wgr.wireinfo
    caloinfo= dcalo.calohits(dummy, lrtracker)
    while len(caloinfo) < 1: # no calo was hit, try again
        angle = random.uniform(-math.pi*0.5+0.1, math.pi*0.5-0.1) # taking vertical out
        sl = math.tan(angle)
        angle2 = random.uniform(-math.pi*0.5+0.1, math.pi*0.5-0.1) # taking vertical out
        slxz = math.tan(angle2)
        dummy = tgen.single_line_manual_with_z(sl, slxz, 0.0, 0.0) # vertex on foil at x=0,y=0
        lrtracker = random.randint(0,1) # pick left or right
        
        cells, radii = wgr.hits(dummy, lrtracker) # left/right tracker half
        info = wgr.wireinfo
        caloinfo= dcalo.calohits(dummy, lrtracker)
    calo_hit_point = dcalo.get_point()

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

    dataStruct.dirx.push_back(dummy.v.x)
    dataStruct.diry.push_back(dummy.v.y)
    dataStruct.dirz.push_back(dummy.v.z)
    dataStruct.pointx.push_back(calo_hit_point.x)
    dataStruct.pointy.push_back(calo_hit_point.y)
    dataStruct.pointz.push_back(calo_hit_point.z)
    dataStruct.charge.push_back(0)
#    dataStruct.pointx.push_back(dummy.p.x) # origin point
#    dataStruct.pointy.push_back(dummy.p.y)
#    dataStruct.pointz.push_back(dummy.p.z)

    counter = 0
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
