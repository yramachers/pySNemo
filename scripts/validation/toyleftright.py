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

Nsims = 1000 # Number of simulated lines

# Set up ROOT data structures for file output and storage
file = root.TFile("/tmp/double_atvertex.tsim","recreate")
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

for i in range(Nsims):
    tgen.double_random_atvertex() # vertex on foil at x=0,y=0
    both = tgen.getLines()
    lines = []

    # enable one line on the left
    lines.append((both[0],0))

    # second line on the right
    lines.append((both[1],1))

    # all hits related truth data in cluster
    cluster = wgr.multi_track_hits(lines)

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

    counter = 0
    for k,val in cluster.iteritems():
        line = lines[k-1][0]  # line3 object
        cells = val[0] # first list in cluster tuple 
        radii = val[1] # as list
        info = val[2]  # as list
        
        dataStruct.dirx.push_back(line.v.x)
        dataStruct.diry.push_back(line.v.y)
        dataStruct.dirz.push_back(line.v.z)
        dataStruct.pointx.push_back(line.p.x)
        dataStruct.pointy.push_back(line.p.y)
        dataStruct.pointz.push_back(line.p.z)

        for w,r,mi in zip(cells,radii,info):
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
    
tree.Write() # write all lines to disk
file.Close()
