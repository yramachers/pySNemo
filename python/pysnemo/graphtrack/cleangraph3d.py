import networkx as nx
import numpy as np
from math import sqrt

class cleangraph(object):
    '''
    Build a graph object from a hit dictionary and clean up the nodes
    and edges by allowing for a radial error, i.e. identify points
    too close according to dr as identical.
    Provides function 'find_path':
    Input: start and end node on graph
    Output: shortest path as list of nodes for fitting and
            path length as euclidean distance between start and end
    '''
    def __init__(self,hitsdict=None, sigma=1.0):
        self.sigma = sigma
        self.graph = nx.Graph() # un-directed empty graph to be filled
        self.fill_graph(hitsdict) # does the filling
        self.hd = hitsdict # local copy for access by functions


    def rphiToxy(self,x,y,r,dr,dphi):
        if (abs(x)>0.0):
            dx = abs(x/(2*r)*(dr-2*r*dphi*y/x))
            dy = abs(r*r/x*(dphi+y*dx/(r*r)))
        else:  # singularity for horizontal lines
            dx=sqrt(0.5)*dr # set absolute minimal error in [cm]
            dy=sqrt(0.5)*dr

        # don't permit too small uncertainties
        if (sqrt(dx**2+dy**2)<=dr):
            dx = sqrt(0.5)*dr # set artificially symmetric
            dy = sqrt(0.5)*dr # set artificially symmetric
        return dx,dy


    def fill_graph(self,hitsdict):
        '''
        First builds the graph simply from all points. All duplicates are
        already removed thanks to the graph package. Remains to map 
        close misses onto each other and clean the nodes and edges this way.
        '''
        tgraph = nx.Graph() # un-directed empty graph to be filled
        for i,val in hitsdict.iteritems():
            for tup in val:
                pair = []
                # tup[0] start ggcell, tup[1] end ggcell, tup[2] tpcontroller
                p = tup[2].pairs 
                perr = tup[2].pair_errors
                # then create points with coords and errors together
                for coords, error in zip(p,perr):
                    x=coords[0][0]-tup[0].xwire
                    y=coords[0][1]-tup[0].ywire
                    r=tup[0].radius
                    dr=error[0][0]
                    dphi=error[0][1]
                    dx,dy = self.rphiToxy(x,y,r,dr,dphi)
                    stpoint=(coords[0][0],coords[0][1],tup[0].zwire,dx,dy,tup[0].dz)

                    x=coords[1][0]-tup[1].xwire
                    y=coords[1][1]-tup[1].ywire
                    r=tup[1].radius
                    dr=error[1][0]
                    dphi=error[1][1]
                    dx,dy = self.rphiToxy(x,y,r,dr,dphi)
                    endpoint=(coords[1][0],coords[1][1],tup[1].zwire,dx,dy,tup[1].dz)
                    pair.append((stpoint,endpoint))
                #print pair
                tgraph.add_edges_from(pair)

        ed = tgraph.edges()
        #print ed
        edgelist = self.adjust_edges(ed) # clean the edge list
        #print edgelist
        self.graph.add_weighted_edges_from(edgelist) # Nodes are built automatically

    def condition(self,t1,t2):
        errx = sqrt(t1[3]**2 + t2[3]**2)
        erry = sqrt(t1[4]**2 + t2[4]**2)
        if ((abs(t1[0]-t2[0])<=self.sigma*errx) and (abs(t1[1]-t2[1])<=self.sigma*erry) and (abs(t1[0]-t2[0])>0.0)):
            return True
        else:
            return False
        
    def average(self,values):
        varr = np.array(values)
        coords = list(np.average(varr,axis=0)[:3])
        err = list(np.sqrt(np.sum(np.square(varr),axis=0))[3:])
        return tuple(coords+err)


    def dist(self,a,b):
        (x1,y1,z1) = a[:3]
        (x2,y2,z2) = b[:3]
        return sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        

    def adjust_edges(self,edgelist=None):
        '''
        allow for radial error by 3 sigma to identify tangent points
        as equal to remove breaks in connections on the rings due to 
        tiny differences in point positions
        '''
        slist=[]
        elist=[]
        for entry in edgelist:
            slist.append(entry[0])
            elist.append(entry[1])

        sindlist = []
        eindlist = []
        mixlist = []
        for n,(s,e) in enumerate(zip(slist,elist)):
            # start on start
            for m,s2 in enumerate(slist):
                if self.condition(s,s2):
                    sindlist.append((n,m))
            # end on end
            for m,e2 in enumerate(elist):
                if self.condition(e,e2):
                    eindlist.append((n,m))
            # start on end
            for m,e2 in enumerate(elist):
                if self.condition(s,e2):
                    mixlist.append((n,m))
        #print 'slist: ',sindlist
        #print 'elist: ',eindlist
        #print 'mix: ',mixlist
        gg = nx.Graph()
        for entry in sindlist:
            gg.add_edge('s.'+str(entry[0]),'s.'+str(entry[1]))

        for entry in eindlist:
            gg.add_edge('e.'+str(entry[0]),'e.'+str(entry[1]))

        for entry in mixlist:
            gg.add_edge('s.'+str(entry[0]),'e.'+str(entry[1]))

        overlap = nx.connected_components(gg)
        for entry in overlap: # list of lists
            signature=[]
            ind = []
            values = []
            for node in entry:
                s,v = node.split('.')
                signature.append(s)
                ind.append(int(v))
                if s == 's':
                    values.append(slist[int(v)])
                elif s == 'e':
                    values.append(elist[int(v)])
            newval = self.average(values)
            # distribute new average for all
            for sig,k in zip(signature,ind):
                if sig == 's':
                    slist[k] = newval
                elif sig == 'e':
                    elist[k] = newval
                    
        # back to edges
        newedgelist = []
        for s,e in zip(slist,elist):
            newedgelist.append((s,e))
        # give the edges a euclidean distance weight
        ed = []
        for val in newedgelist:
            s = val[0]
            e = val[1]
            w = self.dist(s,e)
            ed.append((s,e,w))
        return ed


    def find_path(self,start=(0.0,0.0,0.0),end=(0.0,0.0,0.0)):
        '''
        find the shortest path between start and end nodes on the graph
        using the euclidean metric instead of graph metric
        Output: list of nodes of shortest path and absolute distance
        '''
        if (start==end): # prevent lack of sensical nodes
            return [], 0.0

        if (len(nx.predecessor(self.graph,start,end)) > 0):# if reachable
            #path = nx.shortest_path(self.graph,start,end,'weight')
            path = nx.dijkstra_path(self.graph,start,end,'weight')
            pl = nx.shortest_path_length(self.graph,start,end,'weight')
            return path, pl
        else: # broken graph between start and end; no shortest path
            return [], 0.0
            

