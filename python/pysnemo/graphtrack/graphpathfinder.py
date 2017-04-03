import networkx as nx
import numpy as np
from math import sqrt, atan2, sin, cos

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
        phi = atan2(y,x)
        dx = abs(-r*sin(phi)*dphi+dr*cos(phi))
        if dx<dr:
            dx=dr # minimum absolute error cut
        dy = abs(r*cos(phi)*dphi+dr*sin(phi))
        if dy<dr:
            dy=dr
        return round(dx,4), round(dy,4)


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
                idl = tup[2].idlist
                #print 'idlist in tpcontroller: ',idl
                # then create points with coords and errors together
                for n,(coords, error) in enumerate(zip(p,perr)):
                    id0=idl[n]
                    x=coords[0][0]-tup[0].x
                    y=coords[0][1]-tup[0].y
                    r=tup[0].r
                    dr=error[0][0]
                    dphi=error[0][1]
                    dx,dy = self.rphiToxy(x,y,r,dr,dphi)
                    stpoint=(round(coords[0][0],4),round(coords[0][1],4),round(tup[0].z,4),dx,dy,round(tup[0].sigmaz,4),id0,tup[0].meta_info)

                    id1=idl[n+1]
                    x=coords[1][0]-tup[1].x
                    y=coords[1][1]-tup[1].y
                    r=tup[1].r
                    dr=error[1][0]
                    dphi=error[1][1]
                    dx,dy = self.rphiToxy(x,y,r,dr,dphi)
                    endpoint=(round(coords[1][0],4),round(coords[1][1],4),round(tup[1].z,4),dx,dy,round(tup[1].sigmaz,4),id1,tup[1].meta_info)
                    pair.append((stpoint,endpoint))
                #print pair
                tgraph.add_edges_from(pair)

        ed = tgraph.edges()
        edgelist = self.adjust_edges(ed) # clean the edge list
        #print edgelist
        self.graph.add_weighted_edges_from(edgelist) # Nodes are built automatically

    def condition(self,t1,t2):
        errx = sqrt(t1[3]**2 + t2[3]**2)
        erry = sqrt(t1[4]**2 + t2[4]**2)
        if ((abs(t1[0]-t2[0])<=self.sigma*errx) and (abs(t1[1]-t2[1])<=self.sigma*erry)):
            return True
        else:
            return False
        
    def average(self,values):
        val = [(v[:6]) for v in values]
        ids = [v[6] for v in values]
        mi = [v[7] for v in values]
        idset = set(ids)
        mis = set(mi) # remove multiples
        varr = np.array(val)
        coords = list(np.average(varr,axis=0)[:3])
        err = list(np.sqrt(np.sum(np.square(varr),axis=0))[3:])
        return tuple(coords+err+list(idset)+list(mis))


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
        mixlist2 = []
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
            # end on start
            for m,s2 in enumerate(slist):
                if self.condition(e,s2):
                    mixlist2.append((n,m))
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

        for entry in mixlist2:
            gg.add_edge('e.'+str(entry[0]),'s.'+str(entry[1]))

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
        tol = 1.0e-6
        for val in newedgelist:
            s = val[0]
#            rs=[]
#            for entry in s:
#                rs.append(round(entry,6))
            e = val[1]
#            re=[]
#            for entry in e:
#                re.append(round(entry,6))
            w = self.dist(s,e)
            if w>tol:
                ed.append((tuple(s),tuple(e),round(w,4)))
        return ed


    def find_path(self,start=(0.0,0.0,0.0),end=(0.0,0.0,0.0)):
        '''
        find the shortest path between start and end nodes on the graph
        using the euclidean metric instead of graph metric
        Output: list of nodes of shortest path and absolute distance
        '''
        if (start==end): # prevent lack of sensical nodes
            return [], 0.0

        #print 'find path, start, end: ',start,end
        if (len(nx.predecessor(self.graph,start,end)) > 0):# if reachable
            #path = nx.shortest_path(self.graph,start,end,'weight')
            path = nx.dijkstra_path(self.graph,start,end,'weight')
            pl = nx.shortest_path_length(self.graph,start,end,'weight')
            return path, pl
        else: # broken graph between start and end; no shortest path
            #print 'Broken graph case'
            psdict = nx.single_source_shortest_path(self.graph,start)
            sfromstart = []
            sfromend = []
            for k,entry in psdict.iteritems():
                sfromstart.append(k)
            pedict = nx.single_source_shortest_path(self.graph,end)
            for k,entry in pedict.iteritems():
                sfromend.append(k)
            breakpoints = []
            for k in sfromstart:
                for l in sfromend:
                    stnodlist = psdict[k]
                    endnodlist = pedict[l]
                    if len(stnodlist)>1 and len(endnodlist)>1:
                        for stnod in stnodlist:
                            for endnod in endnodlist:
                                w = self.dist(stnod,endnod)
                                #print 'breakpoint distance: %f'%w
                                if w<19.0: # radius error to arc length 
                                    breakpoints.append((stnod,endnod,w))
            if len(breakpoints)>0:
                self.graph.add_weighted_edges_from(breakpoints)
                if (len(nx.predecessor(self.graph,start,end)) > 0):
                    path = nx.dijkstra_path(self.graph,start,end,'weight')
                    pl = nx.shortest_path_length(self.graph,start,end,'weight')
                    return path, pl
            else:
                return [],0.0
            

