import numpy as np
import networkx as nx
from scipy.spatial import KDTree


class NNSetCluster(object):
    '''
    Nearest neighbour cluster algorithm from xu Zihe (Final Year project
    2014) using basic mono-directional stepping.
    '''
    def __init__(self,data,demonstrator_flag=True):
        '''
        Input: Ring data as list of tracker_hit objects and
               a geometry flag to distinguish between commissioning
               setup (False) and demonstrator (True).

        Output: List of tracker_hit objects and a list of 
                noise tracker_hit objects available with get methods.
        '''
        self.dflag = demonstrator_flag
        self.final_layer = 8

        if demonstrator_flag:
            self.data = sorted(data,key=lambda entry: abs(entry.x)) # sort in x
            self.metainfo = [gc.meta_info for gc in self.data] # total info 
            # split left side(=0) from right side(=1), index [3]
            self.left = []
            self.right = []
            for d in self.data:
                if d.meta_info[3]==0:
                    self.left.append(d)
                else:
                    self.right.append(d)
            # symmetric around x=0 with max at +-405mm
            self.infolist_l = [gc.meta_info for gc in self.left] # meta info tuple
            #print 'data left: ',self.infolist_l
            self.infolist_r = [gc.meta_info for gc in self.right] # meta info tuple
            #print 'data right: ',self.infolist_r
        else:
            sdata = sorted(data,key=lambda entry: entry.x) # sort in x
            sdata.reverse()
            # symmetric around x = 0 (top at +176mm)
            self.data = sdata
            self.metainfo = [gc.meta_info for gc in sdata] # info truth id

        # results container
        self.noise_hits = [] # hold list of tracker_hit objects
        self.track_clusters = [] # holds list of dict


    def getNoise(self):
        return self.noise_hits


    def _buildClusters(self, output, lrn):
        if lrn==0: # left tracker
            infoarray = np.array(self.infolist_l)
            layer_column = infoarray[:,4:]
            dubletlist = infoarray[:,4:].tolist()
            dublets = [tuple(l) for l in dubletlist] # total list

            result = []
            for entry in output: # list of dict
                cl = {}
                for k,v in entry.iteritems(): # three cluster types
                    collection = []
                    for hitset in v: # set of tuples
                        if isinstance(hitset,set):
                            items = []
                            for hit in hitset:
                                indx = dublets.index(hit)
                                items.append(self.left[indx])
                            collection.append(items)
                        else:
                            indx = dublets.index(hitset)
                            collection.append(self.left[indx])
                            
                    cl[k] = collection
                result.append(cl)
            return result
        elif lrn==1: # right tracker
            infoarray = np.array(self.infolist_r)
            layer_column = infoarray[:,4:]
            dubletlist = infoarray[:,4:].tolist()
            dublets = [tuple(l) for l in dubletlist] # total list

            result = []
            for entry in output: # list of dict
                cl = {}
                for k,v in entry.iteritems(): # three cluster types
                    collection = []
                    for hitset in v: # set of tuples
                        if isinstance(hitset,set):
                            items = []
                            for hit in hitset:
                                indx = dublets.index(hit)
                                items.append(self.right[indx])
                            collection.append(items)
                        else:
                            indx = dublets.index(hitset)
                            collection.append(self.right[indx])

                    cl[k] = collection
                result.append(cl)
            return result
        else: # commissioning
            infoarray = np.array(self.metainfo)
            layer_column = infoarray[:,4:]
            dubletlist = infoarray[:,4:].tolist()
            dublets = [tuple(l) for l in dubletlist] # total list
            
            result = []
            for entry in output: # list of dict
                cl = {}
                for k,v in entry.iteritems(): # three cluster types
                    collection = []
                    for hitset in v: # set of tuples
                        if isinstance(hitset,set):
                            items = []
                            for hit in hitset:
                                indx = dublets.index(hit)
                                items.append(self.data[indx])
                            collection.append(items)
                        else:
                            indx = dublets.index(hitset)
                            collection.append(self.data[indx])

                    cl[k] = collection
                result.append(cl)
            return result # list of dict returns


    def getClusters(self):
        return self.track_clusters


    def _setNoise(self):
        for entry in self.track_clusters:
            self.noise_hits.extend(entry[-1]) # the isolated hits


    def start_end(self, graph):
	starters = []
	for entry in graph.nodes():
            if entry[0]==0: # insist on foil
                starters.append(entry)
        #print 'starters 1: ',starters

	degdict = nx.degree(graph)
	# find dangling nodes
	elist_d = []
	elist_c = []
	for k,v in degdict.iteritems():
            if v==1 and k not in starters and not (k[0]==self.final_layer):
                elist_d.append(k)
            if k[0]==self.final_layer: # or calo nodes
                elist_c.append(k)
	#print 'ends dangling: ',elist_d
	#print 'ends  calo: ',elist_c
	return elist_d, elist_c, starters


    def cluster(self, elist, starters, gr):
	clusters=[]
	for end in elist:
		s=set()
		for entry in starters:
                    #print "start, end: ",entry,end
                    #print "graph nodes: ",gr.nodes()
                    if len(gr.nodes())<30:
                        path=nx.all_simple_paths(gr,entry,end)
                    else:
                        path=nx.all_shortest_paths(gr,entry,end)
                    for p in path:
                        s.update(set(p))
                clusters.append(s)
	return clusters


    def check_s2subsets1(self,s1,s2):
	s12 = s1 & s2
	if s12 == s2: # subset test
            return True # take superset
	else:
            return False
        

    def remove_doubles(self,ll):
	remove_set = set()
	for j,entry in enumerate(ll):
            for i in range(j+1,len(ll)):
                if len(ll[j].difference(ll[i]))<1:
                    remove_set.add(j)
        newll = []
	for k,entry in enumerate(ll):
            if k not in remove_set:
                newll.append(entry)
        return newll
    
    
    def component_cluster(self,graph):
        clusters = {}
	ed, ecalo, start = self.start_end(graph)
	if len(start)<1: # no foil hit
            clusters[0] = []
            clusters[1] = []
            clusters[-1] = graph.nodes()
	else:
            fcl = self.cluster(ed,start,graph)
            # index 0 for foil to dangling end clusters
            if len(fcl)>1:
                fcl2 = self.remove_doubles(fcl)
                #print 'removed doubles in fcl '
                newcl = []
                delset = set()
                for j,clset in enumerate(fcl2):
                    for i in range(j+1,len(fcl2)):
                        delflag1 = self.check_s2subsets1(fcl2[i],clset)
                        delflag2 = self.check_s2subsets1(clset,fcl2[i])
                        if delflag1:
                            delset.add(j)
                        elif delflag2:
                            delset.add(i)
                for j,clset in enumerate(fcl2):
                    if j not in delset:
                        newcl.append(clset)
                clusters[0] = newcl
                
            else:
                clusters[0] = fcl

            fccl = self.cluster(ecalo,start,graph)
            # index 1 for foil to calo end clusters
            if len(fccl)>1:
                fccl2 = self.remove_doubles(fccl)
                ##print 'removed doubles in fccl '
                newcl = []
                delset = set()
                for j,clset in enumerate(fccl2):
                    for i in range(j+1,len(fccl2)):
                        delflag1 = self.check_s2subsets1(fccl2[i],clset)
                        delflag2 = self.check_s2subsets1(clset,fccl2[i])
                        if delflag1:
                            delset.add(j)
                        elif delflag2:
                            delset.add(i)

                for j,clset in enumerate(fccl2):
                    if j not in delset:
                        newcl.append(clset)
                clusters[1] = newcl
            else:
                clusters[1] = fccl 
                
            #print 'foil clusters: ',clusters[0]
            #print 'foil, calo clusters: ',clusters[1]
            subdata = graph.nodes()
            for d,item in clusters.iteritems():
                for s in item:
                    setdata = list(s)
                    for entry in setdata:
                        if entry in subdata:
                            subdata.remove(entry)
            clusters[-1] = subdata # index -1 for non-foil data

        return clusters


    def run(self):
        if self.dflag:
            gr = nx.Graph()
            #print 'LEFT Tracker: '
            #print 'Have %d geiger cylinders'%len(self.infolist_l)
            if len(self.infolist_l)>0:
                infoarray = np.array(self.infolist_l)
                layer_column = infoarray[:,4:]
                lclist = infoarray[:,4:].tolist()
                kdt = KDTree(layer_column)
                pairs = kdt.query_pairs(1.5) #set of tuples of pair indices
                for pair in pairs:
                    gr.add_edge(tuple(lclist[pair[0]]),tuple(lclist[pair[1]]))	

                glist = nx.connected_component_subgraphs(gr)
                #print 'Connected components left: ',len(glist)

                cl1 = []
                for i,g in enumerate(glist):
                    cl1.append(self.component_cluster(g)) # list of dict

            else:
                cl1 = []

            gr.clear()
            #print 'RIGHT Tracker: '
            #print 'Have %d geiger cylinders'%len(self.infolist_r)
            if len(self.infolist_r)>0:
                infoarray = np.array(self.infolist_r)
                layer_column = infoarray[:,4:]
                lclist = infoarray[:,4:].tolist()
                kdt = KDTree(layer_column)
                pairs = kdt.query_pairs(1.5) #set of tuples of pair indices
                for pair in pairs:
                    gr.add_edge(tuple(lclist[pair[0]]),tuple(lclist[pair[1]]))	

                glist = nx.connected_component_subgraphs(gr)
                #print 'Connected components right: ',len(glist)

                cl2 = []
                for i,g in enumerate(glist):
                    cl2.append(self.component_cluster(g)) # list of dict

            else:
                cl2 = []

            if len(cl1)>0:
                collection1 = self._buildClusters(cl1,0)
            else:
                collection1 = []

            if len(cl2)>0:
                collection2 = self._buildClusters(cl2,1)
            else:
                collection2 = []
            self.track_clusters = collection1 + collection2 # concatenate
            self._setNoise()
            
        else:
            #print 'Have %d geiger cylinders'%len(self.metainfo)
            infoarray = np.array(self.metainfo)
            layer_column = infoarray[:,4:]
            lclist = infoarray[:,4:].tolist()
            kdt = KDTree(layer_column)
            pairs = kdt.query_pairs(1.5) #set of tuples of pair indices
            gr = nx.Graph()
            for pair in pairs:
                gr.add_edge(tuple(lclist[pair[0]]),tuple(lclist[pair[1]]))	

            glist = nx.connected_component_subgraphs(gr)
            #print 'Connected components in event: ',len(glist)

            output = []
            for i,g in enumerate(glist):
                output.append(self.component_cluster(g)) # list of dict

            if len(output)>0:
                self.track_clusters = self._buildClusters(output,2)
            else:
                self.track_clusters = []
            self._setNoise()
            
