import pysnemo.utility.gg2image as gg2i
from scipy.ndimage import label
import numpy as np
import networkx as nx


class ImageSegmentation(object):
    '''
    Simple tracker pre-clusterer using the tracker wires in each half tracker
    as black/white image pixels (on/off). The scipy ndimage routine 'label' 
    then checks connectedness of structures in that image. Each connected
    structure receives a distinct numerical label. 
    The segmentation
    algorithm then simply checks for splits in these separate structures
    to identify N>1 clusters in connected structures.
    Non-splitting single structures are recovered perfectly.
    Single split structures are identified and treated very 
    conservatively such that true hits in clusters are never lost but
    clusters rather contain too many hits.
    Multiple split structures are also dealt with conservatively but 
    are not necessarily very well split. These have to be re-clustered again 
    with a more capable clusterng algorithm.
    '''
    def __init__(self, data):
        '''
        Input: Ring data as list of tracker_hit objects.

        Output: Dictionary of lists of tracker_hit objects available with get method.
        '''
        self.tracker_layers = 9   # hardcoded tracker wire chamber image size
        self.tracker_rows  = 113
        self.data = data # needed for conversions to/from image



    def run(self):
        '''
        start and steer the processing of raw tracker data to clusters.
        '''
        # image data, take images from left and right half-tracker as ndarrays
        left, right = gg2i.gg_to_image(self.data, self.tracker_rows, self.tracker_layers)

        cll, clr = self.clusterer(left, right) # process the tracker images, return clusters

        self.tracker_clusters = gg2i.images_to_gg(self.data, cll, clr) # returns a dictionary


        
    def getClusters(self):
        return self.tracker_clusters # clusters in form of a dictionary



    def clusterer(self, left, right):
        '''
        main clustering routine, steers what is to happen with each half-tracker.
        '''
        clL = { } # results storage
        clR = { }
        all_left = []
        all_right = []
        row,col = left.shape # same as right.shape hence only once

        anyslope = [[1, 1, 1], [1,1,1], [1,1,1]] # any slope structure


        # check left tracker
        leftlabels, anyleft = label(left,anyslope)
        #print 'initial connected structures left = ',anyleft
        if anyleft>0:
            leftcollection = self._all_connected(leftlabels, anyleft) # get collection of all singly connected segments
            for ncls, left in enumerate(leftcollection): # left is a full image with one structure contained according to label
                all_left.append(self._collectionsplitting(left,leftlabels,0,ncls+1)) # send clusters into list

        # check right tracker
        rightlabels, anyright = label(right,anyslope)
        #print 'initial connected structures right = ',anyright
        if anyright>0: 
            rightcollection = self._all_connected(rightlabels, anyright) # get collection of all singly connected segments
            for ncls, right in enumerate(rightcollection):
                all_right.append(self._collectionsplitting(right,rightlabels,1,ncls+1)) # send clusters into list

        # unravel the list of clusters
        counter = 0
        for cls in all_left:
            for n in cls.keys():
                clL[counter+n] = cls[n]
            counter += n

        counter = 0
        for cls in all_right:
            for n in cls.keys():
                clR[counter+n] = cls[n]
            counter += n

        return clL, clR




    def _all_connected(self, imagelabels, nlabels):
        '''
        Gives a collection of images for each separate connected structure, if any.
        Makes segmentation trivial for all trivial, non-split tracker structures.
        '''
        clusteridx = { }
        collection = []
        for n in range(1,nlabels+1):
            clusteridx[n] = np.where(imagelabels==n) # all indices in imagelabels
        for k in clusteridx.keys():
            copyimage = np.zeros(imagelabels.shape)
            for r,c in zip(clusteridx[k][0],clusteridx[k][1]):
                copyimage[r][c] = 1
            collection.append(copyimage)
        return collection



    def _collectionsplitting(self, data, labeldata, side, clslabel):
        '''
        Main internal routine testing the splitting of connected structures
        '''
        #multiflag = False
        row,col = data.shape # data should be whole image
        anyslope = [[1,1,1], [1,1,1], [1,1,1]] # any slope structure

        split, lstore =  self._splitpoints(data) # check for splits in a connected structure
        #print 'after splitpoints: ',split
        #print lstore
        # case studies
        if split:
            if (2,2) in lstore: # multi split of singly connected structure
                print 'complex clustering in GraphClustering'
                gc = GraphClustering(side)
                return gc.run(data)
            elif (1,2) in lstore and (2,1) in lstore: # crossing tracks
                print 'complex clustering in GraphClustering'
                gc = GraphClustering(side)
                return gc.run(data)

            # single split of connected structure
            for counter,(ll,rr) in enumerate(lstore): # loop over tuples
                if ll<1 or rr<1:
                    continue # don't consider empty data patches, next loop

                if ll<rr:   # find switch location (1,1)->(1,2)
                    loc = counter+1
                    #print 'found location switch at %d'%loc
                    single, nstruct = label(data[:,0:loc],anyslope)
                    #print 'label for single gives %d to %d col position in image: '%(nstruct,loc)
                    multi, nstruct2 = label(data[:,loc:],anyslope)
                    #print 'label for multi gives %d from %d col position in image: '%(nstruct2,loc)

                    if nstruct==1: # just single cluster to merge with split clusters
                        clsingle = self._cluster(single,nstruct,side,0)
                        clmulti = self._cluster(multi,nstruct2,side,loc)
                        cl = self._merge_single_to_multi(clsingle,clmulti)

                    elif nstruct>1: # more fragments to merge to many - re-cluster
                        clsingle = self._cluster(single,nstruct,side,0)
                        clmulti = self._cluster(multi,nstruct2,side,loc)
                        mergedcls = self._merge_single_to_multi(clsingle,clmulti)
                        cl = self._recluster(data,mergedcls,clmulti)

                    else: # can happen for smooth curved tracks
                        clsingle = {1:[]} # empty single, nothing to merge
                        clmulti = self._cluster(multi,nstruct2,side,loc)
                        cl = self._merge_single_to_multi(clsingle,clmulti)
                    break # end of loop in this case

                elif rr<ll:   # go to end and find switch location (2,1)->(1,1)
                    iter = counter
                    while (rr<ll and iter<len(lstore) and ll>0 and rr>0): # prevent empty data patches
                        ll,rr = lstore[iter] # check
                        iter += 1 # next tuple

                    loc = iter-1 # correct for start zero
                    # check final slot
                    ll,rr = lstore[-1]
                    if rr<ll: # no switch until the final column
                        loc = iter
                    #print 'found location switch at %d'%loc
                    single, nstruct = label(data[:,loc:],anyslope)
                    #print 'label for single gives %d from %d col position in image: '%(nstruct,loc)
                    multi, nstruct2 = label(data[:,0:loc],anyslope)
                    #print 'label for multi gives %d to %d col position in image: '%(nstruct2,loc)

                    if nstruct==1: # just single cluster to merge with split clusters
                        clsingle = self._cluster(single,nstruct,side,loc)
                        clmulti = self._cluster(multi,nstruct2,side,0)
                        cl = self._merge_single_to_multi(clsingle,clmulti)

                    elif nstruct>1: # more fragments to merge to many - re-cluster
                        clsingle = self._cluster(single,nstruct,side,loc)
                        clmulti = self._cluster(multi,nstruct2,side,0)
                        mergedcls = self._merge_single_to_multi(clsingle,clmulti)
                        cl = self._recluster(data,mergedcls,clmulti)

                    else: # can happen for smooth curved tracks
                        clsingle = {1:[]} # empty single, nothing to merge
                        clmulti = self._cluster(multi,nstruct2,side,0)
                        cl = self._merge_single_to_multi(clsingle,clmulti)
                    break # end of loop in this case

        else:     # simple clustering with labels
            cl = { }
            clusteridx = { }
            hitinfo = []
            clusteridx[1] = np.where(labeldata==clslabel) # all indices in imagelabels
            for r,c in zip(clusteridx[1][0],clusteridx[1][1]):
                if side<1:
                    mi = (side, r, 8-c) # left tracker
                else:
                    mi = (side, r, c) # right tracker
                hitinfo.append(mi)
            cl[1] = hitinfo # hitinfo structure for hit identification later
            #print 'single cluster[1] has',hitinfo
        return cl



    def _splitpoints(self, data):
        split = False 
        row,col = data.shape # data should be whole image
        labelstore = []
        anyslope = [[1, 1, 1], [1,1,1], [1,1,1]] # any slope structure

        for n in range(1,col):
            # check half-tracker, cut at each column
            leftsplit, nstruct1 = label(data[:,0:n],anyslope)
            rightsplit, nstruct2 = label(data[:,n:],anyslope)
            labelstore.append((nstruct1,nstruct2)) # store as tuple

            if nstruct1>0 and nstruct2>0: # have some data at all in that tracker section
                if nstruct2>nstruct1 or nstruct2<nstruct1: # something split
                    #print 'split %d to %d at column %d'%(nstruct1,nstruct2,n)
                    split = True

        return (split,labelstore)




    def _cluster(self, imagelabels, maxlabels, side, offset):
        '''
        image pixels with label to simplified internal tracker hits
        '''
        hitinfo = []
        clusters = { }
        clusteridx = { }
        #print 'shape image: ',imagelabels.shape
        for n in range(1,maxlabels+1):
            clusteridx[n] = np.where(imagelabels==n) # all indices in imagelabels

        # re-order in tuples for translation to tracker hits
        for n in clusteridx.keys():
            for r,c in zip(clusteridx[n][0],clusteridx[n][1]):
                if side<1:
                    mi = (side, r, 8-c-offset) # left tracker
                else:
                    mi = (side, r, c+offset) # right tracker
                hitinfo.append(mi)
            clusters[n] = hitinfo # hitinfo structure for hit identification later
            #print 'cluster[%d]'%n,' has',hitinfo
            hitinfo = []
        return clusters




    def _merge_single_to_multi(self, single, multi): # expect cluster dictionaries as input
        '''
        Conservatively merge all the single clusters with all separate structures
        such that no true hit in a cluster can be lost. This is where too many 
        hits are inserted in clusters, on purpose. This needs to be sorted later
        in dedicated algorithms.
        '''
        newcls = { } # not overwriting the input
        newlist = []
        copylist = []
        #    print 'in merge_, have %d single cls, %d multi cls.'%(len(single),len(multi))
        for k,v in multi.iteritems():
            for entry in v:
                copylist.append(entry)
            newcls[k] = copylist # copy
            for n in single.keys(): # should be just one but could be two
                for entry in single[n]:
                    newlist.append(entry) # singles entries, merge with all
                newcls[k].extend(newlist)
                newlist = [] # reset
            copylist = [] # clear
            #    for k in newcls.keys():
            #        print 'merged cluster[%d]'%k,' has',newcls[k]
        return newcls



                
    def _recluster(self, data, merged, multi):
        row, col = data.shape
        anyslope = [[1,1,1], [1,1,1], [1,1,1]] # any slope structure
        clusteridx = { }
        compareset1 = set()
        compareset2 = set()
        newcls = { }
        hitinfo = []

        for k in merged.keys(): # for every merged cluster
            choice = 0 # no sub-collection chosen
            copyimage = np.zeros((row,col))
            for s,r,c in merged[k]:
                copyimage[r][c] = 1 # image to be split from merged
            for s,r,c in multi[k]:
                compareset1.add((r,c)) # reference set with multi hits, same size as merged
            #        print 'set to compare to from multi is: ',compareset1
            side = s # doesn't change later
            labeldata, nlabels = label(copyimage,anyslope)  # check for non-connected structures
            #print 're-cluster labels has %d structures.'%nlabels
            if nlabels>1: # separation of previously merged structures possible
                for n in range(1,nlabels+1):
                    clusteridx[n] = np.where(labeldata==n) # all indices in labeldata
                for k2 in clusteridx.keys():
                    for r2,c2 in zip(clusteridx[k2][0],clusteridx[k2][1]):
                        compareset2.add((r2,c2))
                    #                print 'compare set 2 from image is: ',compareset2
                    overlap = compareset1.intersection(compareset2)
                    #                print 'sets overlapping: ',overlap
                    if len(overlap): # not empty intersection
                        choice = k2 # split merged and multi overlap - split chosen
                    compareset2 = set()

                if choice>0: # identified a sub cluster
                    for r2,c2 in zip(clusteridx[choice][0],clusteridx[choice][1]):
                        hitinfo.append((side,r2,c2))
                    newcls[k] = hitinfo
            else:
                newcls[k] = merged[k] # otherwise no reclustering, i.e. copy
            compareset1 = set() # reset
            clusteridx = { }
            hitinfo = []

            #    for k in newcls.keys():
            #        print 'merged cluster[%d]'%k,' has',newcls[k]
        return newcls # re-clustered merged


        

    def _stitch_ends(self, l, m, r):
        '''
        stitching only for the messy case of multi structures merging with multi structures
        again not to loose any true hit in a cluster.
        '''
        w_to_w = set()
        for kl in l.keys():
            for entry in l[kl]: # triplet in list
                s = entry[0]
                row = entry[1]
                column = entry[2]
                for km in m.keys():
                    for i in range(-1,2):# check 3x3 block around 
                        for j in range(-1,2): # the border pixel
                            if (s,row+i,column+j) in m[km]:
                                w_to_w.add((kl,km))
        for tup in w_to_w:
            #print 'left: cluster %d merges with middle cluster %d'%tup
            for entry in l[kl]:
                m[km].append(entry) # form extended cluster for middle clusters from borders

        w_to_w = set()
        for kr in r.keys():
            for entry in r[kr]: # triplet in list
                s = entry[0]
                row = entry[1]
                column = entry[2]
                for km in m.keys():
                    for i in range(-1,2):# check 3x3 block around 
                        for j in range(-1,2): # the border pixel
                            if (s,row+i,column+j) in m[km]:
                                w_to_w.add((kr,km))
        for tup in w_to_w:
            #print 'right: cluster %d merges with middle cluster %d'%tup
            for entry in r[kr]:
                m[km].append(entry) # form extended cluster for middle clusters from borders
        return m




class GraphClustering(object):
    '''
    Tracker image clusterer using the tracker wires in each half tracker
    as black/white image pixels (on/off). Used for complicated track
    structures that go beyond the Image segmentation clusterer.
    '''
    def __init__(self, side):
        '''
        Input: Image data as np.array from Imagesegmentation object.
               Side as in which tracker half the image comes from.

        Output: consistent with collectionsplitting in Imagesegmentation
                a cluster of hitinfo tuples
        '''
        self.side = side


    def run(self, data):
        '''
        Input: Image data as np.array
        Returns clusters
        '''
        store = []
        for col in range(data.shape[1]): # all columns
            nlist = np.where(data[:,col]>0)
            if len(nlist[0]):
                cluster1d = self._oneDcluster(nlist[0])
                store.append(cluster1d)
        cl = self._cluster_withgraph(store)
        print 'in run: before translate:',cl
        return self._translate(store, cl)
        


    def _cluster_withgraph(self, data):
        edges = self._connections(data) # list of pairs of tuples with connections    
        G = nx.Graph()
        G.add_edges_from(edges)

        endslist = self._find_all_deadends(G)
        start, end = self._find_startfinish(edges)
        if len(endslist)>0:
            start.extend(endslist)
            end.extend(endslist)
        print 'Start nodes: ',start
        print 'End nodes: ',end

        cl = { } # final storage
        for s in start:
            for e in end:
                if nx.has_path(G,source=s, target=e):
                    if s != e and s[0] != e[0]: # exclude nodes and identical layers
                        cl[(s,e)] = [p for p in nx.all_shortest_paths(G, source=s, target=e)]
        return cl


    def _translate(self, store, cl):
        # relate graph clusters back to pixel hits
        clusters = { }
        clid = 1
        for k in cl.keys(): # keys are irrelevant, need cluster number as id
            for path in cl[k]: # list of shortest paths in cluster
                collection = []
                for entry in path: # now have nodes
                    layercluster = store[entry[0]-1][entry[1]] # the one-D cluster of pixels
                    for row in layercluster:
                        if self.side<1: # left case
                            collection.append((self.side, row, 10-entry[0]-1))
                        else:
                            collection.append((self.side, row, entry[0]-1))
                clusters[clid] = collection # ready for images to gg translation
                clid += 1
        #print 'in translate: clusters',clusters
        return clusters


    def _oneDcluster(self, data):
        store = []
        row = data[0]
        collection = []
        for entry in data: # array element entry
            if entry-row < 2: # nearest neighbour, ascending order assumed
                collection.append(entry)
                row = entry
            else: # distance>1 = next cluster
                store.append(collection)
                collection = [entry] # first entry in new cluster
                row = entry
        store.append(collection) # safe the final cluster collection
        return store


    def _is_connected(self, alist,blist):
        aset = set(alist)
        asetp1 = set(np.array(alist)+1) # move one unit
        asetm1 = set(np.array(alist)-1) # move one unit
        bset = set(blist)
        if len(aset.intersection(bset))>0: # direct overlap
            return True
        elif len(asetp1.intersection(bset))>0: # neighbour overlap
            return True
        elif len(asetm1.intersection(bset))>0: # neighbour overlap
            return True
        else:
            return False


    def _connections(self, data): # connect 1d cluster lists
        edges = []
        layer = data[0] # first list of nodes
        counter = 2 # second layer, start counting at 1 for layers
        for entry in data[1:]: # next layers with nodes
            for node in layer:
                for nextnode in entry:
                    if self._is_connected(node, nextnode):
                        edges.append(((counter-1, layer.index(node)),(counter, entry.index(nextnode))))
            layer = entry # prepare for next layer
            counter += 1 # next layer
        return edges



    def _find_all_deadends(self, graph):
        ends = []
        for node in graph.nodes():
            if len(graph.neighbors(node)) < 2: # single node end
                if node[0] != 1 and node[0] != 9: # known from start finish finder
                    ends.append(node)
        print 'dead ends found: ', ends
        return ends


    def _find_startfinish(self, edges):
        # find layer=1 and layer=9, the tracker extremes
        starts = set()
        ends = set()
        for s,e in edges: # a list
            if s[0]==1:
                starts.add(s)
            if e[0]==9:
                ends.add(e)
        return list(starts), list(ends)
