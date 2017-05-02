from scipy.ndimage import label
import numpy as np
import pysnemo.utility.gg2image as gg2i


class FilamentCluster(object):
    '''
    Simple tracker pre-clusterer using the tracker wires in each half tracker
    as black/white image pixels (on/off). The scipy ndimage routine 'label' 
    then checks connectedness of structures in that image. Each connected
    structure receives a distinct numerical label. 
    Non-splitting single structures are recovered perfectly.

    If a split in the structure is detected then the segmentation 
    algorithm checks for filaments: 
    Turn structure into a graph, find maximum connected nodes, remove them
    then check again for connectedness of nodes larger than 3, see whether 
    the number of image structures has increased, i.e. filaments were split off.
    If so, stop and declare each filament a cluster with the removed 
    nodes, i.e. image pixels, returned equally to both split structures
    (not clear to whom it belongs, hence to both).
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
        left, right = gg2i.gg_to_left_right_images(self.data, self.tracker_rows, self.tracker_layers)

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
        cl = { }
        pixelstore = [] # complete store of pixels removed and their neighbours
        anyslope = [[1,1,1], [1,1,1], [1,1,1]] # any slope structure

        split =  self._splitpoints(data) # check for splits in a connected structure
        #print 'after splitpoints: ',split
        # case studies: split=0; no split, 
        # =1; simple, use image segmentation, 
        # =2; complex, try filament search
        if split==2:
            '''
            This thins out filament structures in the image by removing 
            maxiumum connected pixels iteratively, leaving isolated 
            structures for the filament clusters. Needs to insert
            the missing pixels afterwards where suitable.
            '''
            broken = False
            #print 'Complex split, trying filament search.'

            # every node is a pixel in the image
            graph = gg2i.image_to_graph(data)
            while not broken:
                remove = [] # max connected node store to remove
                neighbours = [] # remaining neighbour nodes stored
                maxdegree = 0
                for nd in graph.nodes():
                    deg = graph.degree(nd)
                    if deg >= maxdegree:
                        maxdegree = deg
                #print 'maximum degree found: %d'%maxdegree
                if maxdegree < 5: # filaments too thin to find meaningful splitpoints
                    #print 'Filaments too thin, trying image segmentation.'
                    iseg = ImageSegmentation(data, side) # try alternativ with one-side data set
                    iseg.run()
                    cl = iseg.getClusters()
                    return cl # out of function
                    
                # collect nodes to remove
                for nd in graph.nodes():
                    if graph.degree(nd) == maxdegree:
                        remove.append(nd)
                        neighbours.append(graph.neighbors(nd)) # list of list of tuples
                    
                        # remove the nodes with maximum degree
                pixelstore.append((remove, neighbours)) # to put back later
                for nd in remove:
                    graph.remove_node(nd)

                im = gg2i.graph_to_image(graph, self.tracker_rows, self.tracker_layers)
                # check structures
                imagelabels, nlabels = label(im,anyslope)
                if nlabels > 1:
                    cls = self._cluster(imagelabels, nlabels, side)
                    #print 'n labels found: ',nlabels
                    broken = True
            # now put the undecided pixels back into the clusters
            cl = self._restore_pixels(pixelstore, cls)

        elif split==1:
            #print 'Simple split, trying image segmentation.'
            iseg = ImageSegmentation(data, side) # try with one-side data set
            iseg.run()
            cl = iseg.getClusters()
            return cl # out of function

        else:     # simple clustering with labels
            cl = { }
            clusteridx = { }
            hitinfo = []
            clusteridx[1] = np.where(labeldata==clslabel) # all indices in imagelabels
            for r,c in zip(clusteridx[1][0],clusteridx[1][1]):
                mi = (side, r, c) # left tracker
                hitinfo.append(mi)
            cl[1] = hitinfo # hitinfo structure for hit identification later
            #print 'Filament: single cluster[1] has',hitinfo
        return cl




    def _splitpoints(self, data):
        split = 0 # integer cases
        splitstore = split
        row, col = data.shape
        anyslope = [[1, 1, 1], [1,1,1], [1,1,1]] # any slope structure

        for n in range(1,col):
            # check half-tracker, cut at each column
            leftsplit, nstruct1 = label(data[:,0:n],anyslope)
            rightsplit, nstruct2 = label(data[:,n:],anyslope)

            # check cases
            if nstruct1>0 and nstruct2>0: # have some data at all in that tracker section
                if nstruct1==1 and nstruct2==2:
                    split = 1 # simple case split, use image segmentation
                elif nstruct1==2 and nstruct2==1:
                    split = 1 # also vice versa
                elif nstruct1>1 and nstruct2>1: # complex tracker data
                    split = 2 # try filament search instead
                if split > splitstore:
                    splitstore = split # most complex case prevails
                #print 'Fil: split (%d, %d) gives store=%d'%(nstruct1,nstruct2,splitstore)

        return splitstore




    def _cluster(self, imagelabels, maxlabels, side):
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
                mi = (side, r, c) # tracker
                hitinfo.append(mi)
            clusters[n] = hitinfo # hitinfo structure for hit identification later
            #print 'cluster[%d]'%n,' has',hitinfo
            hitinfo = []
        return clusters




    def _restore_pixels(self, pixelstore, cls):
        for tup in pixelstore: # separate removal loops
            pixels = tup[0]
            nnlist = tup[1]

            for k, v in cls.iteritems(): # check each cluster for content
                side = v[0][0] # same for all hitinfo entries
                dubletlist = [(mi[1], mi[2]) for mi in v] # where v is the hitinfo
                for dublet in dubletlist:
                    for nd, nn in zip(pixels, nnlist):
                        if dublet in nn: # at least one of the neighbours survived the pixel removal
                            hit = (side, nd[0], nd[1])
                            if not (hit in v):
                                v.append(hit) # insert node as hitinfo
                                #print 're-inserted info: (%d, %d, %d) in cluster %d'%(side, nd[0], nd[1], k)
        return cls



class ImageSegmentation(object):
    '''
    Alternativ image segmentation to Filament Clusterer for structure
    too thinly connected to take apart as filaments.
    Here the scipy ndimage routine 'label' 
    checks connectedness of structures in that image. Each connected
    structure receives a distinct numerical label. 
    The segmentation
    algorithm then simply checks for splits in these separate structures
    to identify N>1 clusters in connected structures.
    Single split structures are identified and treated very 
    conservatively such that true hits in clusters are never lost but
    clusters rather contain too many hits.
    '''
    def __init__(self, data, side):
        '''
        Input: Ring data as list of tracker_hit objects.

        Output: Dictionary of lists of tracker_hit objects available with get method.
        '''
        self.tracker_layers = 9   # hardcoded tracker wire chamber image size
        self.tracker_rows  = 113
        self.data = data # needed for conversions to/from image
        self.side = side


    def run(self):
        '''
        start and steer the processing of raw tracker data to clusters.
        '''
        # image data, take images from left and right half-tracker as ndarrays
        if isinstance(self.data, np.ndarray): # image received from Filament searcher
            self.tracker_clusters = self.clusterer(self.data) # process the tracker images, return clusters

        else: # full data segementation needed
            image = gg2i.gg_to_single_image(self.data, self.tracker_rows, self.tracker_layers)
            
            self.tracker_clusters = self.clusterer(image) # process the tracker images, return clusters


        
    def getClusters(self):
        return self.tracker_clusters # clusters in form of a dictionary



    def clusterer(self, image):
        '''
        main clustering routine, steers what is to happen with each half-tracker.
        '''
        cls = { } # results storage
        all_side = []

        anyslope = [[1, 1, 1], [1,1,1], [1,1,1]] # any slope structure


        # check left tracker
        labels, anyonside = label(image,anyslope)

        if anyonside>0:
            collection = self._all_connected(labels, anyonside) # get collection of all singly connected segments
            for ncls, onside in enumerate(collection): # left is a full image with one structure contained according to label
                all_side.append(self._collectionsplitting(onside,labels,self.side,ncls+1)) # send clusters into list

        # unravel the list of clusters
        counter = 0
        for cluster in all_side:
            for n in cluster.keys():
                cls[counter+n] = cluster[n]
            counter += n

        return cls




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
        anyslope = [[1,1,1], [1,1,1], [1,1,1]] # any slope structure

        split, lstore =  self._splitpoints(data) # check for splits in a connected structure
        #print 'after splitpoints: ',split
        #print lstore
        # case studies
        if split:
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
                mi = (side, r, c) # left tracker
                hitinfo.append(mi)
            cl[1] = hitinfo # hitinfo structure for hit identification later
            #print 'ImSeg: single cluster[1] has',hitinfo
        return cl



    def _splitpoints(self, data):
        split = False 
        row, col = data.shape
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
                mi = (side, r, c+offset) # left tracker
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
        #print 'in merge_, have %d single cls, %d multi cls.'%(len(single),len(multi))
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
            #for k in newcls.keys():
            #print 'merged cluster[%d]'%k,' has',newcls[k]
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
            
        return newcls # re-clustered merged


        

