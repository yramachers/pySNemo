from pysnemo.io.edm import gcylinder
import ROOT

ROOT.gROOT.ProcessLine(
"struct ClusterStruct{\
   int evtId;\
   int clId;\
   vector<int>*    tracker_id;\
   vector<int>*    tracker_side;\
   vector<int>*    tracker_layer;\
   vector<int>*    tracker_column;\
   vector<int>*    geiger_id;\
   vector<double>* geiger_x;\
   vector<double>* geiger_y;\
   vector<double>* geiger_z;\
   vector<double>* geiger_sigma_z;\
   vector<double>* geiger_r;\
   vector<double>* geiger_sigma_r;\
};");



class clusterToRoot(object):
	'''
	Stream cluster output stored in event to a Root file.
	Can also read the cluster structure back from file.
	keystring: key in event containing a dictionary format 
	           cluster structure (clusterID: list of tuples)
	'''
	def __init__(self, keystring):
            self.key = keystring # force key to access clustered data
	    self.clusterStruct = ROOT.ClusterStruct()
            self.builtTreeHeader = False
            self.totNumEvents = 0
                
        def __repr__(self):
            s = 'Cluster to Root streamer'
            return s

        def createTree(self, outputFile, totNumEvents):
            '''
            Create the TTree for output
            '''

            self.totNumEvents = totNumEvents

            if (not outputFile.IsOpen()):
                print 'ROOT file is not open'
                return None
            else:
                self.outputFile = outputFile
                self.outputFile.cd()
		treename = "cluster_tree"
                if (not outputFile.Get(treename)):
			self.tree = ROOT.TTree(treename,"Cluster data")
			self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get(treename)
                    self.tree.SetDirectory(self.outputFile)
                self.setupClusters()


        def fillTree(self, event):

            self.outputFile.cd()
            if self.builtTreeHeader is False:
                origin = event.getFilename()            
                self.buildTreeHeader(origin, self.totNumEvents)

            self.event_id = event.getEventID()

            if (event.hasKey(self.key)):
                clusterdict = event.getKeyValue(self.key)
                self.StoreClusters(clusterdict)
            else:
                print "Error in pipeline: key (%s) not in event" % self.key


        def writeTree(self):

            self.outputFile.cd()
            self.tree.Write()


	def buildTreeHeader(self, origin,nevents):
            
            if self.builtTreeHeader is False:

                tlist = self.tree.GetUserInfo()
		map = ROOT.TMap() # holds all key,value pairs: only one here
		map.SetName("Event Data Header")

		k = ROOT.TObjString("raw")# where the raw data is
		v = ROOT.TObjString(origin)
		map.Add(k,v)
		k = ROOT.TObjString("cluster_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def setupClusters(self):
		self.clusterStruct.tracker_id = ROOT.std.vector('int')()
		self.clusterStruct.tracker_side = ROOT.std.vector('int')()
		self.clusterStruct.tracker_layer = ROOT.std.vector('int')()
		self.clusterStruct.tracker_column = ROOT.std.vector('int')()
		self.clusterStruct.geiger_id = ROOT.std.vector('int')()
		self.clusterStruct.geiger_x = ROOT.std.vector('double')()
		self.clusterStruct.geiger_y = ROOT.std.vector('double')()
		self.clusterStruct.geiger_z = ROOT.std.vector('double')()
		self.clusterStruct.geiger_sigma_z = ROOT.std.vector('double')()
		self.clusterStruct.geiger_r = ROOT.std.vector('double')()
		self.clusterStruct.geiger_sigma_r = ROOT.std.vector('double')()

		
		self.tree.Branch('evtId', ROOT.AddressOf(self.clusterStruct, 'evtId'), 'evtId/I')
		self.tree.Branch('clId', ROOT.AddressOf(self.clusterStruct, 'clId'), 'clId/I')

		self.tree.Branch('tracker.id', self.clusterStruct.tracker_id)
		self.tree.Branch('tracker.side', self.clusterStruct.tracker_side)
		self.tree.Branch('tracker.layer', self.clusterStruct.tracker_layer)
		self.tree.Branch('tracker.column', self.clusterStruct.tracker_column)
		self.tree.Branch('geiger.id', self.clusterStruct.geiger_id)
		self.tree.Branch('geiger.x', self.clusterStruct.geiger_x)
		self.tree.Branch('geiger.y', self.clusterStruct.geiger_y)
		self.tree.Branch('geiger.z', self.clusterStruct.geiger_z)
		self.tree.Branch('geiger.sigmaz', self.clusterStruct.geiger_sigma_z)
		self.tree.Branch('geiger.r', self.clusterStruct.geiger_r)
		self.tree.Branch('geiger.sigmar', self.clusterStruct.geiger_sigma_r)


	def StoreClusters(self, cldict):
		for k,entry in cldict.iteritems():
			# Store hits in output tree
			self.clusterStruct.evtId = self.event_id
			self.clusterStruct.clId = k

			self.clusterStruct.tracker_id.clear()
			self.clusterStruct.tracker_side.clear()
			self.clusterStruct.tracker_layer.clear()
			self.clusterStruct.tracker_column.clear()
			self.clusterStruct.geiger_id.clear()
			self.clusterStruct.geiger_x.clear()
			self.clusterStruct.geiger_y.clear()
			self.clusterStruct.geiger_z.clear()
			self.clusterStruct.geiger_sigma_z.clear()
			self.clusterStruct.geiger_r.clear()
			self.clusterStruct.geiger_sigma_r.clear()

			for gg in entry:
				mi = gg.meta_info
				self.clusterStruct.geiger_id.push_back(mi[0])
				self.clusterStruct.tracker_id.push_back(mi[1])
				self.clusterStruct.tracker_side.push_back(mi[2])
				self.clusterStruct.tracker_layer.push_back(mi[3])
				self.clusterStruct.tracker_column.push_back(mi[4])
				self.clusterStruct.geiger_x.push_back(gg.xwire)
				self.clusterStruct.geiger_y.push_back(gg.ywire)
				self.clusterStruct.geiger_z.push_back(gg.zwire)
				self.clusterStruct.geiger_sigma_z.push_back(gg.dz)
				self.clusterStruct.geiger_r.push_back(gg.radius)
				self.clusterStruct.geiger_sigma_r.push_back(gg.dr)


			self.tree.Fill()


	def ReturnCluster(self, tree, index):

		cls = self.clusterStruct
		cls.tracker_id = ROOT.std.vector('int')()
		cls.tracker_side = ROOT.std.vector('int')()
		cls.tracker_layer = ROOT.std.vector('int')()
		cls.tracker_column = ROOT.std.vector('int')()
		cls.geiger_id = ROOT.std.vector('int')()
		cls.geiger_x = ROOT.std.vector('double')()
		cls.geiger_y = ROOT.std.vector('double')()
		cls.geiger_z = ROOT.std.vector('double')()
		cls.geiger_sigma_z = ROOT.std.vector('double')()
		cls.geiger_r = ROOT.std.vector('double')()
		cls.geiger_sigma_r = ROOT.std.vector('double')()

		N = tree.GetEntries()

		tree.SetBranchAddress('evtId', ROOT.AddressOf(cls, "evtId") )
		tree.SetBranchAddress('clId', ROOT.AddressOf(cls, "clId") )
		tree.SetBranchAddress('tracker.id', ROOT.AddressOf(cls, "tracker_id") )
		tree.SetBranchAddress('tracker.side', ROOT.AddressOf(cls, "tracker_side") )
		tree.SetBranchAddress('tracker.layer', ROOT.AddressOf(cls, "tracker_layer") )
		tree.SetBranchAddress('tracker.column', ROOT.AddressOf(cls, "tracker_column") )
		tree.SetBranchAddress('geiger.id', ROOT.AddressOf(cls, "geiger_id") )
		tree.SetBranchAddress('geiger.x', ROOT.AddressOf(cls, "geiger_x") )
		tree.SetBranchAddress('geiger.y', ROOT.AddressOf(cls, "geiger_y") )
		tree.SetBranchAddress('geiger.z', ROOT.AddressOf(cls, "geiger_z") )
		tree.SetBranchAddress('geiger.sigmaz', ROOT.AddressOf(cls, "geiger_sigma_z") )
		tree.SetBranchAddress('geiger.r', ROOT.AddressOf(cls, "geiger_r") )
		tree.SetBranchAddress('geiger.sigmar', ROOT.AddressOf(cls, "geiger_sigma_r") )

		cldict = {}

		# First, find out the tree indices which correspond to our
		# event id index. Set the status of all other branches to "off"
		tree.SetBranchStatus("*", 0)
		tree.SetBranchStatus("evtId", 1)

		# Find the entries that correspond to our eventId
		listOfEntries = []
		for i in range(N):
			tree.GetEntry(i)  # cls is filled
			if (cls.evtId == index):
				listOfEntries.append(i)

		# Re-enable all branches, looping over only those entries
		# that match our event index
		tree.SetBranchStatus("*", 1)
		for i in range(len(listOfEntries)):

			tree.GetEntry(listOfEntries[i]) # cls is filled

			data = []
			for j in range(cls.geiger_x.size()):
				x = cls.geiger_x[j]
				y = cls.geiger_y[j]
				z = cls.geiger_z[j]
				dz = cls.geiger_sigma_z[j]
				r = cls.geiger_r[j]
				dr = cls.geiger_sigma_r[j]
				gg = gcylinder(x,y,z,dz,r,dr)
				
				gid = cls.geiger_id[j]
				id = cls.tracker_id[j]
				s = cls.tracker_side[j]
				l = cls.tracker_layer[j]
				c = cls.tracker_column[j]
				gg.set_info(gid,id,s,l,c)
				data.append(gg)
			cldict[cls.clId] = data

		return cldict

