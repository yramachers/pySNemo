from pysnemo.io.edm import tracker_hit
import ROOT
import math

ROOT.gROOT.ProcessLine(
"struct ListStruct{\
   int evtId;\
   vector<int>*    geiger_id;\
   vector<int>*    true_id;\
   vector<int>*    tracker_module;\
   vector<int>*    tracker_side;\
   vector<int>*    tracker_layer;\
   vector<int>*    tracker_column;\
   vector<double>* geiger_x;\
   vector<double>* geiger_y;\
   vector<double>* geiger_z;\
   vector<double>* geiger_sigma_z;\
   vector<double>* geiger_r;\
   vector<double>* geiger_sigma_r;\
};");



class listToRoot(object):
	'''
	Stream tracker_hit list output stored in event to a Root file.
	Can also read the list structure back from file.
	keystring: key in event containing a list of tracker_hits. 
	'''
	def __init__(self, keystring):
            self.key = keystring # force key to access listed data
	    self.listStruct = ROOT.ListStruct()
            self.builtTreeHeader = False
            self.totNumEvents = 0
                
        def __repr__(self):
            s = 'List to Root streamer'
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
		treename = "list_tree"
                if (not outputFile.Get(treename)):
			self.tree = ROOT.TTree(treename,"List data")
			self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get(treename)
                    self.tree.SetDirectory(self.outputFile)
                self.setupLists()


        def fillTree(self, event):

            self.outputFile.cd()
            if self.builtTreeHeader is False:
                origin = event.getFilename()            
                self.buildTreeHeader(origin, self.totNumEvents)

            self.event_id = event.getEventID()

            if (event.hasKey(self.key)):
                listdata = event.getKeyValue(self.key)
                self.StoreLists(listdata)
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
		k = ROOT.TObjString("list_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def setupLists(self):
		self.listStruct.true_id = ROOT.std.vector('int')()
		self.listStruct.tracker_module = ROOT.std.vector('int')()
		self.listStruct.tracker_side = ROOT.std.vector('int')()
		self.listStruct.tracker_layer = ROOT.std.vector('int')()
		self.listStruct.tracker_column = ROOT.std.vector('int')()
		self.listStruct.geiger_id = ROOT.std.vector('int')()
		self.listStruct.geiger_x = ROOT.std.vector('double')()
		self.listStruct.geiger_y = ROOT.std.vector('double')()
		self.listStruct.geiger_z = ROOT.std.vector('double')()
		self.listStruct.geiger_sigma_z = ROOT.std.vector('double')()
		self.listStruct.geiger_r = ROOT.std.vector('double')()
		self.listStruct.geiger_sigma_r = ROOT.std.vector('double')()

		
		self.tree.Branch('evtId', ROOT.AddressOf(self.listStruct, 'evtId'), 'evtId/I')

		self.tree.Branch('true.id', self.listStruct.true_id)
		self.tree.Branch('tracker.module', self.listStruct.tracker_module)
		self.tree.Branch('tracker.side', self.listStruct.tracker_side)
		self.tree.Branch('tracker.layer', self.listStruct.tracker_layer)
		self.tree.Branch('tracker.column', self.listStruct.tracker_column)
		self.tree.Branch('geiger.id', self.listStruct.geiger_id)
		self.tree.Branch('geiger.x', self.listStruct.geiger_x)
		self.tree.Branch('geiger.y', self.listStruct.geiger_y)
		self.tree.Branch('geiger.z', self.listStruct.geiger_z)
		self.tree.Branch('geiger.sigmaz', self.listStruct.geiger_sigma_z)
		self.tree.Branch('geiger.r', self.listStruct.geiger_r)
		self.tree.Branch('geiger.sigmar', self.listStruct.geiger_sigma_r)


	def StoreLists(self, cldata):
		# Store hits in output tree
		self.listStruct.evtId = self.event_id
		
		self.listStruct.true_id.clear()
		self.listStruct.tracker_module.clear()
		self.listStruct.tracker_side.clear()
		self.listStruct.tracker_layer.clear()
		self.listStruct.tracker_column.clear()
		self.listStruct.geiger_id.clear()
		self.listStruct.geiger_x.clear()
		self.listStruct.geiger_y.clear()
		self.listStruct.geiger_z.clear()
		self.listStruct.geiger_sigma_z.clear()
		self.listStruct.geiger_r.clear()
		self.listStruct.geiger_sigma_r.clear()

		for gg in cldata:
			mi = gg.meta_info
			self.listStruct.geiger_id.push_back(mi[0])
			self.listStruct.true_id.push_back(mi[1])
			self.listStruct.tracker_module.push_back(mi[2])
			self.listStruct.tracker_side.push_back(mi[3])
			self.listStruct.tracker_layer.push_back(mi[4])
			self.listStruct.tracker_column.push_back(mi[5])
			self.listStruct.geiger_x.push_back(gg.x)
			self.listStruct.geiger_y.push_back(gg.y)
			self.listStruct.geiger_z.push_back(gg.z)
			self.listStruct.geiger_sigma_z.push_back(gg.sigmaz)
			self.listStruct.geiger_r.push_back(gg.r)
			self.listStruct.geiger_sigma_r.push_back(gg.sigmar)
			
			
		self.tree.Fill()


	def ReturnList(self, tree, index):

		cls = self.listStruct
		cls.true_id = ROOT.std.vector('int')()
		cls.tracker_module = ROOT.std.vector('int')()
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
		tree.SetBranchAddress('true.id', ROOT.AddressOf(cls, "true_id") )
		tree.SetBranchAddress('tracker.module', ROOT.AddressOf(cls, "tracker_module") )
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
		for i in listOfEntries:
			
			tree.GetEntry(i) # cls is filled
			
			data = []
			for j in range(cls.geiger_x.size()):
				x = cls.geiger_x[j]
				y = cls.geiger_y[j]
				z = cls.geiger_z[j]
				dz = cls.geiger_sigma_z[j]
				r = cls.geiger_r[j]
				dr = cls.geiger_sigma_r[j]
				gg = tracker_hit(x,y,z,dz,r,dr)
				if math.isnan(r):
					r = 23.0 # hard coded radius, dr for 
					dr = 1.0 # peripheral hits 
				
				gid = cls.geiger_id[j]
				id = cls.true_id[j]
				m = cls.tracker_module[j]
				s = cls.tracker_side[j]
				l = cls.tracker_layer[j]
				c = cls.tracker_column[j]
				gg.set_info(gid,id,m,s,l,c)
				data.append(gg)

		return data

