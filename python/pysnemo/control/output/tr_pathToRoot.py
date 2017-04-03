import ROOT
import numpy as np
from pysnemo.io.edm import ClusterPaths

ROOT.gROOT.ProcessLine(
"struct DataStruct{\
   int evtId;\
   int clId;\
   int pathId;\
   vector<int>*    geiger_id;\
   vector<int>*    true_id;\
   vector<int>*    tracker_module;\
   vector<int>*    tracker_side;\
   vector<int>*    tracker_layer;\
   vector<int>*    tracker_column;\
   vector<double>* x;\
   vector<double>* y;\
   vector<double>* z;\
   vector<double>* ex;\
   vector<double>* ey;\
   vector<double>* ez;\
};");

class pathToRoot(object):
	'''
	Stream path dictionaries in cluster dictionaries
	stored in event to a Root file.
	Can also read the structure back from file.
	keystring: key in event containing a dictionary format clusters of
	           path candidates in list structure 
		   (list of tuples of (x,y,z,ex,ey,ez))
	'''
	def __init__(self, keystring):
            self.key = keystring # semi-force key to access clustered data
            self.builtTreeHeader = False
            self.totNumEvents = 0
	    self.dataStruct = ROOT.DataStruct()
	
	def __repr__(self):
            s = 'Path to Root streamer'
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
                if (not outputFile.Get("path_tree")):
                    self.tree = ROOT.TTree("path_tree","Path data")
                    self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get("path_tree")
                    self.tree.SetDirectory(self.outputFile)
                self.setupHits()

        
        def fillTree(self, event):

            self.outputFile.cd()
            if self.builtTreeHeader is False:
                origin = event.getFilename()            
                self.buildTreeHeader(origin, self.totNumEvents)

            self.event_id = event.getEventID()

	    if (event.hasKey(self.key)):
                    hits = event.getKeyValue(self.key)
		    self.StoreHits(hits)
	    else:
                    print "Error in pipeline: key (%s) not in event" % self.key
		    


        def writeTree(self):

            self.outputFile.cd()
            self.tree.Write()


	def setupHits(self):

		self.dataStruct.geiger_id = ROOT.std.vector('int')()
		self.dataStruct.true_id = ROOT.std.vector('int')()
		self.dataStruct.tracker_module = ROOT.std.vector('int')()
		self.dataStruct.tracker_side = ROOT.std.vector('int')()
		self.dataStruct.tracker_layer = ROOT.std.vector('int')()
		self.dataStruct.tracker_column = ROOT.std.vector('int')()
		self.dataStruct.x = ROOT.std.vector('double')()
		self.dataStruct.y = ROOT.std.vector('double')()
		self.dataStruct.z = ROOT.std.vector('double')()
		self.dataStruct.ex = ROOT.std.vector('double')()
		self.dataStruct.ey = ROOT.std.vector('double')()
		self.dataStruct.ez = ROOT.std.vector('double')()
		
		self.tree.Branch('evtId', ROOT.AddressOf(self.dataStruct, 'evtId'), 'evtId/I')
		self.tree.Branch('clId', ROOT.AddressOf(self.dataStruct, 'clId'), 'clId/I')
		self.tree.Branch('pathId', ROOT.AddressOf(self.dataStruct, 'pathId'), 'pathId/I')
		self.tree.Branch('geiger_id', self.dataStruct.geiger_id)
		self.tree.Branch('true_id', self.dataStruct.true_id)
		self.tree.Branch('tracker_module', self.dataStruct.tracker_module)
		self.tree.Branch('tracker_side', self.dataStruct.tracker_side)
		self.tree.Branch('tracker_layer', self.dataStruct.tracker_layer)
		self.tree.Branch('tracker_column', self.dataStruct.tracker_column)
		self.tree.Branch('x', self.dataStruct.x)
		self.tree.Branch('y', self.dataStruct.y)
		self.tree.Branch('z', self.dataStruct.z)
		self.tree.Branch('ex', self.dataStruct.ex)
		self.tree.Branch('ey', self.dataStruct.ey)
		self.tree.Branch('ez', self.dataStruct.ez)

		
	def buildTreeHeader(self, origin, nevents):

            if self.builtTreeHeader is False:

		tlist = self.tree.GetUserInfo()
		map = ROOT.TMap() # holds all key,value pairs: only one here
		map.SetName("Event Data Header")

		k = ROOT.TObjString("raw")# where the raw data is
		v = ROOT.TObjString(origin)
		map.Add(k,v)
		k = ROOT.TObjString("path_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def StoreHits(self, hits):
		for k,candidates in hits.iteritems():
			# Store hits in output tree
			self.dataStruct.evtId = self.event_id
			self.dataStruct.clId = k
			for l,paths in candidates.iteritems():
				self.dataStruct.pathId = l
				self.dataStruct.geiger_id.clear()
				self.dataStruct.true_id.clear()
				self.dataStruct.tracker_module.clear()
				self.dataStruct.tracker_side.clear()
				self.dataStruct.tracker_layer.clear()
				self.dataStruct.tracker_column.clear()
				self.dataStruct.x.clear()
				self.dataStruct.y.clear()
				self.dataStruct.z.clear()
				self.dataStruct.ex.clear()
				self.dataStruct.ey.clear()
				self.dataStruct.ez.clear()
		
				for entry in paths:
					mi = entry[-1]
					self.dataStruct.geiger_id.push_back(mi[0])
					self.dataStruct.true_id.push_back(mi[1])
					self.dataStruct.tracker_module.push_back(mi[2])
					self.dataStruct.tracker_side.push_back(mi[3])
					self.dataStruct.tracker_layer.push_back(mi[4])
					self.dataStruct.tracker_column.push_back(mi[5])
					self.dataStruct.x.push_back(entry[0])
					self.dataStruct.y.push_back(entry[1])
					self.dataStruct.z.push_back(entry[2])
					self.dataStruct.ex.push_back(entry[3])
					self.dataStruct.ey.push_back(entry[4])
					self.dataStruct.ez.push_back(entry[5])

				self.tree.Fill()


	def ReturnPaths(self, tree, index):
		ds = self.dataStruct
		ds.geiger_id = ROOT.std.vector('int')()
		ds.true_id = ROOT.std.vector('int')()
		ds.tracker_module = ROOT.std.vector('int')()
		ds.tracker_side = ROOT.std.vector('int')()
		ds.tracker_layer = ROOT.std.vector('int')()
		ds.tracker_column = ROOT.std.vector('int')()
		ds.x = ROOT.std.vector('double')()
		ds.y = ROOT.std.vector('double')()
		ds.z = ROOT.std.vector('double')()
		ds.ex = ROOT.std.vector('double')()
		ds.ey = ROOT.std.vector('double')()
		ds.ez = ROOT.std.vector('double')()

		N = tree.GetEntries()

		tree.SetBranchAddress('evtId', ROOT.AddressOf(ds, "evtId") )
		tree.SetBranchAddress('clId', ROOT.AddressOf(ds, "clId") )
		tree.SetBranchAddress('pathId', ROOT.AddressOf(ds, "pathId") )
		tree.SetBranchAddress('geiger_id', ROOT.AddressOf(ds, "geiger_id") )
		tree.SetBranchAddress('true_id', ROOT.AddressOf(ds, "true_id") )
		tree.SetBranchAddress('tracker_module', ROOT.AddressOf(ds, "tracker_module") )
		tree.SetBranchAddress('tracker_side', ROOT.AddressOf(ds, "tracker_side") )
		tree.SetBranchAddress('tracker_layer', ROOT.AddressOf(ds, "tracker_layer") )
		tree.SetBranchAddress('tracker_column', ROOT.AddressOf(ds, "tracker_column") )
		tree.SetBranchAddress('x', ROOT.AddressOf(ds, "x") )
		tree.SetBranchAddress('y', ROOT.AddressOf(ds, "y") )
		tree.SetBranchAddress('z', ROOT.AddressOf(ds, "z") )
		tree.SetBranchAddress('ex', ROOT.AddressOf(ds, "ex") )
		tree.SetBranchAddress('ey', ROOT.AddressOf(ds, "ey") )
		tree.SetBranchAddress('ez', ROOT.AddressOf(ds, "ez") )

		cldict = {}

		# First, find out the tree indices which correspond to our
		# event id index. Set the status of all other branches to "off"
		tree.SetBranchStatus("*", 0)
		tree.SetBranchStatus("evtId", 1)

		# Find the entries that correspond to our eventId
		listOfEntries = []
		for i in range(N):
			tree.GetEntry(i)  # self.dataStruct is filled
			if (ds.evtId == index):
				listOfEntries.append(i)

		# Re-enable all branches, looping over only those entries
		# that match our event index
		tree.SetBranchStatus("*", 1)

		clid = set()
		candidates = {}
		for i in range(len(listOfEntries)):
			tree.GetEntry(listOfEntries[i]) # self.dataStruct is filled
			clid.add(ds.clId)
			path = []
			for xc,yc,zc,ex,ey,ez,gid,tid,tm,ts,tl,tc in zip(ds.x,ds.y,ds.z,ds.ex,ds.ey,ds.ez,ds.geiger_id,ds.true_id,ds.tracker_module,ds.tracker_side,ds.tracker_layer,ds.tracker_column):
				minfo = (gid,tid,tm,ts,tl,tc)
				path.append((xc,yc,zc,ex,ey,ez,minfo))

			candidates[(ds.clId,ds.pathId)] = path

		result = []
		for entry in clid: # unique clId numbers
			cp = ClusterPaths(entry)
#			cand = {}
			for tup,p in candidates.iteritems():
				if tup[0]==entry:
					cp.add_path(tup[1], p)
#					cand[tup[1]] = p
#			cldict[entry] = cand
			result.append(cp)

#		return cldict
		return result # return list of ClusterPaths object
	

