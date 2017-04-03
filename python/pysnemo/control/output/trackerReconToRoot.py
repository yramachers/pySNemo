import ROOT
from numpy import array


ROOT.gROOT.ProcessLine(
"struct treconStruct{\
   int evtId;\
   int clId;\
   vector<int>* pathId;\
   vector<double>* fy;\
   vector<double>* fyup;\
   vector<double>* fydown;\
   vector<double>* fz;\
   vector<double>* fzup;\
   vector<double>* fzdown;\
   vector<double>* cy;\
   vector<double>* cyup;\
   vector<double>* cydown;\
   vector<double>* cz;\
   vector<double>* czup;\
   vector<double>* czdown;\
   vector<bool>* bpflag;\
   vector<int>* charge;\
   vector<vector<double> >* angles;\
};");

class trackerReconToRoot(object):
	'''
	Stream track reconstruction output stored in event to a Root file.
	keystring: key in event containing a list format 
	           hits structure (list of tuples)
	'''
	def __init__(self, keystring):
            self.key = keystring
            self.builtTreeHeader = False
            self.totNumEvents = 0
	    self.treconStruct = ROOT.treconStruct()
	
	def __repr__(self):
            s = 'Track Reconstruction to Root streamer'
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
                if (not outputFile.Get("trecon_tree")):
                    self.tree = ROOT.TTree("trecon_tree", "Tracker recon data")
                    self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get("trecon_tree")
                    self.tree.SetDirectory(self.outputFile)
                self.setupTValInfo()

        
        def fillTree(self, event):

            self.outputFile.cd()
            if self.builtTreeHeader is False:
                origin = event.getFilename()            
                self.buildTreeHeader(origin, self.totNumEvents)

            self.event_id = event.getEventID()

            if (self.key is not None):
                if (event.hasKey(self.key)):
                    tvalDict = event.getKeyValue(self.key)
		    self.StoreTValInfo(tvalDict)
                else:
                    print "Error in pipeline: key (%s) not in event" % self.key


        def writeTree(self):

            self.outputFile.cd()
            self.tree.Write()


	def setupTValInfo(self):

		self.treconStruct.pathId = ROOT.std.vector('int')()
		self.treconStruct.fy = ROOT.std.vector('double')()
		self.treconStruct.fyup = ROOT.std.vector('double')()
		self.treconStruct.fydown = ROOT.std.vector('double')()
		self.treconStruct.fz = ROOT.std.vector('double')()
		self.treconStruct.fzup = ROOT.std.vector('double')()
		self.treconStruct.fzdown = ROOT.std.vector('double')()
		self.treconStruct.cy = ROOT.std.vector('double')()
		self.treconStruct.cyup = ROOT.std.vector('double')()
		self.treconStruct.cydown = ROOT.std.vector('double')()
		self.treconStruct.cz = ROOT.std.vector('double')()
		self.treconStruct.czup = ROOT.std.vector('double')()
		self.treconStruct.czdown = ROOT.std.vector('double')()
		self.treconStruct.bpflag = ROOT.std.vector('bool')()
		self.treconStruct.charge = ROOT.std.vector('int')()
		self.treconStruct.angles = ROOT.std.vector(ROOT.std.vector('double'))()
                self.ti = ROOT.std.vector('double')()

		self.tree.Branch('evtId', ROOT.AddressOf(self.treconStruct, 'evtId'), 'evtId/I')
		self.tree.Branch('clId', ROOT.AddressOf(self.treconStruct, 'clId'), 'clId/I')
		self.tree.Branch('pathId', self.treconStruct.pathId)
		self.tree.Branch('fy', self.treconStruct.fy)
		self.tree.Branch('fyup', self.treconStruct.fyup)
		self.tree.Branch('fydown', self.treconStruct.fydown)
		self.tree.Branch('fz', self.treconStruct.fz)
		self.tree.Branch('fzup', self.treconStruct.fzup)
		self.tree.Branch('fzdown', self.treconStruct.fzdown)
		self.tree.Branch('cy', self.treconStruct.cy)
		self.tree.Branch('cyup', self.treconStruct.cyup)
		self.tree.Branch('cydown', self.treconStruct.cydown)
		self.tree.Branch('cz', self.treconStruct.cz)
		self.tree.Branch('czup', self.treconStruct.czup)
		self.tree.Branch('czdown', self.treconStruct.czdown)
		self.tree.Branch('bpflag', self.treconStruct.bpflag)
		self.tree.Branch('charge', self.treconStruct.charge)
		self.tree.Branch('angles', self.treconStruct.angles)
		
	def buildTreeHeader(self, origin, nevents):

            if self.builtTreeHeader is False:

		tlist = self.tree.GetUserInfo()
		map = ROOT.TMap() # holds all key,value pairs: only one here
		map.SetName("Event Data Header")

		k = ROOT.TObjString("raw")# where the raw data is
		v = ROOT.TObjString(origin)
		map.Add(k,v)
		k = ROOT.TObjString("trecon_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def StoreTValInfo(self, treconDict):
		for k,candidates in treconDict.iteritems():
			# Store trecon info in output tree
			self.treconStruct.evtId = self.event_id
			self.treconStruct.clId = k

			self.treconStruct.pathId.clear()
			self.treconStruct.fy.clear()
			self.treconStruct.fyup.clear()
			self.treconStruct.fydown.clear()
			self.treconStruct.fz.clear()
			self.treconStruct.fzup.clear()
			self.treconStruct.fzdown.clear()
			self.treconStruct.cy.clear()
			self.treconStruct.cyup.clear()
			self.treconStruct.cydown.clear()
			self.treconStruct.cz.clear()
			self.treconStruct.czup.clear()
			self.treconStruct.czdown.clear()
			self.treconStruct.bpflag.clear()
			self.treconStruct.angles.clear()
			self.treconStruct.charge.clear()
			for l,path in candidates.iteritems():
				self.treconStruct.pathId.push_back(l)
				rfoil = path[0] # tuple
				rcalo = path[1] # tuple
				bp = path[2] # bool
				phi = path[3] # list of double

				self.treconStruct.fy.push_back(rfoil[0])
				self.treconStruct.fyup.push_back(rfoil[1])
				self.treconStruct.fydown.push_back(rfoil[2])
				self.treconStruct.fz.push_back(rfoil[3])
				self.treconStruct.fzup.push_back(rfoil[4])
				self.treconStruct.fzdown.push_back(rfoil[5])
				self.treconStruct.cy.push_back(rcalo[0])
				self.treconStruct.cyup.push_back(rcalo[1])
				self.treconStruct.cydown.push_back(rcalo[2])
				self.treconStruct.cz.push_back(rcalo[3])
				self.treconStruct.czup.push_back(rcalo[4])
				self.treconStruct.czdown.push_back(rcalo[5])
				self.treconStruct.bpflag.push_back(bp)
				self.treconStruct.charge.push_back(0)
				
				self.ti.clear()
				for a in phi:
					self.ti.push_back(a)
				self.treconStruct.angles.push_back(self.ti)
				
			self.tree.Fill()
