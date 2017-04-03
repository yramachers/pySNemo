import ROOT
from numpy import array


ROOT.gROOT.ProcessLine(
"struct tvalidationStruct{\
   int evtId;\
   int nClusters;\
   vector<int>* mainId;\
   vector<double>* fresy;\
   vector<double>* fresyup;\
   vector<double>* fresydown;\
   vector<double>* fresz;\
   vector<double>* freszup;\
   vector<double>* freszdown;\
   vector<double>* cresy;\
   vector<double>* cresyup;\
   vector<double>* cresydown;\
   vector<double>* cresz;\
   vector<double>* creszup;\
   vector<double>* creszdown;\
};");

class trackvalidationToRoot(object):
	'''
	Stream track validation output stored in event to a Root file.
	keystring: key in event containing a list format 
	           hits structure (list of tuples)
	'''
	def __init__(self, keystring):
            self.key = keystring
            self.builtTreeHeader = False
            self.totNumEvents = 0
	    self.tvalidationStruct = ROOT.tvalidationStruct()
	
	def __repr__(self):
            s = 'Track Validation to Root streamer'
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
                if (not outputFile.Get("tvalidation_tree")):
                    self.tree = ROOT.TTree("tvalidation_tree", "TValidation data")
                    self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get("tvalidation_tree")
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

		self.tvalidationStruct.mainId = ROOT.std.vector('int')()
		self.tvalidationStruct.fresy = ROOT.std.vector('double')()
		self.tvalidationStruct.fresyup = ROOT.std.vector('double')()
		self.tvalidationStruct.fresydown = ROOT.std.vector('double')()
		self.tvalidationStruct.fresz = ROOT.std.vector('double')()
		self.tvalidationStruct.freszup = ROOT.std.vector('double')()
		self.tvalidationStruct.freszdown = ROOT.std.vector('double')()
		self.tvalidationStruct.cresy = ROOT.std.vector('double')()
		self.tvalidationStruct.cresyup = ROOT.std.vector('double')()
		self.tvalidationStruct.cresydown = ROOT.std.vector('double')()
		self.tvalidationStruct.cresz = ROOT.std.vector('double')()
		self.tvalidationStruct.creszup = ROOT.std.vector('double')()
		self.tvalidationStruct.creszdown = ROOT.std.vector('double')()

		self.tree.Branch('evtId', ROOT.AddressOf(self.tvalidationStruct, 'evtId'), 'evtId/I')
		self.tree.Branch('nClusters', ROOT.AddressOf(self.tvalidationStruct, 'nClusters'), 'nClusters/I')
		self.tree.Branch('mainId', self.tvalidationStruct.mainId)
		self.tree.Branch('fresy', self.tvalidationStruct.fresy)
		self.tree.Branch('fresyup', self.tvalidationStruct.fresyup)
		self.tree.Branch('fresydown', self.tvalidationStruct.fresydown)
		self.tree.Branch('fresz', self.tvalidationStruct.fresz)
		self.tree.Branch('freszup', self.tvalidationStruct.freszup)
		self.tree.Branch('freszdown', self.tvalidationStruct.freszdown)
		self.tree.Branch('cresy', self.tvalidationStruct.cresy)
		self.tree.Branch('cresyup', self.tvalidationStruct.cresyup)
		self.tree.Branch('cresydown', self.tvalidationStruct.cresydown)
		self.tree.Branch('cresz', self.tvalidationStruct.cresz)
		self.tree.Branch('creszup', self.tvalidationStruct.creszup)
		self.tree.Branch('creszdown', self.tvalidationStruct.creszdown)
		
	def buildTreeHeader(self, origin, nevents):

            if self.builtTreeHeader is False:

		tlist = self.tree.GetUserInfo()
		map = ROOT.TMap() # holds all key,value pairs: only one here
		map.SetName("Event Data Header")

		k = ROOT.TObjString("raw")# where the raw data is
		v = ROOT.TObjString(origin)
		map.Add(k,v)
		k = ROOT.TObjString("tvalidation_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def StoreTValInfo(self, tvalidationDict):

		# Store tvalidation info in output tree
		self.tvalidationStruct.evtId = self.event_id
		self.tvalidationStruct.nClusters = len(tvalidationDict)
			
		self.tvalidationStruct.mainId.clear()
		self.tvalidationStruct.fresy.clear()
		self.tvalidationStruct.fresyup.clear()
		self.tvalidationStruct.fresydown.clear()
		self.tvalidationStruct.fresz.clear()
		self.tvalidationStruct.freszup.clear()
		self.tvalidationStruct.freszdown.clear()
		self.tvalidationStruct.cresy.clear()
		self.tvalidationStruct.cresyup.clear()
		self.tvalidationStruct.cresydown.clear()
		self.tvalidationStruct.cresz.clear()
		self.tvalidationStruct.creszup.clear()
		self.tvalidationStruct.creszdown.clear()

		for k,val in tvalidationDict.iteritems():
			rfoil = val[0]
			rcalo = val[1]

			self.tvalidationStruct.mainId.push_back(k)
			self.tvalidationStruct.fresy.push_back(rfoil[0])
			self.tvalidationStruct.fresyup.push_back(rfoil[1])
			self.tvalidationStruct.fresydown.push_back(rfoil[2])
			self.tvalidationStruct.fresz.push_back(rfoil[3])
			self.tvalidationStruct.freszup.push_back(rfoil[4])
			self.tvalidationStruct.freszdown.push_back(rfoil[5])
			self.tvalidationStruct.cresy.push_back(rcalo[0])
			self.tvalidationStruct.cresyup.push_back(rcalo[1])
			self.tvalidationStruct.cresydown.push_back(rcalo[2])
			self.tvalidationStruct.cresz.push_back(rcalo[3])
			self.tvalidationStruct.creszup.push_back(rcalo[4])
			self.tvalidationStruct.creszdown.push_back(rcalo[5])

		self.tree.Fill()

