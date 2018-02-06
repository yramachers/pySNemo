import ROOT
from pysnemo.io.edm import Particle
from pysnemo.utility.helix import HelixFit

ROOT.gROOT.ProcessLine(
"struct ParticleEventStorage{\
  int evtId;\
  int nofparticles;\
  int nofgammas;\
  std::vector<int>* cluster_id;\
  std::vector<int>* path_id;\
  std::vector<double>* charge;\
  std::vector<double>* chargeerr;\
  std::vector<double>* vertex_y;\
  std::vector<double>* vertex_z;\
  std::vector<double>* foil_dir_x;\
  std::vector<double>* foil_dir_y;\
  std::vector<double>* foil_dir_z;\
  std::vector<int>* calo_associated;\
  std::vector<int>* calo_type;\
  std::vector<double>* calo_energy;\
  std::vector<double>* calo_sigma_energy;\
  std::vector<double>* calo_time;\
  std::vector<double>* calo_sigma_time;\
  std::vector<int>* calo_side;\
  std::vector<int>* calo_column;\
  std::vector<int>* calo_row;\
  std::vector<int>* calo_wall;\
  std::vector<double>* calo_loc_x;\
  std::vector<double>* calo_loc_y;\
  std::vector<double>* calo_loc_z;\
};");

class ptdToRoot(object):
	'''
	Stream full reconstruction output stored in event to a Root file.
        Storage only, no reader required.
	keystring: key in event containing a list format 
	           hits structure (list of tuples)
	'''
	def __init__(self, keystring, alt_geom = False):
            self.key = keystring
            self.builtTreeHeader = False
            self.totNumEvents = 0
	    self.alt_geom = alt_geom
	    self.ptdStruct = ROOT.ParticleEventStorage()
	
	def __repr__(self):
            s = 'Full Reconstruction to Root streamer'
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
                if (not outputFile.Get("ptd_tree")):
                    self.tree = ROOT.TTree("ptd_tree", "Reconstruction data")
                    self.tree.SetDirectory(self.outputFile)
                else:
                    self.tree = outputFile.Get("ptd_tree")
                    self.tree.SetDirectory(self.outputFile)
                self.setupTree()

        
        def fillTree(self, event):

            self.outputFile.cd()
            if self.builtTreeHeader is False:
                origin = event.getFilename()            
                self.buildTreeHeader(origin, self.totNumEvents)

            self.event_id = event.getEventID()

            if (self.key is not None):
                if (event.hasKey(self.key)):
                    tvalDict = event.getKeyValue(self.key)
		    self.StoreData(tvalDict)
                else:
                    print "Error in pipeline: key (%s) not in event" % self.key


        def writeTree(self):

            self.outputFile.cd()
            self.tree.Write()


	def setupTree(self):

            self.ptdStruct.cluster_id = ROOT.std.vector('int')()
            self.ptdStruct.path_id = ROOT.std.vector('int')()
            self.ptdStruct.charge = ROOT.std.vector('double')()
            self.ptdStruct.chargeerr = ROOT.std.vector('double')()
            self.ptdStruct.vertex_y = ROOT.std.vector('double')()
            self.ptdStruct.vertex_z = ROOT.std.vector('double')()
            self.ptdStruct.foil_dir_x = ROOT.std.vector('double')()
            self.ptdStruct.foil_dir_y = ROOT.std.vector('double')()
            self.ptdStruct.foil_dir_z = ROOT.std.vector('double')()
            self.ptdStruct.calo_associated = ROOT.std.vector('int')()
            self.ptdStruct.calo_type = ROOT.std.vector('int')()
            self.ptdStruct.calo_energy = ROOT.std.vector('double')()
            self.ptdStruct.calo_sigma_energy = ROOT.std.vector('double')()
            self.ptdStruct.calo_time = ROOT.std.vector('double')()
            self.ptdStruct.calo_sigma_time = ROOT.std.vector('double')()
            self.ptdStruct.calo_side = ROOT.std.vector('int')()
            self.ptdStruct.calo_column = ROOT.std.vector('int')()
            self.ptdStruct.calo_row = ROOT.std.vector('int')()
            self.ptdStruct.calo_wall = ROOT.std.vector('int')()
            self.ptdStruct.calo_loc_x = ROOT.std.vector('double')()
            self.ptdStruct.calo_loc_y = ROOT.std.vector('double')()
            self.ptdStruct.calo_loc_z = ROOT.std.vector('double')()

            self.tree.Branch('eventId', ROOT.AddressOf(self.ptdStruct, 'evtId'), 'evtId/I')
            self.tree.Branch('nofparticles', ROOT.AddressOf(self.ptdStruct, 'nofparticles'), 'nofparticles/I')
            self.tree.Branch('nofgammas', ROOT.AddressOf(self.ptdStruct, 'nofgammas'), 'nofgammas/I')
            self.tree.Branch('cluster_id', self.ptdStruct.cluster_id)
            self.tree.Branch('path_id', self.ptdStruct.path_id)
            self.tree.Branch('charge', self.ptdStruct.charge)
            self.tree.Branch('charge_err', self.ptdStruct.chargeerr)
            self.tree.Branch('vertex_y', self.ptdStruct.vertex_y)
            self.tree.Branch('vertex_z', self.ptdStruct.vertex_z)
            self.tree.Branch('foil_dir_x', self.ptdStruct.foil_dir_x)
            self.tree.Branch('foil_dir_y', self.ptdStruct.foil_dir_y)
            self.tree.Branch('foil_dir_z', self.ptdStruct.foil_dir_z)
            self.tree.Branch('calo_associated', self.ptdStruct.calo_associated)
            self.tree.Branch('calo_type', self.ptdStruct.calo_type)
            self.tree.Branch('calo_energy', self.ptdStruct.calo_energy)
            self.tree.Branch('calo_sigma_energy', self.ptdStruct.calo_sigma_energy)
            self.tree.Branch('calo_time', self.ptdStruct.calo_time)
            self.tree.Branch('calo_sigma_time', self.ptdStruct.calo_sigma_time)
            self.tree.Branch('calo_side', self.ptdStruct.calo_side)
            self.tree.Branch('calo_column', self.ptdStruct.calo_column)
            self.tree.Branch('calo_row', self.ptdStruct.calo_row)
            self.tree.Branch('calo_wall', self.ptdStruct.calo_wall)
            self.tree.Branch('calo_loc_x', self.ptdStruct.calo_loc_x)
            self.tree.Branch('calo_loc_y', self.ptdStruct.calo_loc_y)
            self.tree.Branch('calo_loc_z', self.ptdStruct.calo_loc_z)

            
        def buildTreeHeader(self, origin, nevents):

            if self.builtTreeHeader is False:

		tlist = self.tree.GetUserInfo()
		map = ROOT.TMap() # holds all key,value pairs: only one here
		map.SetName("Event Data Header")

                if isinstance(origin,list):
                    for fn in origin:
                        k = ROOT.TObjString("data"+str(origin.index(fn)))# where the raw data is
                        v = ROOT.TObjString(fn)
                        map.Add(k,v)
                else:
                    k = ROOT.TObjString("raw")# where the raw data is
                    v = ROOT.TObjString(origin)
                    map.Add(k,v)
		k = ROOT.TObjString("ptd_tree") # what type of data
		v = ROOT.TObjString(self.key) # what was the processing key
		map.Add(k,v)

		nevts = ROOT.TParameter(int)("nevents",nevents) # storing
		# numerical parameter is easier like this

		tlist.AddFirst(map) # should be the first entry in TList
		tlist.Add(nevts) # appended to list
		# there is space for more objects.

                self.builtTreeHeader = True


	def StoreData(self, ptdData):
            if ptdData is None:
		    return
	    # count first
	    p = 0
	    g = 0
            for entry in ptdData:
                if isinstance(entry,Particle):
                    p += 1 # a particle
#		    print 'Particle: ',entry
                else:
#		    print 'Iso Gamma: ',entry
                    g += 1 # an isolated calo hit
                    
	    self.ptdStruct.cluster_id.clear()
	    self.ptdStruct.path_id.clear()
	    self.ptdStruct.charge.clear()
	    self.ptdStruct.chargeerr.clear()
	    self.ptdStruct.vertex_y.clear()
	    self.ptdStruct.vertex_z.clear()
	    self.ptdStruct.foil_dir_x.clear()
	    self.ptdStruct.foil_dir_y.clear()
	    self.ptdStruct.foil_dir_z.clear()
	    self.ptdStruct.calo_associated.clear()
	    self.ptdStruct.calo_type.clear()
	    self.ptdStruct.calo_energy.clear()
	    self.ptdStruct.calo_sigma_energy.clear()
	    self.ptdStruct.calo_time.clear()
	    self.ptdStruct.calo_sigma_time.clear()
	    self.ptdStruct.calo_side.clear()
	    self.ptdStruct.calo_column.clear()
	    self.ptdStruct.calo_row.clear()
	    self.ptdStruct.calo_wall.clear()
	    self.ptdStruct.calo_loc_x.clear()
	    self.ptdStruct.calo_loc_y.clear()
	    self.ptdStruct.calo_loc_z.clear()

	    self.ptdStruct.evtId = self.event_id
	    self.ptdStruct.nofparticles = p
	    self.ptdStruct.nofgammas = g
            for particle in ptdData:
                # Store ptd info in output tree

		if isinstance(particle,Particle): # found particles
			pid = particle.id
			self.ptdStruct.cluster_id.push_back(pid[0])
			self.ptdStruct.path_id.push_back(pid[1])

			fvtx = particle.foilvertex
			self.ptdStruct.vertex_y.push_back(fvtx[0])
			self.ptdStruct.vertex_z.push_back(fvtx[1])

			fitter = particle.fitter
			side = particle.id[1]
			if isinstance(fitter,HelixFit):
				if side < 0: # left tracker, invert fitted charge and momentum
					self.ptdStruct.charge.push_back(-fitter.fitcharge)
					self.ptdStruct.foil_dir_x.push_back(-fitter.fitmomentum[0])
					self.ptdStruct.foil_dir_y.push_back(-fitter.fitmomentum[1])
					self.ptdStruct.foil_dir_z.push_back(fitter.fitmomentum[2])
				else:
					self.ptdStruct.charge.push_back(fitter.fitcharge)
					self.ptdStruct.foil_dir_x.push_back(fitter.fitmomentum[0])
					self.ptdStruct.foil_dir_y.push_back(fitter.fitmomentum[1])
					self.ptdStruct.foil_dir_z.push_back(fitter.fitmomentum[2])
				self.ptdStruct.chargeerr.push_back(fitter.fitcharge_error)
			elif isinstance(fitter, list): # broken line solution
				self.ptdStruct.charge.push_back(0.0)
				self.ptdStruct.chargeerr.push_back(0.0)
				self.ptdStruct.foil_dir_x.push_back(1.0)
				vertexline = fitter[0]
				linepar = vertexline[0]
				if len(linepar)>2:
					self.ptdStruct.foil_dir_y.push_back(linepar[1])# slope xy
					self.ptdStruct.foil_dir_z.push_back(linepar[3])# slope xz
				else:
					self.ptdStruct.foil_dir_y.push_back(linepar[1])# slope xy
					self.ptdStruct.foil_dir_z.push_back(linepar[1])# slope xy
			else:
				self.ptdStruct.charge.push_back(0.0)
				self.ptdStruct.chargeerr.push_back(0.0)
				self.ptdStruct.foil_dir_x.push_back(1.0)
				linepar = fitter[0]
				if len(linepar)>2:
					self.ptdStruct.foil_dir_y.push_back(linepar[1])# slope xy
					self.ptdStruct.foil_dir_z.push_back(linepar[3])# slope xz
				else:
					self.ptdStruct.foil_dir_y.push_back(linepar[1])# slope xy
					self.ptdStruct.foil_dir_z.push_back(linepar[1])# slope xy

			chit = particle.calo_hit
			calomi = chit.meta_info
			type = calomi[1]
			side = calomi[3] # the tracker side
			wall = calomi[6] # the wall
			self.ptdStruct.calo_associated.push_back(1)
			self.ptdStruct.calo_type.push_back(type)
			self.ptdStruct.calo_energy.push_back(chit.energy)
			self.ptdStruct.calo_sigma_energy.push_back(chit.sigma_e)
			self.ptdStruct.calo_time.push_back(chit.time)
			self.ptdStruct.calo_sigma_time.push_back(chit.sigma_t)
			self.ptdStruct.calo_side.push_back(side)
			self.ptdStruct.calo_column.push_back(calomi[4])
			self.ptdStruct.calo_row.push_back(calomi[5])
			self.ptdStruct.calo_wall.push_back(wall)

			cvtx = particle.calovertex
			# Main Wall
			if type==0:
				# main calo sitting at x = +- 43.5cm, or +-27.5738cm
				if side == 0:
					if self.alt_geom:
						self.ptdStruct.calo_loc_x.push_back(-275.738)
					else:
						self.ptdStruct.calo_loc_x.push_back(-435.0)
					self.ptdStruct.calo_loc_y.push_back(cvtx[0])
					self.ptdStruct.calo_loc_z.push_back(cvtx[1])
				else:
					if self.alt_geom:
						self.ptdStruct.calo_loc_x.push_back(275.738)
					else:
						self.ptdStruct.calo_loc_x.push_back(435.0)
					self.ptdStruct.calo_loc_y.push_back(cvtx[0])
					self.ptdStruct.calo_loc_z.push_back(cvtx[1])
					
			# X Wall
			elif type==1:
				# xwall calo sitting at y = +- 250.55cm
				if wall == 0:
					self.ptdStruct.calo_loc_y.push_back(-2505.5)
					self.ptdStruct.calo_loc_x.push_back(cvtx[0])
					self.ptdStruct.calo_loc_z.push_back(cvtx[1])
				else:
					self.ptdStruct.calo_loc_y.push_back(2505.5)
					self.ptdStruct.calo_loc_x.push_back(cvtx[0])
					self.ptdStruct.calo_loc_z.push_back(cvtx[1])
						
			# Gamma Veto
			elif type==2:
				# gveto calo sitting at z = +- 155.0cm
				if wall == 0:
					self.ptdStruct.calo_loc_z.push_back(-1550.0)
					self.ptdStruct.calo_loc_x.push_back(cvtx[0])
					self.ptdStruct.calo_loc_y.push_back(cvtx[1])
				else:
					self.ptdStruct.calo_loc_z.push_back(1550.0)
					self.ptdStruct.calo_loc_x.push_back(cvtx[0])
					self.ptdStruct.calo_loc_y.push_back(cvtx[1])

		else: # non-associated calo hit
			calomi = particle.meta_info
			self.ptdStruct.cluster_id.push_back(calomi[0])
			self.ptdStruct.path_id.push_back(-1)
			self.ptdStruct.vertex_y.push_back(0.0)
			self.ptdStruct.vertex_z.push_back(0.0)
			self.ptdStruct.charge.push_back(0.0)
			self.ptdStruct.chargeerr.push_back(0.0)
			self.ptdStruct.foil_dir_x.push_back(0.0)
			self.ptdStruct.foil_dir_y.push_back(0.0)
			self.ptdStruct.foil_dir_z.push_back(0.0)
			
			type = calomi[1]
			side = calomi[3] # the tracker side
			wall = calomi[6] # the wall
			self.ptdStruct.calo_associated.push_back(0)

			self.ptdStruct.calo_type.push_back(type)
			self.ptdStruct.calo_energy.push_back(particle.energy)
			self.ptdStruct.calo_sigma_energy.push_back(particle.sigma_e)
			self.ptdStruct.calo_time.push_back(particle.time)
			self.ptdStruct.calo_sigma_time.push_back(particle.sigma_t)
			self.ptdStruct.calo_side.push_back(side)
			self.ptdStruct.calo_column.push_back(calomi[4])
			self.ptdStruct.calo_row.push_back(calomi[5])
			self.ptdStruct.calo_wall.push_back(wall)
			self.ptdStruct.calo_loc_y.push_back(0.0)
			self.ptdStruct.calo_loc_x.push_back(0.0)
			self.ptdStruct.calo_loc_z.push_back(0.0)
			
	    self.tree.Fill()
                        
