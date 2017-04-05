# pySNemoReaderDecorator -  reads data from disk and provides
#                           functionality to present an event in a cleaner
#                           dictionary based way.
#
# Copyright (c) 2013 by YR
# Copyright (c) 2011 by The University of Warwick
#
# pySNemo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pySNemo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with pySNemo.  If not, see <http://www.gnu.org/licenses/>.

import ROOT

ROOT.gROOT.ProcessLine(
"struct ToySimDataStruct{\
   vector<double>* dirx;\
   vector<double>* diry;\
   vector<double>* dirz;\
   vector<double>* pointx;\
   vector<double>* pointy;\
   vector<double>* pointz;\
   vector<double>* radius;\
   vector<double>* wirex;\
   vector<double>* wirey;\
   vector<double>* wirez;\
   vector<int>*    grid_id;\
   vector<int>*    grid_side;\
   vector<int>*    grid_layer;\
   vector<int>*    grid_column;\
};");

ROOT.gROOT.ProcessLine(
"struct flsimDataStruct{\
   int             truetracker_nohits;\
   vector<int>*    truetracker_id;\
   vector<int>*    truetracker_module;\
   vector<int>*    truetracker_side;\
   vector<int>*    truetracker_layer;\
   vector<int>*    truetracker_column;\
   vector<double>* truetracker_time;\
   vector<double>* truetracker_xstart;\
   vector<double>* truetracker_ystart;\
   vector<double>* truetracker_zstart;\
   vector<double>* truetracker_xstop;\
   vector<double>* truetracker_ystop;\
   vector<double>* truetracker_zstop;\
   vector<int>*    truetracker_trackid;\
   vector<int>*    truetracker_parenttrackid;\
   int             tracker_nohits;\
   vector<int>*    tracker_id;\
   vector<int>*    tracker_module;\
   vector<int>*    tracker_side;\
   vector<int>*    tracker_layer;\
   vector<int>*    tracker_column;\
   vector<double>* tracker_x;\
   vector<double>* tracker_y;\
   vector<double>* tracker_z;\
   vector<double>* tracker_sigmaz;\
   vector<double>* tracker_r;\
   vector<double>* tracker_sigmar;\
   vector<int>*    tracker_truehitid;\
   int             truecalo_nohits;\
   vector<int>*    truecalo_id;\
   vector<int>*    truecalo_type;\
   vector<int>*    truecalo_module;\
   vector<int>*    truecalo_side;\
   vector<int>*    truecalo_column;\
   vector<int>*    truecalo_row;\
   vector<int>*    truecalo_wall;\
   vector<double>* truecalo_time;\
   vector<double>* truecalo_x;\
   vector<double>* truecalo_y;\
   vector<double>* truecalo_z;\
   vector<double>* truecalo_energy;\
   int             calo_nohits;\
   vector<int>*    calo_id;\
   vector<int>*    calo_type;\
   vector<int>*    calo_module;\
   vector<int>*    calo_side;\
   vector<int>*    calo_column;\
   vector<int>*    calo_row;\
   vector<int>*    calo_wall;\
   vector<double>* calo_time;\
   vector<double>* calo_sigmatime;\
   vector<double>* calo_energy;\
   vector<double>* calo_sigmaenergy;\
};");

ROOT.gROOT.ProcessLine(
"struct TruthDataStruct{\
   double          truevertex_x;\
   double          truevertex_y;\
   double          truevertex_z;\
   double          truevertex_time;\
   int             trueparticle_noparticles;\
   vector<int>*    trueparticle_id;\
   vector<int>*    trueparticle_type;\
   vector<double>* trueparticle_kinenergy;\
   vector<double>* trueparticle_px;\
   vector<double>* trueparticle_py;\
   vector<double>* trueparticle_pz;\
   vector<double>* trueparticle_time;\
};");


class flsimreader(object):
    '''
    Simple file reader of flsimulate output root files with 
    selected data fields as specified in the flsimDataStruct.

    Takes root file pointer (a TFile object) and returns immediately
    the tree (hard-coded name: "SimData") and a flsimDataStruct object.
    Filling the data structure takes place in the decorator when 
    constructing an Event.
    '''
    def __init__(self, file):
        self.tree = ROOT.TTree()
        self.ds0 = ROOT.flsimDataStruct()
        self.ds1 = ROOT.TruthDataStruct()
        self.tree = file.Get("SimData")
        self.tree.SetDirectory(file)


class tsimreader(object):
    '''
    Simple file reader of toysimulation output root files  
    as specified in the ToySimDataStruct.

    Takes root file pointer (a TFile object) and returns immediately
    the tree (hard-coded name: "hit_tree") and a ToySimDataStruct object.
    Filling the data structure takes place in the decorator when 
    constructing an Event.
    '''
    def __init__(self, file):
        self.ds = ROOT.ToySimDataStruct()
        self.tree = file.Get("hit_tree")


