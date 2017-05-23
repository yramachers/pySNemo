# pySNemoReaderDecorator -  reads data from disk and provides
#                           functionality to present an event in a cleaner
#                           dictionary based way.
#
# Copyright (c) 2013 by YR
# Copyright (c) 2011 by The University of Warwick
#
# This file is part of pySNemo.
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

import re
import math
import pysnemo.io.edm as EDM
from pysnemo.io.io_formats import flsimreader
from pysnemo.io.io_formats import tsimreader

from pysnemo.control.output.tr_clusterToRoot import clusterToRoot
from pysnemo.control.output.tr_pathToRoot import pathToRoot
from pysnemo.control.output.tr_listToRoot import listToRoot

import ROOT

"""
flsimulate file reader decorators
---------------------------------

The flsim FileReader object provided by the io_formats package simply
provides raw access to the data collections in an flsim root output file.

This module provides some simple decorator objects to present this data in
more usable forms.

Toysimulator file reader decorators
-----------------------------------
The ToySim FileReader object provided by the io_formats package simply
provides raw access to the data collections in a toysimulation root 
output file.

This module provides some simple decorator objects to present this data in
more usable forms.

"""


class flsimRun(object):
    """basic flsimFileReader decorator
    Must be a ROOT file containing TTree.

    It constructs a dictionary based event from the data so that later 
    algorithms may easily extract needed raw and truth data, as well as add 
    their own data to the event. The events can thus be used in a pipeline
    type analysis

    It converts the entries in a flsimulate root
    output file to an event structure:

    event -> dict:
        key : "raw"
            key : "gg" -> tracker hits
            key : "ggtruth" -> truth info on tracker hits from simulation
            key : "calo_hits" -> calorimeter hits
            key : "calotruth" -> calo truth hits
            key : "vertex" -> vertex info
            key : "true_particle" -> truth info about simulated particle

    Individual events can be extracted via the get_event method, and all events
    can be iterated over via the events generator method.
    """
    def __init__(self, filename):
        """construct an instance to read file filename.
        """
        self._file = ROOT.TFile(filename,"open")
        self._file.cd()
        # get tree and wrap content structure in datastruct for filling
        r  = flsimreader(self._file)
        self.tree = r.tree
        self.datastruct = r.ds0
        self.datastruct_truth = r.ds1


    def number_of_events(self):
        """return number of events in the file
        """
        return self.tree.GetEntries()


    def get_event(self, index):
        """return the event at position index
        Functions build selected raw data from 
        self.datastruct into the event dictionary.
        """
        event = {}
        self._setup_structure()
        self.tree.GetEntry(index) # now self.datastruct is filled

        event["raw"] = {}
        if (self.datastruct.tracker_nohits>0) :
            event["raw"]["gg"] = self.__build_observed_gg()

        if (self.datastruct.truetracker_nohits>0) :
            event["raw"]["ggtruth"] = self.__build_truth_gg()

        if (self.datastruct.calo_nohits>0) :
            event["raw"]["calo_hits"] = self.__build_calo_hits()

        if (self.datastruct.truecalo_nohits>0) :
            event["raw"]["calotruth"] = self.__build_calo_truth()

        event["raw"]["vertex"] = self.__build_vertex()
        event["raw"]["trueparticle"] = self.__build_trueparticle()

        return event


    def events(self):
        """generator yielding all events in the file in sequence
        """
        for i in range(self.number_of_events()):
            yield self.get_event(i)


    def _setup_structure(self):
        '''
        set up structure of flsim output
        '''
        self.datastruct.truetracker_nohits = 0;
        self.datastruct.truetracker_id = ROOT.std.vector('int')()
        self.datastruct.truetracker_module = ROOT.std.vector('int')()
        self.datastruct.truetracker_side = ROOT.std.vector('int')()
        self.datastruct.truetracker_layer = ROOT.std.vector('int')()
        self.datastruct.truetracker_column = ROOT.std.vector('int')()
        self.datastruct.truetracker_trackid = ROOT.std.vector('int')()
        self.datastruct.truetracker_parenttrackid = ROOT.std.vector('int')()
        self.datastruct.truetracker_time = ROOT.std.vector('double')()
        self.datastruct.truetracker_xstart = ROOT.std.vector('double')()
        self.datastruct.truetracker_ystart = ROOT.std.vector('double')()
        self.datastruct.truetracker_zstart = ROOT.std.vector('double')()
        self.datastruct.truetracker_xstop = ROOT.std.vector('double')()
        self.datastruct.truetracker_ystop = ROOT.std.vector('double')()
        self.datastruct.truetracker_zstop = ROOT.std.vector('double')()

        self.datastruct.truecalo_nohits = 0;
        self.datastruct.truecalo_id = ROOT.std.vector('int')()
        self.datastruct.truecalo_type = ROOT.std.vector('int')()
        self.datastruct.truecalo_module = ROOT.std.vector('int')()
        self.datastruct.truecalo_side = ROOT.std.vector('int')()
        self.datastruct.truecalo_column = ROOT.std.vector('int')()
        self.datastruct.truecalo_row = ROOT.std.vector('int')()
        self.datastruct.truecalo_wall = ROOT.std.vector('int')()
        self.datastruct.truecalo_time = ROOT.std.vector('double')()
        self.datastruct.truecalo_x = ROOT.std.vector('double')()
        self.datastruct.truecalo_y = ROOT.std.vector('double')()
        self.datastruct.truecalo_z = ROOT.std.vector('double')()
        self.datastruct.truecalo_energy = ROOT.std.vector('double')()
        
        self.datastruct.tracker_nohits = 0;
        self.datastruct.tracker_id = ROOT.std.vector('int')()
        self.datastruct.tracker_truehitid = ROOT.std.vector('int')()
        self.datastruct.tracker_module = ROOT.std.vector('int')()
        self.datastruct.tracker_side = ROOT.std.vector('int')()
        self.datastruct.tracker_layer = ROOT.std.vector('int')()
        self.datastruct.tracker_column = ROOT.std.vector('int')()
        self.datastruct.tracker_x = ROOT.std.vector('double')()
        self.datastruct.tracker_y = ROOT.std.vector('double')()
        self.datastruct.tracker_z = ROOT.std.vector('double')()
        self.datastruct.tracker_sigmaz = ROOT.std.vector('double')()
        self.datastruct.tracker_r = ROOT.std.vector('double')()
        self.datastruct.tracker_sigmar = ROOT.std.vector('double')()

        self.datastruct.calo_nohits = 0;
        self.datastruct.calo_id = ROOT.std.vector('int')()
        self.datastruct.calo_type = ROOT.std.vector('int')()
        self.datastruct.calo_module = ROOT.std.vector('int')()
        self.datastruct.calo_side = ROOT.std.vector('int')()
        self.datastruct.calo_column = ROOT.std.vector('int')()
        self.datastruct.calo_row = ROOT.std.vector('int')()
        self.datastruct.calo_wall = ROOT.std.vector('int')()
        self.datastruct.calo_time = ROOT.std.vector('double')()
        self.datastruct.calo_sigmatime = ROOT.std.vector('double')()
        self.datastruct.calo_energy = ROOT.std.vector('double')()
        self.datastruct.calo_sigmaenergy = ROOT.std.vector('double')()

        self.datastruct_truth.trueparticle_noparticles = 0;
        self.datastruct_truth.trueparticle_id = ROOT.std.vector('int')()
        self.datastruct_truth.trueparticle_type = ROOT.std.vector('int')()
        self.datastruct_truth.trueparticle_px = ROOT.std.vector('double')()
        self.datastruct_truth.trueparticle_py = ROOT.std.vector('double')()
        self.datastruct_truth.trueparticle_pz = ROOT.std.vector('double')()
        self.datastruct_truth.trueparticle_time = ROOT.std.vector('double')()
        self.datastruct_truth.trueparticle_kinenergy = ROOT.std.vector('double')()

        self.datastruct_truth.truevertex_x = 0.0;
        self.datastruct_truth.truevertex_y = 0.0;
        self.datastruct_truth.truevertex_z = 0.0;
        self.datastruct_truth.truevertex_time = 0.0;
        
        self.tree.SetBranchAddress('truetracker.nohits', ROOT.AddressOf(self.datastruct, "truetracker_nohits") )
        self.tree.SetBranchAddress('truetracker.id', ROOT.AddressOf(self.datastruct, "truetracker_id") )
        self.tree.SetBranchAddress('truetracker.module', ROOT.AddressOf(self.datastruct, "truetracker_module") )
        self.tree.SetBranchAddress('truetracker.side', ROOT.AddressOf(self.datastruct, "truetracker_side") )
        self.tree.SetBranchAddress('truetracker.layer', ROOT.AddressOf(self.datastruct, "truetracker_layer") )
        self.tree.SetBranchAddress('truetracker.column', ROOT.AddressOf(self.datastruct, "truetracker_column") )
        self.tree.SetBranchAddress('truetracker.time', ROOT.AddressOf(self.datastruct, "truetracker_time") )
        self.tree.SetBranchAddress('truetracker.xstart', ROOT.AddressOf(self.datastruct, "truetracker_xstart") )
        self.tree.SetBranchAddress('truetracker.ystart', ROOT.AddressOf(self.datastruct, "truetracker_ystart") )
        self.tree.SetBranchAddress('truetracker.zstart', ROOT.AddressOf(self.datastruct, "truetracker_zstart") )
        self.tree.SetBranchAddress('truetracker.xstop', ROOT.AddressOf(self.datastruct, "truetracker_xstop") )
        self.tree.SetBranchAddress('truetracker.ystop', ROOT.AddressOf(self.datastruct, "truetracker_ystop") )
        self.tree.SetBranchAddress('truetracker.zstop', ROOT.AddressOf(self.datastruct, "truetracker_zstop") )
        self.tree.SetBranchAddress('truetracker.trackid', ROOT.AddressOf(self.datastruct, "truetracker_trackid") )
        self.tree.SetBranchAddress('truetracker.parenttrackid', ROOT.AddressOf(self.datastruct, "truetracker_parenttrackid") )

        self.tree.SetBranchAddress('tracker.nohits', ROOT.AddressOf(self.datastruct, "tracker_nohits") )
        self.tree.SetBranchAddress('tracker.id', ROOT.AddressOf(self.datastruct, "tracker_id") )
        self.tree.SetBranchAddress('tracker.truehitid', ROOT.AddressOf(self.datastruct, "tracker_truehitid") )
        self.tree.SetBranchAddress('tracker.module', ROOT.AddressOf(self.datastruct, "tracker_module") )
        self.tree.SetBranchAddress('tracker.side', ROOT.AddressOf(self.datastruct, "tracker_side") )
        self.tree.SetBranchAddress('tracker.layer', ROOT.AddressOf(self.datastruct, "tracker_layer") )
        self.tree.SetBranchAddress('tracker.column', ROOT.AddressOf(self.datastruct, "tracker_column") )
        self.tree.SetBranchAddress('tracker.x', ROOT.AddressOf(self.datastruct, "tracker_x") )
        self.tree.SetBranchAddress('tracker.y', ROOT.AddressOf(self.datastruct, "tracker_y") )
        self.tree.SetBranchAddress('tracker.z', ROOT.AddressOf(self.datastruct, "tracker_z") )
        self.tree.SetBranchAddress('tracker.sigmaz', ROOT.AddressOf(self.datastruct, "tracker_sigmaz") )
        self.tree.SetBranchAddress('tracker.r', ROOT.AddressOf(self.datastruct, "tracker_r") )
        self.tree.SetBranchAddress('tracker.sigmar', ROOT.AddressOf(self.datastruct, "tracker_sigmar") )

        self.tree.SetBranchAddress('truecalo.nohits', ROOT.AddressOf(self.datastruct, "truecalo_nohits") )
        self.tree.SetBranchAddress('truecalo.id', ROOT.AddressOf(self.datastruct, "truecalo_id") )
        self.tree.SetBranchAddress('truecalo.type', ROOT.AddressOf(self.datastruct, "truecalo_type") )
        self.tree.SetBranchAddress('truecalo.module', ROOT.AddressOf(self.datastruct, "truecalo_module") )
        self.tree.SetBranchAddress('truecalo.side', ROOT.AddressOf(self.datastruct, "truecalo_side") )
        self.tree.SetBranchAddress('truecalo.column', ROOT.AddressOf(self.datastruct, "truecalo_column") )
        self.tree.SetBranchAddress('truecalo.row', ROOT.AddressOf(self.datastruct, "truecalo_row") )
        self.tree.SetBranchAddress('truecalo.wall', ROOT.AddressOf(self.datastruct, "truecalo_wall") )
        self.tree.SetBranchAddress('truecalo.time', ROOT.AddressOf(self.datastruct, "truecalo_time") )
        self.tree.SetBranchAddress('truecalo.x', ROOT.AddressOf(self.datastruct, "truecalo_x") )
        self.tree.SetBranchAddress('truecalo.y', ROOT.AddressOf(self.datastruct, "truecalo_y") )
        self.tree.SetBranchAddress('truecalo.z', ROOT.AddressOf(self.datastruct, "truecalo_z") )
        self.tree.SetBranchAddress('truecalo.energy', ROOT.AddressOf(self.datastruct, "truecalo_energy") )

        self.tree.SetBranchAddress('calo.nohits', ROOT.AddressOf(self.datastruct, "calo_nohits") )
        self.tree.SetBranchAddress('calo.id', ROOT.AddressOf(self.datastruct, "calo_id") )
        self.tree.SetBranchAddress('calo.type', ROOT.AddressOf(self.datastruct, "calo_type") )
        self.tree.SetBranchAddress('calo.module', ROOT.AddressOf(self.datastruct, "calo_module") )
        self.tree.SetBranchAddress('calo.side', ROOT.AddressOf(self.datastruct, "calo_side") )
        self.tree.SetBranchAddress('calo.column', ROOT.AddressOf(self.datastruct, "calo_column") )
        self.tree.SetBranchAddress('calo.row', ROOT.AddressOf(self.datastruct, "calo_row") )
        self.tree.SetBranchAddress('calo.wall', ROOT.AddressOf(self.datastruct, "calo_wall") )
        self.tree.SetBranchAddress('calo.time', ROOT.AddressOf(self.datastruct, "calo_time") )
        self.tree.SetBranchAddress('calo.sigmatime', ROOT.AddressOf(self.datastruct, "calo_sigmatime") )
        self.tree.SetBranchAddress('calo.energy', ROOT.AddressOf(self.datastruct, "calo_energy") )
        self.tree.SetBranchAddress('calo.sigmaenergy', ROOT.AddressOf(self.datastruct, "calo_sigmaenergy") )

        self.tree.SetBranchAddress('truevertex.x', ROOT.AddressOf(self.datastruct_truth, "truevertex_x") )
        self.tree.SetBranchAddress('truevertex.y', ROOT.AddressOf(self.datastruct_truth, "truevertex_y") )
        self.tree.SetBranchAddress('truevertex.z', ROOT.AddressOf(self.datastruct_truth, "truevertex_z") )
        self.tree.SetBranchAddress('truevertex.time', ROOT.AddressOf(self.datastruct_truth, "truevertex_time") )

        self.tree.SetBranchAddress('trueparticle.noparticles', ROOT.AddressOf(self.datastruct_truth, "trueparticle_noparticles") )
        self.tree.SetBranchAddress('trueparticle.id', ROOT.AddressOf(self.datastruct_truth, "trueparticle_id") )
        self.tree.SetBranchAddress('trueparticle.type', ROOT.AddressOf(self.datastruct_truth, "trueparticle_type") )
        self.tree.SetBranchAddress('trueparticle.px', ROOT.AddressOf(self.datastruct_truth, "trueparticle_px") )
        self.tree.SetBranchAddress('trueparticle.py', ROOT.AddressOf(self.datastruct_truth, "trueparticle_py") )
        self.tree.SetBranchAddress('trueparticle.pz', ROOT.AddressOf(self.datastruct_truth, "trueparticle_pz") )
        self.tree.SetBranchAddress('trueparticle.time', ROOT.AddressOf(self.datastruct_truth, "trueparticle_time") )
        self.tree.SetBranchAddress('trueparticle.kinenergy', ROOT.AddressOf(self.datastruct_truth, "trueparticle_kinenergy") )

        self.tree.SetBranchStatus("*", 0)
        # activate only the wanted branches 
        self.tree.SetBranchStatus("truetracker.nohits", 1) 
        self.tree.SetBranchStatus("truetracker.id", 1) 
        self.tree.SetBranchStatus("truetracker.module", 1) 
        self.tree.SetBranchStatus("truetracker.side", 1) 
        self.tree.SetBranchStatus("truetracker.layer", 1) 
        self.tree.SetBranchStatus("truetracker.column", 1) 
        self.tree.SetBranchStatus("truetracker.time", 1) 
        self.tree.SetBranchStatus("truetracker.xstart", 1) 
        self.tree.SetBranchStatus("truetracker.ystart", 1) 
        self.tree.SetBranchStatus("truetracker.zstart", 1) 
        self.tree.SetBranchStatus("truetracker.xstop", 1) 
        self.tree.SetBranchStatus("truetracker.ystop", 1) 
        self.tree.SetBranchStatus("truetracker.zstop", 1) 
        self.tree.SetBranchStatus("truetracker.trackid", 1) 
        self.tree.SetBranchStatus("truetracker.parenttrackid", 1) 
        self.tree.SetBranchStatus("tracker.nohits", 1) 
        self.tree.SetBranchStatus("tracker.id", 1) 
        self.tree.SetBranchStatus("tracker.module", 1) 
        self.tree.SetBranchStatus("tracker.side", 1) 
        self.tree.SetBranchStatus("tracker.layer", 1) 
        self.tree.SetBranchStatus("tracker.column", 1) 
        self.tree.SetBranchStatus("tracker.x", 1) 
        self.tree.SetBranchStatus("tracker.y", 1) 
        self.tree.SetBranchStatus("tracker.z", 1) 
        self.tree.SetBranchStatus("tracker.sigmaz", 1) 
        self.tree.SetBranchStatus("tracker.r", 1) 
        self.tree.SetBranchStatus("tracker.sigmar", 1) 
        self.tree.SetBranchStatus("tracker.truehitid", 1) 
        self.tree.SetBranchStatus("truecalo.nohits", 1) 
        self.tree.SetBranchStatus("truecalo.id", 1) 
        self.tree.SetBranchStatus("truecalo.type", 1) 
        self.tree.SetBranchStatus("truecalo.x", 1) 
        self.tree.SetBranchStatus("truecalo.y", 1) 
        self.tree.SetBranchStatus("truecalo.z", 1) 
        self.tree.SetBranchStatus("truecalo.time", 1) 
        self.tree.SetBranchStatus("truecalo.energy", 1) 
        self.tree.SetBranchStatus("truecalo.module", 1) 
        self.tree.SetBranchStatus("truecalo.side", 1) 
        self.tree.SetBranchStatus("truecalo.wall", 1) 
        self.tree.SetBranchStatus("truecalo.column", 1) 
        self.tree.SetBranchStatus("truecalo.row", 1) 
        self.tree.SetBranchStatus("calo.nohits", 1) 
        self.tree.SetBranchStatus("calo.id", 1) 
        self.tree.SetBranchStatus("calo.module", 1) 
        self.tree.SetBranchStatus("calo.side", 1) 
        self.tree.SetBranchStatus("calo.column", 1) 
        self.tree.SetBranchStatus("calo.row", 1) 
        self.tree.SetBranchStatus("calo.wall", 1) 
        self.tree.SetBranchStatus("calo.time", 1) 
        self.tree.SetBranchStatus("calo.sigmatime", 1) 
        self.tree.SetBranchStatus("calo.energy", 1) 
        self.tree.SetBranchStatus("calo.sigmaenergy", 1) 
        self.tree.SetBranchStatus("calo.type", 1) 
        self.tree.SetBranchStatus("truevertex.x", 1) 
        self.tree.SetBranchStatus("truevertex.y", 1) 
        self.tree.SetBranchStatus("truevertex.z", 1) 
        self.tree.SetBranchStatus("truevertex.time", 1) 
        self.tree.SetBranchStatus("trueparticle.noparticles", 1) 
        self.tree.SetBranchStatus("trueparticle.id", 1) 
        self.tree.SetBranchStatus("trueparticle.type", 1) 
        self.tree.SetBranchStatus("trueparticle.px", 1) 
        self.tree.SetBranchStatus("trueparticle.py", 1) 
        self.tree.SetBranchStatus("trueparticle.pz", 1) 
        self.tree.SetBranchStatus("trueparticle.time", 1) 
        self.tree.SetBranchStatus("trueparticle.kinenergy", 1) 
        

    def __build_observed_gg(self):
        """return list of geiger cylinders.
        """
        hits = []
        for i in range(self.datastruct.tracker_nohits):
            x = self.datastruct.tracker_x[i]
            y = self.datastruct.tracker_y[i]
            z = self.datastruct.tracker_z[i]
            dz = self.datastruct.tracker_sigmaz[i]
            r = self.datastruct.tracker_r[i]
            dr = self.datastruct.tracker_sigmar[i]
            if math.isnan(r):
                r = 23.0 # hard coded radius, dr for 
                dr = 1.0 # peripheral hits - extra large radius
            gg = EDM.tracker_hit(x,y,z,dz,r,dr)

            id = self.datastruct.tracker_id[i]
            tid = self.datastruct.tracker_truehitid[i]
            m = self.datastruct.tracker_module[i]
            s = self.datastruct.tracker_side[i]
            l = self.datastruct.tracker_layer[i]
            c = self.datastruct.tracker_column[i]
            gg.set_info(id,tid,m,s,l,c)
            hits.append(gg)
            
        return hits


    def __build_truth_gg(self):
        """return list truth info for geiger cylinders.
        """
        hits = []
        for i in range(self.datastruct.truetracker_time.size()):
            t  = self.datastruct.truetracker_time[i]
            x0 = self.datastruct.truetracker_xstart[i]
            y0 = self.datastruct.truetracker_ystart[i]
            z0 = self.datastruct.truetracker_zstart[i]
            x1 = self.datastruct.truetracker_xstop[i]
            y1 = self.datastruct.truetracker_ystop[i]
            z1 = self.datastruct.truetracker_zstop[i]
            gg = EDM.gcylinder_truth(t,x0,y0,z0,x1,y1,z1)

            id = self.datastruct.truetracker_id[i]
            m = self.datastruct.truetracker_module[i]
            s = self.datastruct.truetracker_side[i]
            l = self.datastruct.truetracker_layer[i]
            c = self.datastruct.truetracker_column[i]
            trid = self.datastruct.truetracker_trackid[i]
            pid = self.datastruct.truetracker_parenttrackid[i]
            gg.set_info(id,m,s,l,c,trid,pid)
            hits.append(gg)
            
        return hits


    def __build_calo_hits(self):
        """return list of calorimeter hit objects
        """
        data = []

        for i in range(self.datastruct.calo_nohits):
            t = self.datastruct.calo_time[i]
            dt = self.datastruct.calo_sigmatime[i]
            e = self.datastruct.calo_energy[i]
            de = self.datastruct.calo_sigmaenergy[i]
            ch = EDM.calo_hit(t,dt,e,de)

            id = self.datastruct.calo_id[i]
            type = self.datastruct.calo_type[i]
            m = self.datastruct.calo_module[i]
            s = self.datastruct.calo_side[i]
            c = self.datastruct.calo_column[i]
            r = self.datastruct.calo_row[i]
            w = self.datastruct.calo_wall[i]
            ch.set_info(id,type,m,s,c,r,w)
            data.append(ch)

        return data


    def __build_calo_truth(self):
        """return list of calorimeter truth hit objects
        """
        data = []

        for i in range(self.datastruct.truecalo_time.size()):
            t = self.datastruct.truecalo_time[i]
            x = self.datastruct.truecalo_x[i]
            y = self.datastruct.truecalo_y[i]
            z = self.datastruct.truecalo_z[i]
            e = self.datastruct.truecalo_energy[i]
            ch = EDM.calo_truth_hit(t,x,y,z,e)

            id = self.datastruct.truecalo_id[i]
            type = self.datastruct.truecalo_type[i]
            m = self.datastruct.truecalo_module[i]
            s = self.datastruct.truecalo_side[i]
            c = self.datastruct.truecalo_column[i]
            r = self.datastruct.truecalo_row[i]
            w = self.datastruct.truecalo_wall[i]
            ch.set_info(id,type,m,s,c,r,w)
            data.append(ch)

        return data


    def __build_vertex(self):
        """return list of vertex objects
        """
        data = []

        t = self.datastruct_truth.truevertex_time
        x = self.datastruct_truth.truevertex_x
        y = self.datastruct_truth.truevertex_y
        z = self.datastruct_truth.truevertex_z
        
        v = EDM.truevertex(x,y,z,t)
        data.append(v)
        return data


    def __build_trueparticle(self):
        """return list of true particle objects
        """
        data = []
        count = self.datastruct_truth.trueparticle_noparticles
        for i in range(count):
            id = self.datastruct_truth.trueparticle_id[i]
            type = self.datastruct_truth.trueparticle_type[i]
            t = self.datastruct_truth.trueparticle_time[i]
            px = self.datastruct_truth.trueparticle_px[i]
            py = self.datastruct_truth.trueparticle_py[i]
            pz = self.datastruct_truth.trueparticle_pz[i]
            en = self.datastruct_truth.trueparticle_kinenergy[i]

            p = EDM.trueparticle(id,type,px,py,pz,t,en)
            data.append(p)
        return data



class ToySimRun(object):
    """TotSim file decorator.
    Must be a ROOT file containing TTree.

    It constructs a dictionary based event from the data so that later
    algorithms may esaily extract needed raw and truth data, as well as
    add their own data to the event. The events can thus be used in a
    pipeline type analysis.

    It converts the entries in a Toysimulation root
    output file to an event structure:

    event -> dict:
        key : "raw" -> dict
            key : "gg" -> tracker hits

    NOTE: Since the toy MC events generated by TrackGen involve no
    physics processes, the PDG code of each track is set to 0.

    Individual events can be extracted via the get_event method, and all
    events can be iterated over via the events generator method.
    """
    def __init__(self, filename):
        """Construct an instance to read file filename
        """
        self._file = ROOT.TFile(filename,"open")
        self._file.cd()
        # get tree and wrap content structure in datastruct for filling
        r  = tsimreader(self._file)
        self.tree = r.tree
        self.datastruct = r.ds


    def number_of_events(self):
        """Return number of events in the file
        """
        return self.tree.GetEntries()


    def get_event(self, index):
        """Return the event at position index
        """
        event = { }
        self.datastruct.dirx = ROOT.std.vector('double')()
        self.datastruct.diry = ROOT.std.vector('double')()
        self.datastruct.dirz = ROOT.std.vector('double')()
        self.datastruct.pointx = ROOT.std.vector('double')()
        self.datastruct.pointy = ROOT.std.vector('double')()
        self.datastruct.pointz = ROOT.std.vector('double')()
        self.datastruct.breakpointx = ROOT.std.vector('double')()
        self.datastruct.breakpointy = ROOT.std.vector('double')()
        self.datastruct.bpangle = ROOT.std.vector('double')()
        self.datastruct.radius = ROOT.std.vector('double')()
        self.datastruct.wirex = ROOT.std.vector('double')()
        self.datastruct.wirey = ROOT.std.vector('double')()
        self.datastruct.wirez = ROOT.std.vector('double')()
        self.datastruct.grid_id = ROOT.std.vector('int')()
        self.datastruct.grid_side = ROOT.std.vector('int')()
        self.datastruct.grid_layer = ROOT.std.vector('int')()
        self.datastruct.grid_column = ROOT.std.vector('int')()
        self.datastruct.break_layer = ROOT.std.vector('int')()

        self.tree.SetBranchAddress('dirx', ROOT.AddressOf(self.datastruct, "dirx") )
        self.tree.SetBranchAddress('diry', ROOT.AddressOf(self.datastruct, "diry") )
        self.tree.SetBranchAddress('dirz', ROOT.AddressOf(self.datastruct, "dirz") )
        self.tree.SetBranchAddress('pointx', ROOT.AddressOf(self.datastruct, "pointx") )
        self.tree.SetBranchAddress('pointy', ROOT.AddressOf(self.datastruct, "pointy") )
        self.tree.SetBranchAddress('pointz', ROOT.AddressOf(self.datastruct, "pointz") )
        self.tree.SetBranchAddress('breakpointx', ROOT.AddressOf(self.datastruct, "breakpointx") )
        self.tree.SetBranchAddress('breakpointy', ROOT.AddressOf(self.datastruct, "breakpointy") )
        self.tree.SetBranchAddress('bpangle', ROOT.AddressOf(self.datastruct, "bpangle") )
        self.tree.SetBranchAddress('radius', ROOT.AddressOf(self.datastruct, "radius") )
        self.tree.SetBranchAddress('wirex', ROOT.AddressOf(self.datastruct, "wirex") )
        self.tree.SetBranchAddress('wirey', ROOT.AddressOf(self.datastruct, "wirey") )
        self.tree.SetBranchAddress('wirez', ROOT.AddressOf(self.datastruct, "wirez") )
        self.tree.SetBranchAddress('grid_id', ROOT.AddressOf(self.datastruct, "grid_id") )
        self.tree.SetBranchAddress('grid_side', ROOT.AddressOf(self.datastruct, "grid_side") )
        self.tree.SetBranchAddress('grid_layer', ROOT.AddressOf(self.datastruct, "grid_layer") )
        self.tree.SetBranchAddress('grid_column', ROOT.AddressOf(self.datastruct, "grid_column") )
        self.tree.SetBranchAddress('break_layer', ROOT.AddressOf(self.datastruct, "break_layer") )


        self.tree.GetEntry(index) # now self.datastruct is filled

        event["raw"] = {}
        event["raw"]["gg"] = self.__build_observed_gg()
        event["raw"]["truthsim"] = self.__build_truth()

        return event
    
    def events(self):
        """generator yielding all events in the file in sequence
        """
        for i in range(0, self.number_of_events()):
            yield self.get_event(i)


    def __build_observed_gg(self):
        """return list of tracker hit data.
        """
        hits = []
        for i in range(self.datastruct.wirex.size()):
            x = self.datastruct.wirex[i]
            y = self.datastruct.wirey[i]
            z = self.datastruct.wirez[i]
            dz = 7.0 # [mm]
            r = self.datastruct.radius[i]
            dr = 0.9 # [mm]
            gg = EDM.tracker_hit(x,y,z,dz,r,dr) # from edm

            id = self.datastruct.grid_id[i]
            tid = 0
            m = 0
            s = self.datastruct.grid_side[i]
            l = self.datastruct.grid_layer[i]
            c = self.datastruct.grid_column[i]
            gg.set_info(id,tid,m,s,l,c)

            hits.append(gg)
            
        return hits


    def __build_truth(self):
        """return list of truth objects used to create tracker hit data.
        """
        hits = []
        bpcontainer = []
        for i in range(self.datastruct.break_layer.size()):
            bpx = self.datastruct.breakpointx[i]
            bpy = self.datastruct.breakpointy[i]
            bpa = self.datastruct.bpangle[i]
            bpl = self.datastruct.break_layer[i]
            bpcontainer.append((bpx,bpy,bpa,bpl))

        for i in range(self.datastruct.dirx.size()):
            dirx = self.datastruct.dirx[i]
            diry = self.datastruct.diry[i]
            dirz = self.datastruct.dirz[i]
            pointx = self.datastruct.pointx[i]
            pointy = self.datastruct.pointy[i]
            pointz = self.datastruct.pointz[i]
            
            par = (dirx, diry, dirz, pointx, pointy, pointz)
            gg = EDM.toytruth(i,par,bpcontainer) # from edm

            hits.append(gg)
                
        return hits



        
            
class EventRun(object):
    """basic EventFileReader
    Must be a ROOT file containing TTrees with an event dictionary key.
    
    It constructs a dictionary based event from the data so that later 
    algorithms may easily extract processsed data, as well as add 
    their own data to the event. The events can thus be used in a pipeline
    type analysis

    Individual events can be extracted via the get_event method, and all events
    can be iterated over via the events generator method.
    """
    def __init__(self, filename):
        """construct an instance to read file filename.
        """
        self._file = ROOT.TFile(filename,"open")
        self._file.cd()
        self._listTree = ROOT.TTree()
        self._pathTree = ROOT.TTree()
        self._clusterTree = ROOT.TTree()

        self.cluster_flag = False
        self.path_flag = False
        self.list_flag = False
        self.nEvents = 0 # All trees should have the same number of events

        if (self._file.FindObjectAny("cluster_tree")):
            self._clusterTree = self._file.Get("cluster_tree")
            self._clusterTree.SetDirectory(self._file)
            self.clusterTList = self._clusterTree.GetUserInfo()
            if self.clusterTList.GetEntries()>0: # empty UserInfo case
                self.nEvents = self.clusterTList.FindObject("nevents").GetVal()
            self.cluster_flag = True

        if (self._file.FindObjectAny("path_tree")):
            self._pathTree = self._file.Get("path_tree")
            self._pathTree.SetDirectory(self._file)
            self.pathTList = self._pathTree.GetUserInfo()
            if self.pathTList.GetEntries()>0: # empty UserInfo case
                self.nEvents = self.pathTList.FindObject("nevents").GetVal()
            self.path_flag = True

        if (self._file.FindObjectAny("list_tree")):
            self._listTree = self._file.Get("list_tree")
            self._listTree.SetDirectory(self._file)
            self.listTList = self._listTree.GetUserInfo()
            if self.listTList.GetEntries()>0: # empty UserInfo case
                self.nEvents = self.listTList.FindObject("nevents").GetVal()
            self.list_flag = True



    def number_of_events(self):
        """return number of events in the file
        """
        return self.nEvents

    def get_event(self, index):
        """return the event at position index
        """
        event = {}

        if (self.cluster_flag):
            keytemp = self.clusterTList.First().GetValue("cluster_tree")
            s = str(keytemp) # which outbox key built the cluster
            event[s] = self._build_cluster_dict(index,s)

        if (self.path_flag):
            keytemp = self.pathTList.First().GetValue("path_tree")
            s = str(keytemp) # which outbox key built the path list
            event[s] = self._build_path_list(index,s)

        if (self.list_flag):
            keytemp = self.listTList.First().GetValue("list_tree")
            s = str(keytemp) # which outbox key built the gcylinder list
            event[s] = self._build_geiger_list(index,s)

        
        return event


    def _build_cluster_dict(self, index, s):
        """return cluster dictionary from tree
        """
        stream = clusterToRoot(s)
        cld = stream.ReturnCluster(self._clusterTree, index) # writing module knows how to read
        return cld


    def _build_path_list(self, index, s):
        """return sequence of paths from tree
        """
        stream = pathToRoot(s)
        hl = stream.ReturnPaths(self._pathTree, index) # writing module knows how to read
        return hl


    def _build_geiger_list(self, index, s):
        """return sequence of gcylinder from tree
        """
        stream = listToRoot(s)
        hl = stream.ReturnList(self._listTree, index) # writing module knows how to read
        return hl


    def events(self):
        """generator yielding all events in the file in sequence
        """
        for i in range(self.number_of_events()):
            yield self.get_event(i)





def create_reader_instance(filename):
	"""Return a reader instance appropriate
	for the file provided.

	.root - flsimRun decorator
	.tsim - ToySimRun decorator
	"""
        if isinstance(filename, list): # only for [.root, .evt] filelist
            readerlist = []
            rootname = filename[0]
            fname = filename[1]
            if re.search('.root$', rootname):
		readerlist.append(flsimRun(rootname))
                if re.search('.evt$', fname):
                    readerlist.append(EventRun(fname))
                    return readerlist
                else:
                    return None
            else:
                return None
        else:
            if re.search('.tsim$', filename):
		return ToySimRun(filename)
            elif re.search('.root$', filename):
		return flsimRun(filename)
            elif re.search('.evt$', filename):
		return EventRun(filename)
            else:
		return None

