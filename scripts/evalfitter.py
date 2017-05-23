#!/usr/bin/env python
import networkx as nx
import numpy as np
from scipy.spatial import KDTree

# control modules
from pysnemo.control.event import Event, EventLoop
from pysnemo.control.pipeline import Pipeline
# reconstruction modules
from pysnemo.cluster.control import FilamentService
from pysnemo.cellular_automaton.control import CAService
from pysnemo.ca_glue.control import GlueService, RoadService
from pysnemo.graphtrack.control import TrackingService
from pysnemo.fitter.control import FittingService
from pysnemo.utility.interval import Interval



def print_event_info(event):
	'''
	Any function (or function object) that takes a single
	event-like object as an argument can act as an event
	processor. This one just prints the filename and event
	number.
	'''
	print event
        print "available Keys:"
        for k in event.getKeys():
            print k
	

def check_sweep(event):
	cluster_dict = event.getKeyValue('sweeping_out') # returns dict
	cluster_list = event.getKeyValue('sweeping_out_noise') # returns list
	if len(cluster_dict)<1 and len(cluster_list)<1: # empty CA
		cls = event.getKeyValue('cluster_out') # returns pre-cluster
		event.setKeyValue('sweeping_out', cls) # results for writer
		return event



def connected(noiselist):
	gr = nx.Graph()
	infolist = [trh.meta_info for trh in noiselist] # last in meta info
	infoarray = np.array(infolist)
	#print 'noise info array: ',infoarray # after clustering, no side issue
	layer_column = infoarray[:,4:] # layer number and column number
	lclist = infoarray[:,4:].tolist()
	kdt = KDTree(layer_column)
	pairs = kdt.query_pairs(1.5) #set of tuples of pair indices
	for pair in pairs:
		gr.add_edge(tuple(lclist[pair[0]]),tuple(lclist[pair[1]]))	
		
	gg = nx.connected_component_subgraphs(gr)
	return gg


def noise_to_cluster(event):
	''' 
	This helper function takes the noise list from clustering
	and adds the noise events in as a new cluster if useful.
	'''
	cluster_dict = event.getKeyValue('sweeping_out') # returns dict
	noise = event.getKeyValue('sweeping_out_noise') # returns list
	print 'Sweeping has %d clusters'%len(cluster_dict)
	print 'Sweeping has %d noise hits'%len(noise)
	if len(noise)>=3: # anything else should remain noise for now
		infolist = [trh.meta_info for trh in noise]
		infoarray = np.array(infolist)
		layer_column = infoarray[:,4:] # layer number and column number
		lclist = infoarray[:,4:].tolist()
		maxkey = cluster_dict.keys()[-1]
		k = maxkey
		glist = connected(noise)
		for gr in glist:
			trhlist = []
			if len(gr)>2: # allow a bunch of 3 to be a cluster
				for info in gr.nodes():
					linfo = list(info)
					if linfo in lclist:
						trhlist.append(noise[lclist.index(linfo)])
				k += 1
				cluster_dict[k]=trhlist # add cluster
#				print 'added cluster from noise, len=',len(gr)
		if k>maxkey: # cluster dict changed
			# overwrite existing key
			event.setKeyValue('sweeping_out', cluster_dict)
			return event
		# else do nothing



def print_track_info(event):
	''' 
	This function prints the number of tracks and some info
	from graphtrack.
	'''
	data = event.getKeyValue('track_out') # returns dict
	if isinstance(data,dict):
		print 'TRACK: Found %d clusters'%len(data)
		for k,val in data.iteritems():
			print 'track key = %s' % str(k)
			print 'TRACK: Found %d track candidates'%len(val)
			for k2,v2 in val.iteritems():
				print 'TRACK: key = %d' % k2
				for entry in v2:
					print entry
				print ''	
	elif isinstance(data,list): # or list of clusterpath objects from file
		print 'TRACK: Found %d clusters'%len(data)
		for entry in data:
			print 'cluster key = %s' % entry.id
			print 'TRACK: Found %d track candidates'%len(entry.paths)
		



def print_fitter(event):
	''' 
	This function prints the fitter output
	'''
	fittuple = event.getKeyValue('fit_out') # returns dict
	for tup in fittuple:
		print 'FITTER: Found %d results'%len(fitlist)
		for entry in tup:
			for pnumber, val in entry.fitterpaths.iteritems(): # val a list of fitters
				for item in val:
					if isinstance(item,tuple):
						print '\nCluster key: %d, Candidates key = %s, Chisq = %f' % (entry.id,pnumber,item[2])
						for i in range(len(item[0])):
							print 'Fit parameter %d: %f +- %f'%(i,item[0][i],item[1][i])
					else: # a HelixFit object
						print item


def print_validation(event):
	''' 
	This function prints the comparison of fit parameter with truth input
	'''
	d = event.getKeyValue('raw')
	if 'truthsim' in d:
		truth = d['truthsim']
		for toytruth in truth:
			print toytruth # helix:(px,py,charge,refx,refy,refz)

	fittuple = event.getKeyValue('fit_out') # returns dict
	for tup in fittuple:
		for entry in tup:
			for pnumber, val in entry.fitterpaths.iteritems(): # val a list of fitters
				for item in val:
					if isinstance(item,tuple):
						bpangles = item[4]
						bestfit = item[0] # a list
						if bpangles and bestfit: # a broken line fit results
							err = item[1]     # a list
							chi2 = item[2]    # a number
							id = (entry.id,pnumber) # id tuple
							Iic = Interval(bestfit[0]-3*err[0], bestfit[0]+3*err[0]) # allow 3 sigma interval
							Islope = Interval(bestfit[1]-3*err[1], bestfit[1]+3*err[1]) # allow 3 sigma interval
							for toytruth in truth:
								par = toytruth.getParameters()
								slope = par[1] # in x-y plane
								ic = par[4] # on y-axis at x=0
								bplist = toytruth.bplist
								if slope in Islope and ic in Iic:
									print 'Broken line validation found: with chi2=%f at (cluster %d, path %d)'%(chi2,id[0],id[1])
									print 'Broken line: angles ',item
									print 'Truth: angle ',bplist # list of tuples
						elif bestfit: # not empty
								err = item[1]     # a list
								chi2 = item[2]    # a number
								id = (entry.id,pnumber) # id tuple
								#Iic = Interval(bestfit[0]-3*err[0], bestfit[0]+3*err[0]) # allow 3 sigma interval
								Islope = Interval(bestfit[1]-3*err[1], bestfit[1]+3*err[1]) # allow 3 sigma interval
								for toytruth in truth:
									par = toytruth.getParameters()
									slope = par[1] # in x-y plane
									#ic = par[4] # on y-axis at x=0
									if slope in Islope:
										print 'Line validation found: with chi2=%f at (cluster %d, path %d)'%(chi2,id[0],id[1])
										print 'best slope = %f +- %f'%(bestfit[1],err[1])
										print 'best ic    = %f +- %f'%(bestfit[0],err[0])
					else: # a HelixFit object
						bestfit = item.fitmomentum     # a list
						err = item.fitmomentum_errors  # a list
						chi2 = item.chi2               # a number
						refpos = item.par              # a list
						referr = item.errors           # a list
						id = (entry.id,pnumber)      # id tuple
						Irefx = Interval(refpos[0]-3*referr[0], refpos[0]+3*referr[0]) # allow 3 sigma interval
						Irefy = Interval(refpos[1]-3*referr[1], refpos[1]+3*referr[1])
						Irefz = Interval(refpos[2]-3*referr[2], refpos[2]+3*referr[2])
						Ipx = Interval(bestfit[0]-3*err[0], bestfit[0]+3*err[0]) # allow 3 sigma interval
						Ipy = Interval(bestfit[1]-3*err[1], bestfit[1]+3*err[1])

						#print 'Helix fit: '
						#print item
						for toytruth in truth:
							par = toytruth.getParameters()
							px = par[0]
							py = par[1]
							refx = par[3]
							refy = par[4]
							refz = par[5]
							if refx in Irefx and refy in Irefy and refz in Irefz:
								print 'Helix reference point fit OK'

								if px in Ipx and py in Ipy:
									print 'Helix validation found: with chi2=%f at (cluster %d, path %d)'%(chi2,id[0],id[1])
									print item



if __name__ == '__main__':
	# Specify the file we want to run over
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/double_atvertex.tsim'
	readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/single_vertex_helix.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/double_Vvertex.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/single_breakp.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/triple_Sshaped.tsim'


	# Instantiate a Filament service object
	fil = FilamentService(None,'cluster_out')

	# set up CA Service
	ca1 = CAService('cluster_out', 'ca1_out', 0, 46.0, 2.0, True) # x
	ca1.configure(cell_radius=63.0)     # adapt to snemo
	ca1.configure(cell_radius_max=126.0)
	ca1.configure(slope_tolerance=0.1)
	ca2 = CAService('cluster_out', 'ca2_out', 1, 46.0, 2.0, True) # y
	ca2.configure(cell_radius=63.0)     # adapt to snemo
	ca2.configure(cell_radius_max=126.0)
	ca1.configure(slope_tolerance=0.1)
	
	# set up Glue Service
	glue = GlueService('ca1_out','ca2_out','glue_out')
	
	# set up Road Service
	sweep = RoadService('glue_out','sweeping_out')

	# Instantiate a tracking service
	track = TrackingService('sweeping_out','track_out',sigma=1.0)

	# Instantiate a Fitting service
	fitter = FittingService('track_out','fit_out',Bfield=0.0025)
	#fitter = FittingService('track_out','fit_out',Bfield=0.0)


	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    fil,
			    ca1,
			    ca2,
			    glue,
			    sweep,
			    check_sweep,
			    noise_to_cluster,
			    track,
			    fitter,
			    print_validation)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=10)

	# Run the event loop
	loop.run()
