#!/usr/bin/env python
import networkx as nx
import numpy as np
from scipy.spatial import KDTree

# control modules
import logging
import logging.handlers
from pysnemo.control.stream_event import Event, EventLoopToFile
from pysnemo.control.pipeline import Pipeline
from pysnemo.control.output.ptdstorage import ptdToRoot
# reconstruction modules
from pysnemo.cluster.control import FilamentService
from pysnemo.cellular_automaton.control import CAService
from pysnemo.ca_glue.control import GlueService, RoadService
from pysnemo.graphtrack.control import TrackingService
from pysnemo.fitter.control import FittingService
from pysnemo.reconstruction.control import FitExtrapolatorService
from pysnemo.reconstruction.control import ConsolidatorService
from pysnemo.io.edm import Particle



def setupLog(fname):
	'''
	Setting up logging to file for the pipeline using the standard
	library 'logging' module in python.
	'''
	logger = logging.getLogger('eventloop')
	logger.setLevel(logging.DEBUG)

	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')

	# Add the log message handler to the logger
	LOG_FILENAME = fname
	handler = logging.handlers.RotatingFileHandler(
              LOG_FILENAME, maxBytes=20000, backupCount=3)
	handler.setFormatter(formatter)

	# add the handlers to the logger
	logger.addHandler(handler)

	logger.info('*** Start ***')

	
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
	

def print_truth(event):
	''' 
	This function prints the comparison of fit parameter with truth input
	'''
	d = event.getKeyValue('raw')
	if 'truthsim' in d:
		truth = d['truthsim']
		for toytruth in truth:
			print toytruth # helix:(px,py,charge,refx,refy,refz)


def print_calo(event):
	''' 
	This function prints the calorimeter data
	'''
	d = event.getKeyValue('raw')
	if 'calo_hits' in d:
		truth = d['calo_hits']
		for toytruth in truth:
			print toytruth


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
	#print 'Sweeping has %d clusters'%len(cluster_dict)
	#print 'Sweeping has %d noise hits'%len(noise)
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
		print 'FITTER: Found %d results'%len(fittuple)
		for entry in tup:
			for pnumber, val in entry.fitterpaths.iteritems(): # val a list of fitters
				for item in val:
					if isinstance(item,tuple):
						print '\nCluster key: %d, Candidates key = %s, Chisq = %f' % (entry.id,pnumber,item[2])
						for i in range(len(item[0])):
							print 'Fit parameter %d: %f +- %f'%(i,item[0][i],item[1][i])
					else: # a HelixFit object
						print item


def print_consolidator(event):
	''' 
	This function prints the consolidator output
	'''
	final = event.getKeyValue('final_out') # returns dict
	if final is not None:
		sump = 0
		sumc = 0
		for entry in final:
			if isinstance(entry,Particle):
				sump += 1
			else:
				sumc += 1
			print entry
		print '%d particle(s), %d calo hits'%(sump,sumc)




if __name__ == '__main__':
	# log all processes in file
	logname = 'testcluster.log'
	setupLog(logname)

	# Specify the file we want to run over
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/helix_calo.tsim'
	readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/multiscatter_calo.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/illumination.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/Vvertex_calo.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/leftright_calo.tsim'

	outfile = '/tmp/validation.evt'

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
	#fitter = FittingService('track_out','fit_out',Bfield=0.0025)
	fitter = FittingService('track_out','fit_out',Bfield=0.0)

	# Instantiate a FitExtrapolator
	extra = FitExtrapolatorService('fit_out','extra_out')
	final = ConsolidatorService('extra_out','final_out', 1) # takes cut off

	# writer
	writer = ptdToRoot('final_out') 

	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_truth,
			    print_calo,
			    fil,
			    ca1,
			    ca2,
			    glue,
			    sweep,
			    check_sweep,
			    noise_to_cluster,
			    track,
			    fitter,
			    extra,
			    final,
			    print_consolidator)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoopToFile(readfile, outfile, operation=pipeline, write_operation=[writer], first_event=0, last_event=999)

	# Run the event loop
	loop.run()
