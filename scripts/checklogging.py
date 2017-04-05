#!/usr/bin/env python

# control modules
import logging
import logging.handlers
from pysnemo.control.event import Event, EventLoop
from pysnemo.control.pipeline import Pipeline
# reconstruction modules
from pysnemo.cluster.control import ImageSegService
from pysnemo.cellular_automaton.control import CAService
from pysnemo.ca_glue.control import GlueService, RoadService


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
	

def print_number_of_hits(event):
	''' 
	This function prints the number of geiger cells in the
	raw hit selection of the event as read from disk.
	'''
	d = event.getKeyValue('raw')
	hits = d['gg']
	print 'Number of geiger entries: ',len(hits)


def print_number_of_clusters(event):
	''' 
	This function prints the number of clusters for key cluster_out.
	'''
	cluster_dict = event.getKeyValue('cluster_out') # returns dict
	print 'ImSeg: Found %d clusters'%len(cluster_dict)
	for k,val in cluster_dict.iteritems():
		print 'No of entries in cluster = %d' % len(val)
		for hit in val: # list tracker_hits
			print hit.meta_info
				

def print_number_of_entries(event):
	''' 
	This function prints the number of entries in all available keys 
	'''
	d = event.getKeyValue('raw')
	print 'Number of dictionaries: ',len(d)
	if 'gg' in d:
		hits = d['gg']
		print 'Number of gg entries: ',len(hits)
	if 'ggtruth' in d:
		hits = d['ggtruth']
		print 'Number of gg truth entries: ',len(hits)
	if 'calo_hits' in d:
		hits = d['calo_hits']
		print 'Number of calo entries: ',len(hits)
	if 'calotruth' in d:
		hits = d['calotruth']
		print 'Number of calo truth entries: ',len(hits)
	if 'muonpaddles' in d:
		hits = d['muonpaddles']
		print 'Number of muon paddle entries: ',len(hits)
	v = d['vertex']
	print 'Number of vertex entries: ',len(v)
	pa = d['trueparticle']
	print 'Number of prim particles entries: ',len(pa)

 

def check_sweep(event):
	cluster_dict = event.getKeyValue('sweeping_out') # returns dict
	cluster_list = event.getKeyValue('sweeping_out_noise') # returns list
	if len(cluster_dict)<1 and len(cluster_list)<1: # empty CA
		cls = event.getKeyValue('cluster_out') # returns pre-cluster
		event.setKeyValue('sweeping_out', cls) # results for writer
		return event



def print_sweeping(event):
	''' 
	This function prints the number of clusters after sweeping.
	'''
	cluster_dict = event.getKeyValue('sweeping_out') # returns dict
	print 'Found %d clusters after sweeping'%len(cluster_dict)
	for k,val in cluster_dict.iteritems():
		print 'SWEEP: cluster key = %s, No of hits in cluster = %d' % (k,len(val))
		for n in val:
			print n.meta_info
	cluster_list = event.getKeyValue('sweeping_out_noise') # returns list
	print 'Found %d noise events after sweeping'%len(cluster_list)
	#for n in cluster_list:
		#print n



if __name__ == '__main__':
	# log all processes in file
	logname = 'testcluster.log'
	setupLog(logname)

	# Specify the file we want to run over
	readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo_se82_1000.root'
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo106_e1MeV.root'


	# Instantiate a dbscan reconstruction object
	imseg = ImageSegService(None,'cluster_out')

	# set up CA Service
	ca1 = CAService('cluster_out', 'ca1_out', 0, 46.0, 2.0, True) # x
	ca1.configure(cell_radius=63.0)     # adapt to snemo
	ca1.configure(cell_radius_max=126.0)
	ca2 = CAService('cluster_out', 'ca2_out', 1, 46.0, 2.0, True) # y
	ca2.configure(cell_radius=63.0)     # adapt to snemo
	ca2.configure(cell_radius_max=126.0)
	
	# set up Glue Service
	glue = GlueService('ca1_out','ca2_out','glue_out')
	
	# set up Road Service
	sweep = RoadService('glue_out','sweeping_out')

	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_number_of_hits,
			    imseg,
			    print_number_of_clusters,
			    ca1,
			    ca2,
			    glue,
			    sweep,
			    check_sweep,
			    print_sweeping)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=1)

	# Run the event loop
	loop.run()
