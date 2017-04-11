#!/usr/bin/env python

# control modules
from pysnemo.control.event import Event, EventLoop
from pysnemo.control.pipeline import Pipeline
# reconstruction modules
from pysnemo.cluster.control import ImageSegService
from pysnemo.cellular_automaton.control import CAService
from pysnemo.ca_glue.control import GlueService, RoadService


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
				

def print_gg_entries(event):
	''' 
	This function prints the number of entries in all available keys 
	'''
	d = event.getKeyValue('raw')
	print 'Number of dictionaries: ',len(d)
	if 'gg' in d:
		hits = d['gg']
		for hit in hits:
			print hit.meta_info
	if 'truthsim' in d:
		truth = d['truthsim']
		for toytruth in truth:
			print toytruth


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
	# Specify the file we want to run over
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo_se82_1000.root'
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo106_e1MeV.root'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/python/scripts/validation/double_atvertex.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/double_Vvertex.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/single_vertex_helix.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/single_breakp.tsim'


	# Instantiate an image segmentation object
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
			    print_gg_entries,
			    imseg,
			    print_number_of_clusters,
			    ca1,
			    ca2,
			    glue,
			    sweep,
			    check_sweep,
			    print_sweeping)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=10)

	# Run the event loop
	loop.run()
