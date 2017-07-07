#!/usr/bin/env python

# control modules
from pysnemo.control.event import Event, EventLoop
from pysnemo.control.pipeline import Pipeline
# reconstruction modules

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
	

def print_gg_entries(event):
	''' 
	This function prints the number of entries in all available keys 
	'''
	d = event.getKeyValue('raw')
	print 'Number of dictionaries: ',len(d)
	if 'gg' in d:
		hits = d['gg']
		for hit in hits:
			print hit


def print_calo_entries(event):
	''' 
	This function prints the number of entries in all available keys 
	'''
	d = event.getKeyValue('raw')
	print 'Number of dictionaries: ',len(d)
	if 'calo_hits' in d:
		hits = d['calo_hits']
		for hit in hits:
			print hit


if __name__ == '__main__':
	# Specify the file we want to run over
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo106_se82_0nu.root'
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/DemoTestData/demo_se82_1000.root'
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/Data/se82_1000.root'
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/Data/tracker_muons_10.root'
	readfile  = 'validation/2helix_calo.tsim'

	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_gg_entries,
			    print_calo_entries)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=10)

	# Run the event loop, which processes the entire file
	loop.run()
