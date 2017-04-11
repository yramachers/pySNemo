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
	

def print_tracker_info(event):
	''' 
	This function prints tracker information
	'''
	d = event.getKeyValue('raw')
	if 'gg' in d:
		hits = d['gg']
		print 'Number of gg entries: ',len(hits)
		for cell in hits:
			print 'wire: (%f, %f, %f)'% (cell.x,cell.y,cell.z)
			print 'errors r,z: (%f, %f)'% (cell.sigmar,cell.sigmaz)
			print 'radius: %f'% (cell.r)
			mi = cell.meta_info
			for i in mi:
				print 'Info: %d'%i
	else:
		print 'No gg entries.'
		

def print_tracker_truth_info(event):
	''' 
	This function prints tracker truth information
	'''
	d = event.getKeyValue('raw')
	if 'ggtruth' in d:
		hits = d['ggtruth']
		print 'Number of truth gg entries: ',len(hits)
		for cell in hits:
			print 'start: (%f, %f, %f)'% (cell.x0,cell.y0,cell.z0)
			print 'stop: (%f, %f, %f)'% (cell.x1,cell.y1,cell.z1)
			print 'time: %f'% (cell.time)
			mi = cell.meta_info
			for i in mi:
				print 'Info: %d'%i
	else:
		print 'No gg truth entries.'
		


if __name__ == '__main__':
	# Specify the file we want to run over
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/Data/Se82bb0nu.root'

	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_tracker_info,
			    print_tracker_truth_info)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=0)

	# Run the event loop, which processes the entire file
	loop.run()
