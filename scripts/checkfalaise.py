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
		

def print_muonpaddle_info(event):
	''' 
	This function prints muon paddle truth information
	'''
	d = event.getKeyValue('raw')
	if 'muonpaddles' in d:
		hits = d['muonpaddles']
		print 'Number of muonpaddle entries: ',len(hits)
		for cell in hits:
			print 'start: (%f, %f, %f)'% (cell.x0,cell.y0,cell.z0)
			print 'stop: (%f, %f, %f)'% (cell.x1,cell.y1,cell.z1)
			print 'start time, stop time: %f, %f'% (cell.starttime,cell.stoptime)
			print 'energy deposit: %f'% (cell.energy)
			mi = cell.meta_info
			print 'muon Id: %d'%mi
	else:
		print 'No muonpaddle entries.'
		

if __name__ == '__main__':
	# Specify the file we want to run over
	#readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/Data/Se82bb0nu.root'
	readfile  = '/storage/epp2/phsdaq/Sandbox/Falaise/trunk/workdir/Data/tracker_muons_10.root'

	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_tracker_info,
			    print_tracker_truth_info,
			    print_muonpaddle_info)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=0)

	# Run the event loop, which processes the entire file
	loop.run()
