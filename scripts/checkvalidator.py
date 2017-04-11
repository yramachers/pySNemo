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


if __name__ == '__main__':
	# Specify the file we want to run over
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/single_breakp.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/double_atvertex.tsim'
	#readfile  = '/home/epp/phsdaq/Code/pySNemo/workdir/validation/validationData/double_Vvertex.tsim'


	# Build a pipeline
	pipeline = Pipeline(print_event_info,
			    print_number_of_entries,
			    print_gg_entries)

	# Create an event loop, giving it the filename and the function to run
	loop = EventLoop(readfile, operation=pipeline, first_event=0, last_event=10)

	# Run the event loop, which processes the entire file
	loop.run()
