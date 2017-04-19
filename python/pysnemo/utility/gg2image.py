import networkx as nx
import numpy as np
from scipy.spatial import KDTree

def gg_to_image(data, tracker_rows, tracker_columns):
	''' 
	Create an image from geiger counter data in the tracker.
	Assumes track_hit data as list in data.
	'''
	left = np.zeros((tracker_rows,tracker_columns)) # left tracker
	right = np.zeros((tracker_rows,tracker_columns)) # left tracker
	for hit in data:    # assumes track_hit data as list in data
		mi = hit.meta_info
		side  = mi[3]
		layer = mi[4]
		row   = mi[5]
		if side < 1: # left of foil
			left[row,8-layer]=1
		else:
			right[row,layer]=1
	return left, right  # return as 2 numpy arrays



def gg_to_single_image(data, tracker_rows, tracker_columns):
	''' 
	Create an image from geiger counter data in the tracker.
	Assumes track_hit data as list in data.
	'''
	im = np.zeros((tracker_rows,tracker_columns)) # left tracker
	for hit in data:    # assumes track_hit data as list in data
		mi = hit.meta_info
		layer = mi[4]
		row   = mi[5]
		im[row,layer]=1
	return im  # return as numpy array



def gg_to_left_right_images(data, tracker_rows, tracker_columns):
	''' 
	Create two images from geiger counter data in the tracker.
	Assumes track_hit data as list in data. No mirroring of left image.
	'''
	left = np.zeros((tracker_rows,tracker_columns)) # left tracker
	right = np.zeros((tracker_rows,tracker_columns)) # left tracker
	for hit in data:    # assumes track_hit data as list in data
		mi = hit.meta_info
		side  = mi[3]
		layer = mi[4]
		row   = mi[5]
		if side < 1: # left of foil
			left[row,layer]=1
		else:
			right[row,layer]=1
	return left, right  # return as 2 numpy arrays



def images_to_gg(data, clleft, clright):
	''' 
	Return proper geiger hits from hitinfo clusters created in image segmentation
	'''
	kleft = len(clleft)
	kright = len(clright)
	clusters = { }

	for n in range(kleft+kright):
		clusters[n+1]=[]
	hitlist = []
	for hit in data:
		mi = hit.meta_info
		side  = mi[3]
		layer = mi[4]
		row   = mi[5]
		test = (side, row, layer) # check hit against all cluster lists
		for k,v in clleft.iteritems():
			if test in v:
				hitlist.append((k,side,hit))
		for k,v in clright.iteritems():
			if test in v:
				hitlist.append((k,side,hit))
	for k,s,h in hitlist:
		if s<1: # left
			clusters[k].append(h)
		else: # right
			clusters[kleft+k].append(h)
	return clusters




def image_to_graph(data): # image data as input
    graph = nx.Graph() # un-directed empty graph to be filled

    pixels = np.transpose(np.nonzero(data)) # nonzero pixels, row,col format
    kdt = KDTree(pixels)
    pairs = kdt.query_pairs(1.5) #set of tuples of pair indices incl diagonal
    for pair in pairs:
        graph.add_edge(tuple(pixels[pair[0]]),tuple(pixels[pair[1]]))
    return graph



def graph_to_image(gr, tracker_rows, tracker_columns): # graph as input
    image = np.zeros((tracker_rows,tracker_columns))
    for nd in gr.nodes(): # nodes are the pixels
        image[nd[0],nd[1]] = 1
    return image


