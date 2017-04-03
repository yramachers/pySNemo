import numpy as np

def gg_to_image(data, tracker_rows, tracker_columns):
	''' 
	Create an image from geiger counter data in the tracker.
	Assumes track_hit data as list in data.
	'''
	folded = tracker_rows * tracker_columns
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

