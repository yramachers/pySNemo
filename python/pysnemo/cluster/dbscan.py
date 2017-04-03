# dbscan - implementation of DBSCAN density based clustering
#
# Copyright (c) 2010, 2011 Ben Morgan <Ben.Morgan@warwick.ac.uk>
# Copyright (c) 2010, 2011 The University of Warwick
#
# This file is part of pySNemo.
#
# pySNemo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pySNemo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with pySNemo.  If not, see <http://www.gnu.org/licenses/>.

"""
DBSCAN Density Based Clustering
===============================

This module implements the DBSCAN (Density Based Spatial Clustering for
Applications with Noise) spatial clustering algorithm:

    A density-based algorithm for discovering clusters in large spatial
    databases with noise
    M. Ester, H.P. Kriegel, S. Joerg and X. Xu
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.71.1980

The implementation follows the pseudo-algorithm outlined in Ester et al,
though also see a potentially simpler definition here:

http://en.wikipedia.org/wiki/DBSCAN

DBSCAN is implemented as a basic function taking the dataset to be
clustered together with the DBSCAN parameters 'epsilon' and 'min_points'.

The function returns a dictionary mapping cluster id to a list of
point indices in the original dataset.

"""

import pysnemo.cluster.dbscan_algorithm as DBS
from numpy import array, where

def dbprocess(data, clusterpar, eps):
    nclass, ttype, e, boolans = DBS.dbscan(data,clusterpar,eps)
    nn=nclass.max()
    return nn, nclass

def dbscan(X, epsilon, min_points):
    """return list of clusters of pointset found by DBSCAN.

    Parameters:
    -----------
    X : iterable of iterables
        Collection of points to be clustered. 
    epsilon : float
        Distance definining neighbourhood around each point.
    min_points : integer
        Minimum number of points required in neighbourhood
        
    Returns:
    --------
    clusters : dict
        Mapping between cluster id and numpy array of indices in dataset X

    """

    data = array(X) # Note that dbscan_algorithm only works with numpy arrays
    nclusters, nclass = dbprocess(data,min_points,epsilon) #for dbscan

    # build output cluster dictionary
    clusters = {}
    for i in range(1,int(nclusters)+1):
        ind = where(i==nclass)
        clusters[i] = ind # clusters numbered from 1 to nclusters

    return clusters


