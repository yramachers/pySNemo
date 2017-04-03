# pca_transformation.py - Utility for pca transformation calculations
#
# Copyright (c) 2012 YR
# Copyright (c) 2011 The University of Warwick
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
pca_transformation: utility for calculating pca
-----------------------------------------------

"""
from numpy import array, linalg, transpose, dot
from scipy import cov

# PCA related helper functions
def PCA_EigenVectors(hits):
    '''
    Input expects hits as a numpy array.
    Utility function:
    from utilities package but here only requests the full set
    of eigenvectors in order to transform the data.
    '''
    # takes data as numpy array
    data_array = transpose(hits)
    etranspose = []

    # Get eigenvalues and eigenvectors
    if (len(data_array) > 0):
        eigenval, eigenvec = linalg.eig(cov(data_array))
        # Transpose eigenvec to return to dataset
        etranspose = transpose(eigenvec)

    return etranspose

def PCA_EigenVectors_Values(fullhits):
    '''
    Input expects hits as a list of lists
    Utility function:
    from utilities package but here only requests the full set
    of eigenvectors in order to transform the data.
    '''
    X1 = array([row[:3] for row in fullhits]) # voxel data only
    # takes data as numpy array
    data_array = transpose(X1)
    # Get eigenvalues and eigenvectors
    eigenval = []
    etranspose = []

    if (len(data_array) > 0):
        eigenval, eigenvec = linalg.eig(cov(data_array))
        # Transpose eigenvec to return to dataset
        etranspose = transpose(eigenvec)

    return eigenval, etranspose


# the PCA transformation function for N>=3 Dim hits
def pca_transformation(fullhits):
    '''
    Analysis: perform pca and transform the data along
    major axis. Then hand back transformed hits
    Input expects hits as list of tuples
    '''

    l = []
    if (len(fullhits) > 0):

        X1 = array([row[:3] for row in fullhits]) # voxel data only
        metadata = [row[3:] for row in fullhits] # all the attached meta data

        ev = PCA_EigenVectors(X1) # get the 3D set of eigenvectors
        # transform
        tX = dot(X1,(ev.T))
        # attach meta data again
        for vox, meta in zip(tX,metadata):
            l.append(tuple(list(vox)+list(meta)))

        # sort along major axis
        l.sort(key=lambda tup: tup[0]) # in-place sort for coordinate 0
        
    return l, ev


def inverse_pca(hits, ev):
    '''
    Analysis: perform the inverse pca transformation.
              Then hand back transformed hits.
    Input expects hits as list of tuples
    '''
    l = []
    if (len(hits) > 0):

        X1 = array([row[:3] for row in hits]) # voxel data only
        metadata = [row[3:] for row in hits] # all the attached meta data

        # inverse transform
        tX = dot(X1,ev)

        # attach meta data again
        for vox, meta in zip(tX,metadata):
            l.append(tuple(list(vox)+list(meta)))
    return l

