from pysnemo.utility.rawgeigerhits import raw_hits
from pysnemo.graphtrack.graphpathfinder import cleangraph

class PointCloud(object):
    '''
    Point Cloud out of a collection of tracker_hit objects
    '''
    def __init__(self,data, sigma=3.0, lattice=44.0):
        '''
        Input: Ring data as list of tracker_hit objects.
               Sigma tolerance to merge close tangent points into one
               lattice constant to get the geometry units and unit cell
        
        Output: List of tracker_hit objects available with run() method.
        '''
        self.data = data
        self.sigma = sigma
        self.lattice = lattice


    def run(self):
        rh = raw_hits(self.sigma, self.lattice, self.data)
        hitsdict = rh.make_hitdictionary()

        gg=cleangraph(hitsdict,self.sigma) # sigma tolerance
        #nd=gg.graph.nodes()
        #return nd
        ed=gg.graph.edges()
        return ed
