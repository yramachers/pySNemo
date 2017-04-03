from math import pi
class Interval(object):
    """
    Represents an interval. 
    Defined as half-open interval [start,end), which includes the start 
    position but not the end.
    Start and end do not have to be numeric types. 
    """
    
    
    def __init__(self, start, end):
        """Construct, start must be <= end."""
        if start > end:
            raise ValueError('Start (%s) must not be greater than end (%s)' % (start, end))
        self._start = start
        self._end = end
        
         
    start = property(fget=lambda self: self._start, doc="The interval's start")
    end = property(fget=lambda self: self._end, doc="The interval's end")
     

    def __str__(self):
        return '[%s,%s)' % (self.start, self.end)
    
    
    def __repr__(self):
        return '[%s,%s)' % (self.start, self.end)
    
    
    def __cmp__(self, other):
        if None == other:
            return 1
        start_cmp = cmp(self.start, other.start)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.end, other.end)


    def checkforpihalf(self,d):
        sign = 1
        if (d < -pi/2.0):
            d += pi
            sign = -1
        elif (d > pi/2.0):
            d -= pi
            sign = -1
        return sign,d


    def intersection(self, other):
        """Intersection. @return: An empty intersection if there is none."""
        if self > other:
            other, self = self, other
        if self.end <= other.start:
            return Interval(self.start, self.start)
        return Interval(other.start, self.end)


    def hull(self, other):
        """@return: Interval containing both self and other."""
        if self > other:
            other, self = self, other
        return Interval(self.start, other.end)
    

    def overlap(self, other):
        """@return: True iff self intersects other."""
        if self > other:
            other, self = self, other
        return self.end > other.start
         

    def angle_overlap(self, other):
        """@return: True iff self intersects other on a ring in phi."""
        # bring to right hemisphere
        swapped1=0
        swapped2=0
        signl,temp_start = self.checkforpihalf(self.start)
        signr,temp_end = self.checkforpihalf(self.end)
        if (signl<0 and signr<0):
            if (temp_end<=temp_start):
                temp = Interval(temp_end,temp_start)
            else:
                temp = Interval(temp_start,temp_end)
            swapped1=1
        signol,o_start = self.checkforpihalf(other.start)
        signor,o_end = self.checkforpihalf(other.end)
        if (signol<0 and signor<0):
            if (o_end<=o_start):
                tempo = Interval(o_end,o_start)
            else:
                tempo = Interval(o_start,o_end)
            swapped2=1
        if (swapped1 and swapped2):
            return temp.overlap(tempo)
        elif (swapped1 and not swapped2):
            return 0
        elif (swapped2 and not swapped1):
            return 0
        else: 
            return self.overlap(other)

    def __contains__(self, item):
        """@return: True iff item in self."""
        return self.start <= item and item < self.end
         
    def has(self, item):
        """@return: True iff item in self, inclusive."""
        return self.start <= item and item <= self.end
         

    def subset(self, other):
        """@return: True iff self is subset of other."""
        return self.start >= other.start and self.end <= other.end
         

    def empty(self):
        """@return: True iff self is empty."""
        return self.start == self.end

    def midinterval(self):
        """ returns the mid interval value for numeric types."""
        return (self.end + self.start)*0.5

    def separation(self, other):
        """@return: The distance between self and other."""
        if self > other:
            other, self = self, other
        if self.end > other.start:
            return 0
        else:
            return other.start - self.end
        
## end of http://code.activestate.com/recipes/576816/ }}}

