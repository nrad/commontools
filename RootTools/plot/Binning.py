''' Class to implement variable binnings
'''

#Standard imports
import numbers

class Binning:

    def __init__(self, l):
        if ( type(l)==type([]) or type(l)==type(()) ) and ( not len(l)==3 or not all( isinstance(x, numbers.Number) for x in l ) ):
            raise ValueError( "Specified %r. Need list or tuple of len 3. For variable bin-widths use Binning.thresholds" % l )
        self.binning = l
        self.binning_is_explicit = False

    @classmethod
    def fromThresholds(cls, thresholds):
        res = cls([1,0,1])
        if not all( isinstance(x, numbers.Number) for x in thresholds ):
            raise ValueError( "Need list or of numbers, got %r" % thresholds )
        res.binning = thresholds
        res.binning.sort()
        res.binning_is_explicit = True
        
        return res
