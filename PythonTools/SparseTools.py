import ROOT
import functools
import uproot
import numpy as np

ROOT.gROOT.ProcessLine(".L drawsparse.C")


##
## move this to PythonTools
def default_kwargs(**defaultKwargs):
    """
    wrapper to set default keyword argument for the function
    https://stackoverflow.com/questions/15301999/default-arguments-with-args-and-kwargs
    """
    
    def actual_decorator(fn):
        @functools.wraps(fn)
        def g(*args, **kwargs):
            defaultKwargs.update(kwargs)
            return fn(*args, **defaultKwargs)
        return g
    return actual_decorator




##
##  TODO: switch to using root_numpy
def fillTH2FromList(x,y, bins, name ="Histo", func=None):
    h = ROOT.TH2D(name,name,*bins)
    ntot = len(x)
    assert ntot == len(y)
    
    for ibin in range(1,ntot):
        xval = x[ibin]
        yval = y[ibin]
        if func:
            xval, yval = func(xval, yval)
        h.Fill(xval,yval)
    return h


def getTH2FromList(x,y,bins, name='hist', title='hist'):
    import root_numpy
    """
        uses np.histogram2d to convert the values in x,y arrays to 2d histogram as np.array
        requires root_numpy, 
        careful that the Entries will not be the number of data points, but the number of 'bins'
        If this is important, use fillTH2DFromList (much slower)
    """
    nx, xmin, xmax, ny, ymin, ymax = bins
    nph = np.histogram2d(x, y, (nx,ny), ( (xmin,xmax), (ymin,ymax))  )
    h2d = ROOT.TH2F( name, title, *bins )
    h2d = root_numpy.array2hist( nph[0], h2d )
    return h2d






def getTHnSparse(di, keys_bins, name="MyTHnSparse"):
    """
        di: a dict of the uproot arrays
        key_bins: keys of the di to be used, and the corresponding bins.
                 np.array([
                            ('x', (100,0,10)),
                            ('y', (300,0,10)),
                            ('z', (200,-10,10)),
                 ])
                 
        
        also adds attributes
    
    """
    
    from array import array
    keys_bins = np.array( keys_bins )
    bins_all = keys_bins.T[1]
    names = keys_bins.T[0]
    
    nbins = np.vstack( bins_all ).T[0].astype('int')
    xmin  = np.vstack( bins_all ).T[1]
    xmax  = np.vstack( bins_all ).T[2]

    vals = np.array([di[x] for x in names])
    
    thn = ROOT.THnSparseD(name, name, len(nbins), array('i',nbins), array('d',xmin), array('d',xmax) )
    for vals in np.array(vals.T):
        thn.Fill( array('d', vals) )
    for i,name in enumerate(names):
        thn.GetAxis(i).SetName(name)
        thn.GetAxis(i).SetTitle(name)
        setattr(thn, name, thn.GetAxis(i) )   # just for shortcut, is lost if Clone()'ed
        setattr(thn, f"index_{name}", i)      # 
    return thn
   


    
    
    
    
class SparseUp():
    
    def __init__(self, name):
        """
            upr: uproot tree
        
        """
        #self.upr = upr # dont want to store the upr
        self.name = name
        
        
    def get_arrays(self, upr, keys, basename="", idx=None):
        """
            Get arrays of tree branches,
            
            
            keys: list of branch identifiers
            branchname will be f"{basename}_{key}" if basename specified
            
            idx: array of indicies. if given, it will be used to mask the output array 
        
        """
        
        di   = {}
        #upr  = upr if upr else self.upr
        
        col_base = f"{basename}_" if basename else basename
        cols = upr.arrays([f"{col_base}{k}" for k in keys], outputtype=tuple)
        for k,c in zip(keys,cols):
            if type(idx) == type(None):
                di[k]=c
            else:
                #di[k]=c[idx.astype(int)]
                if not idx.dtype in ( np.dtype('bool'), np.dtype('int') ):
                    raise TypeError("index array (idx) has have dtype int or bool. but it was %s. You can use idx.astype('int' or 'bool')"%idx.dtype)
                di[k]=c[idx]
        return di


def arrayIsNone(array):
    if type(array) in ( type(None), ):
        return True
    if hasattr(array, 'any'):
        return not array.any()






