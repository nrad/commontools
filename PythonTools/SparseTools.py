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







def getTHnSparse(di, keys_bins, name="MyTHnSparse", weight=None):
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

    thn = ROOT.THnSparseD(name, name, len(nbins), array('i',nbins), array('d',xmin), array('d',xmax) )
    thn.Sumw2()
    if not di:
        return thn

    vals = np.array([di[x] for x in names])

    for vals in np.array(vals.T):
        if weight:
            thn.Fill( array('d', vals), weight )
        else:
            thn.Fill( array('d', vals) )
    for i,name in enumerate(names):
        thn.GetAxis(i).SetName(name)
        thn.GetAxis(i).SetTitle(name)
        setattr(thn, name, thn.GetAxis(i) )   # just for shortcut, is lost if Clone()'ed
        setattr(thn, f"index_{name}", i)      #
    return thn
 


    
    
    
    
class SparseUp():
    
    def __init__(self, name, keys_bins=()):
        """
            upr: uproot tree
            name: just a name :)
            keys_bins should look like:
                 keys_bins = np.array([
                            ('x', (100,0,10)),
                            ('y', (300,0,10)),
                            ('z', (200,-10,10)),
                 ])
            
        
        """
        #self.upr = upr # dont want to store the upr
        self.name = name
        self.keys_bins = keys_bins
        self.keys = list(zip(*keys_bins))[0] if len(keys_bins) else ()
        self.bins = list(zip(*keys_bins))[1] if len(keys_bins) else ()
        
    def get_arrays(self, upr, keys=None, basename="", idx=None):
        """
            Get arrays of tree branches,
            
            
            keys: list of branch identifiers
            branchname will be f"{basename}_{key}" if basename specified
            
            idx: array of indicies. if given, it will be used to mask the output array 
        
        """
        
        di   = {}
        keys = keys if keys else self.keys
        #upr  = upr if upr else self.upr
        
        col_base = f"{basename}_" if basename else basename
        col_names = [f"{col_base}{k}" for k in keys]
        #cols = upr.arrays(col_names, outputtype=tuple)
        cols = upr.arrays(col_names, library="np")
        #for k,c in cols.items(): #lzip(keys,cols):
        for k in col_names: #lzip(keys,cols):
            c = cols[k]
            if type(idx) == type(None):
                di[k]=c
            else:
                #di[k]=c[idx.astype(int)]
                #if not idx.dtype in ( np.dtype('bool'), np.dtype('int') ):
                #    raise TypeError("index array (idx) has have dtype int or bool. but it was %s. You can use idx.astype('int' or 'bool')"%idx.dtype)
                di[k]=c[idx]
     
        self.di = di
        return self.di


    def get_thn(self, di=None, bins=None, weight=None, name=None):
        di = di if di else self.di
        bins = bins if bins else self.keys_bins 
        name = name if name else self.name
        thn = getTHnSparse(di, keys_bins=bins, weight=weight, name=name)
        fix_up_thn(thn)
        self.thn = thn
        return thn

    def getProj(self, x,y=None, name=None, thn=None):
        thn = thn if thn else self.thn
        if not y:
            proj = thn.Projection(getattr(thn,"index_%s"%x))
        else:
            proj = thn.Projection(getattr(thn,"index_%s"%x), getattr(thn,"index_%s"%y))
        name = name if name else f"{x}" + ("" if not y else f"_vs_{y}" )
        proj.SetName(name)
        proj.SetTitle(name) 
        return proj

    def setcut(self, ax, min_max=()):
        thn = self.thn
        getattr(thn,ax).SetRange(  *min_max )
        return thn

def arrayIsNone(array):
    if type(array) in ( type(None), ):
        return True
    if hasattr(array, 'any'):
        return not array.any()



def fix_up_thn(thn):
    axes_labels = dict( [ (i,x.GetTitle()) for i,x in enumerate( thn.GetListOfAxes() ) ] )
    for iax, axname in axes_labels.items():
        setattr(thn, axname, thn.GetAxis(iax) )
        setattr(thn, f"index_{axname}", iax)
    setattr(thn, 'axes_labels', axes_labels )
    setattr(thn, 'axes_idx',  {v:k for k,v in axes_labels.items() })





