''' What is a plot?
'''

# Standard imports
import ROOT
from math import sqrt

#RootTools
from RootTools.plot.PlotBase import PlotBase
from RootTools.core.TreeVariable import ScalarTreeVariable

def addOverFlowBin1D( histo, addOverFlowBin = None):

    if not any( isinstance(histo, o) for o in [ROOT.TH1, ROOT.TProfile]):
        raise NotImplementedError( "Can add overflow bin only to 1D histos. Got %r" %histo )

    if addOverFlowBin is not None and not hasattr(histo, 'overflowApplied'):
        if addOverFlowBin.lower() == "upper" or addOverFlowBin.lower() == "both":
            nbins = histo.GetNbinsX()
            histo.SetBinContent(nbins , histo.GetBinContent(nbins) + histo.GetBinContent(nbins + 1))
            histo.SetBinError(nbins , sqrt(histo.GetBinError(nbins)**2 + histo.GetBinError(nbins + 1)**2))
        if addOverFlowBin.lower() == "lower" or addOverFlowBin.lower() == "both":
            histo.SetBinContent(1 , histo.GetBinContent(0) + histo.GetBinContent(1))
            histo.SetBinError(1 , sqrt(histo.GetBinError(0)**2 + histo.GetBinError(1)**2))
        histo.overflowApplied = True

class Plot( PlotBase ):

    defaultStack           = None
    defaultAttribute       = None
    defaultBinning         = None
    defaultName            = None
    defaultSelectionString = None
    defaultWeight          = None
    defaultHistoClass      = ROOT.TH1F
    defaultTexX            = ""
    defaultTexY            = "Number of Events"
    defaultAddOverFlowBin  = None

    @staticmethod
    def setDefaults(stack = None, attribute = None, binning = None, name = None, selectionString = None, weight = None, histo_class = ROOT.TH1F,
                 texX = "", texY = "Number of events", addOverFlowBin = None):
        Plot.defaultStack           = stack
        Plot.defaultAttribute       = attribute
        Plot.defaultBinning         = binning
        Plot.defaultName            = name
        Plot.defaultSelectionString = selectionString
        #Plot.defaultWeight          = staticmethod(weight) if not isinstance( weight, (list, tuple)) else [map(staticmethod, w2) for w2 in weight ]
        Plot.defaultWeight          = weight #if not isinstance( weight, (list, tuple)) else [map(staticmethod, w2) for w2 in weight ]
        Plot.defaultHistoClass      = histo_class
        Plot.defaultTexX            = texX
        Plot.defaultTexY            = texY
        Plot.defaultAddOverFlowBin  = addOverFlowBin

    def __init__(self, stack = None, attribute = None, binning = None, name = None, selectionString = None, weight = None, histo_class = None,
                 texX = None, texY = None, addOverFlowBin = None, read_variables = []):
        ''' A plot needs a
        'stack' of Sample instances, e.g. [[mc1, mc2, ...], [data], [signal1, signal2,...]], 
        'attribute' Can be string -> getattr( event, attribute),  Variable instance 
        'selectionString' to be used on top of each samples selectionString, 
        'weight' function, 
        'histo_class', e.g. ROOT.TH1F or ROOT.TProfile1D
        'texX', 'texY' labels for x and y axis 
        ''' 

        if isinstance(attribute, (list, tuple)):
            attributes = attribute
        elif attribute is None:
            attributes = Plot.defaultAttribute
        else:
            attributes = [ attribute ]
   
        super(Plot, self).__init__( \
            stack           = stack            if stack           is not None else Plot.defaultStack,
            selectionString = selectionString  if selectionString is not None else Plot.defaultSelectionString,
            weight          = weight           if weight          is not None else Plot.defaultWeight,
            texX            = texX             if texX            is not None else Plot.defaultTexX,
            texY            = texY             if texY            is not None else Plot.defaultTexY,
            name            = name             if name            is not None else Plot.defaultName if Plot.defaultName is not None else attribute.name,
            read_variables  = read_variables,
            attributes      = attributes 
        )

        self.binning         = binning          if binning         is not None else Plot.defaultBinning
        self.histo_class     = histo_class      if histo_class     is not None else Plot.defaultHistoClass
        self.addOverFlowBin  = addOverFlowBin   if addOverFlowBin  is not None else Plot.defaultAddOverFlowBin

    @classmethod
    def fromHisto(cls, name, histos, texX = defaultTexX, texY = defaultTexY):
        res = cls(stack=None, name=name, attribute=None, binning=None, selectionString = None, weight = None, histo_class = None,\
            texX = texX, texY = texY)
        res.histos = histos
        return res
