"""
https://root-forum.cern.ch/t/automatic-tlegend-positioning/23844

"""

import ROOT
#import itertools.chain



# Some convenience function to easily iterate over the parts of the collections


# Needed if importing this script from another script in case TMultiGraphs are used
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)


# Start a bit right of the Yaxis and above the Xaxis to not overlap with the ticks
def_start, def_stop = 0.13, 0.89
def_x_width, def_y_width = 0.3, 0.2

##
##https://root-forum.cern.ch/t/from-user-to-ndc/9943/5
##double GetNDC(double x) {
##  gPad->Update();//this is necessary!
##  return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
## }##


def fixLegend(leg):
    leg.SetBorderSize(0)
    leg.SetShadowColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    #leg.SetMargin(0.05)


def transform_to_user(canvas, x1, y1, x2, y2):
    """
    Transforms from Pad coordinates to User coordinates.

    This can probably be replaced by using the built-in conversion commands.
    """
    xstart = canvas.GetX1()
    xlength = canvas.GetX2() - xstart
    xlow = xlength * x1 + xstart
    xhigh = xlength * x2 + xstart
    if canvas.GetLogx():
        xlow = 10**xlow
        xhigh = 10**xhigh

    ystart = canvas.GetY1()
    ylength = canvas.GetY2() - ystart
    ylow = ylength * y1 + ystart
    yhigh = ylength * y2 + ystart
    if canvas.GetLogy():
        ylow = 10**ylow
        yhigh = 10**yhigh

    return xlow, ylow, xhigh, yhigh


def overlap_h(hist, x1, y1, x2, y2):
    xlow = hist.FindFixBin(x1)
    xhigh = hist.FindFixBin(x2)
    for i in range(xlow, xhigh + 1):
        val = hist.GetBinContent(i)
        # Values
        if y1 <= val <= y2:
            return True
        # Errors
        if val + hist.GetBinErrorUp(i) > y1 and val - hist.GetBinErrorLow(i) < y2:
            # print "Overlap with histo", hist.GetName(), "at bin", i
            return True
    return False


def overlap_rect(rect1, rect2):
    """Do the two rectangles overlap?"""
    if rect1[0] > rect2[2] or rect1[2] < rect2[0]:
        return False
    if rect1[1] > rect2[3] or rect1[3] < rect2[1]:
        return False
    return True

def overlap_g(graph, x1, y1, x2, y2):
    x_values = list(graph.GetX())
    y_values = list(graph.GetY())
    x_err = list(graph.GetEX()) or [0] * len(x_values)
    y_err = list(graph.GetEY()) or [0] * len(y_values)

    for x, ex, y, ey in zip(x_values, x_err, y_values, y_err):
        # Could maybe be less conservative
        if overlap_rect((x1, y1, x2, y2), (x - ex, y - ey, x + ex, y + ey)):
            # print "Overlap with graph", graph.GetName(), "at point", (x, y)
            return True
    return False

#@def_start, def_stop = 0.13, 0.89
#def_x_width, def_y_width = 0.4, 0.3


def places(x_min, y_min, x_max, y_max, x_width, y_width, n=10):
    """Scans from top left to bottom right
    (x_min, y_max-y_width, x_min+x_width, y_width)
    (x_min, y_max-y_width, x_min+x_width, y_width)

    """
    import numpy as np
    for width_fact in np.linspace(1,0.5,10):
        x_width_ = x_width*width_fact
        y_width_ = y_width*width_fact
        for y_min_ in np.linspace(y_max-y_width, y_min, n):
            for x_min_ in np.linspace(x_min, x_max-x_width, n):
                loc = (x_min_, y_min_, x_min_+x_width_, y_min_+y_width)
                loc = [round(x,2) for x in loc]
                #print(loc)
                yield loc
        
        
def place_legend_(canvas=None, x1=None, y1=None, x2=None, y2=None, x_width=None, y_width=None, header="", option="LP", return_leg=True):
    # If position is specified, use that
    canvas = canvas if canvas else ROOT.gPad
    if all(x is not None for x in (x1, x2, y1, y2)):
        return canvas.BuildLegend(x1, y1, x2, y2, header, option)
    # Make sure all objects are correctly registered
    canvas.Update()
    
    def_x1  = round( 0.01 +  ( canvas.GetUxmin()-canvas.GetX1() ) / (canvas.GetX2()-canvas.GetX1()) , 3)
    def_x2  = round(-0.01 +  ( canvas.GetUxmax()-canvas.GetX1() ) / (canvas.GetX2()-canvas.GetX1()) , 3)
    def_y1  = round( 0.01 +  ( canvas.GetUymin()-canvas.GetY1() ) / (canvas.GetY2()-canvas.GetY1()) , 3)
    def_y2  = round(-0.01 +  ( canvas.GetUymax()-canvas.GetY1() ) / (canvas.GetY2()-canvas.GetY1()) , 3)

    x1 = x1 if x1 else def_x1
    x2 = x2 if x2 else def_x2
    y1 = y1 if y1 else def_y1
    y2 = y2 if y2 else def_y2
    
    x_width = x_width if x_width else def_x_width
    y_width = y_width if y_width else def_y_width

    
    PLACES = places(def_x1, def_y1,  def_x2, def_y2, x_width, y_width )

    # Build a list of objects to check for overlaps
    objects = []
    for x in canvas.GetListOfPrimitives():
        if isinstance(x, ROOT.TH1) or isinstance(x, ROOT.TGraph):
            objects.append(x)
        elif isinstance(x, ROOT.THStack) or isinstance(x, ROOT.TMultiGraph):
            objects.extend(x)

    # As a fallback, use the default values, taken from TCanvas::BuildLegend
    args = (0.5, 0.67, 0.88, 0.88, header, option)
    for place in PLACES:
        #place_user = canvas.PadtoU(*place)
        #print(place)
        place_user = transform_to_user(canvas, *place)
        # Make sure there are no overlaps
        if any(obj.Overlap(*place_user) for obj in objects):
            continue
        args = (place[0], place[1], place[2], place[3], header, option)
        #print(args)
        break
    #print("chosen: ", args)
    return canvas.BuildLegend(*args) if return_leg else args


def place_legend(*args,**kwargs):
    leg = place_legend_(*args,**kwargs)
    fixLegend(leg)
    return leg
#def place_legend(canvas, x1=None, y1=None, x2=None, y2=None, x_width=None, y_width=None, header="", option="LP", return_leg=True):
#    # If position is specified, use that
#    if all(x is not None for x in (x1, x2, y1, y2)):
#        return canvas.BuildLegend(x1, y1, x2, y2, header, option)
#
#
#    start   = x1      if x1      else def_start 
#    stop    = x2      if x2      else def_stop 
#    x_width = x_width if x_width else def_x_width
#    y_width = y_width if y_width else def_y_width
#
#    PLACES = [(start, stop - y_width, start + x_width, stop),  # top left opt
#              (start, start, start + x_width, start + y_width),  # bottom left opt
#              (stop - x_width, stop - y_width, stop, stop),  # top right opt
#              (stop - x_width, start, stop, start + y_width),  # bottom right opt
#              (stop - x_width, 0.5 - y_width / 2, stop, 0.5 + y_width / 2),  # right
#              (start, 0.5 - y_width / 2, start + x_width, 0.5 + y_width / 2)]  # left
#    # Make sure all objects are correctly registered
#    canvas.Update()
#
#    # Build a list of objects to check for overlaps
#    objects = []
#    for x in canvas.GetListOfPrimitives():
#        if isinstance(x, ROOT.TH1) or isinstance(x, ROOT.TGraph):
#            objects.append(x)
#        elif isinstance(x, ROOT.THStack) or isinstance(x, ROOT.TMultiGraph):
#            objects.extend(x)
#
#    for place in PLACES:
#        #place_user = canvas.PadtoU(*place)
#        place_user = transform_to_user(canvas, *place)
#        # Make sure there are no overlaps
#        if any(obj.Overlap(*place_user) for obj in objects):
#            continue
#        args = (place[0], place[1], place[2], place[3], header, option)
#        return fixLegend(canvas.BuildLegend(*args)) if return_leg else args
#        #return canvas.BuildLegend(place[0], place[1], place[2], place[3], header, option)
#    # As a fallback, use the default values, taken from TCanvas::BuildLegend
#    args = (0.5, 0.67, 0.88, 0.88, header, option)
#    return fixLegend(canvas.BuildLegend(*args)) if return_leg else args
#    #return canvas.BuildLegend(0.5, 0.67, 0.88, 0.88, header, option)




# Monkey patch ROOT objects to make it all work
ROOT.THStack.__iter__ = lambda self: iter(self.GetHists())
ROOT.TMultiGraph.__iter__ = lambda self: iter(self.GetListOfGraphs())
ROOT.TH1.Overlap = overlap_h
ROOT.TGraph.Overlap = overlap_g
#ROOT.TPad.PadtoU = transform_to_user
ROOT.TPad.PlaceLegend = place_legend
