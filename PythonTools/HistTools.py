from PythonTools.u_float import u_float
from PythonTools.ROOTHelpers import drawRatioPlot
from PythonTools.Legend import fixLegend
from PythonTools.standard import drawLatex, getAndSetEList, getHistoFromSample, saveCanvas, th1Func, th2Func
from PythonTools.StatHelpers import *

import os
import ROOT
#import Workspace.DegenerateStopAnalysis.tools.defaults as defaults
#import Workspace.DegenerateStopAnalysis.tools.tdrstyle as tdrstyle
from copy import deepcopy
import uuid
import numpy as np
from uncertainties import ufloat
import math



#from replot import defaults

#################
##
##    Get Plots
##
##################





def getStackFromHists(histList,sName=None,scale=None, normalize=False, transparency=False):
  #print "::::::::::::::::::::::::::::::::::::::::::: Getting stack" , sNam
  if not sName:
    from PythonTools.NavidTools import uniqueName
    #sName = "stack_%s"%uniqueHash()
    sName = uniqueName("stack_")
  stk=ROOT.THStack(sName,sName)

  if transparency:
    alphaBase=0.80
    alphaDiff=0.70
    alphas=[alphaBase-i*alphaDiff/len(histList) for i in range(len(histList)) ]
    #print alphas
    #print histList

  for i, hist in enumerate(histList):
    #h = hist.Clone()
    h = hist

    #  h.ClearUnderflowAndOverflow()  remove for efficiecy plots
    if scale:
      print ("    Scaling: ", sName if sName else [ hist.GetName(), hist.GetTitle() ])
      h.Scale(scale)
    if normalize:
      if h.Integral():
        h.Scale(1/h.Integral()) 
      else:
        print("Histogram Integral is zero, can't normalize",  sName if sName else [ hist.GetName(), hist.GetTitle()]) 
    if transparency:
      h.SetFillColorAlpha(h.GetFillColor(), alphas[i])
    stk.Add(h)
  return stk


def getStackTot(stack):
    mc_hist = stack.GetHists().Last().Clone(stack.GetName()+"_tot_"  )
    mc_hist.Reset()
    mc_hist.Merge( stack.GetHists() )
    mc_hist.SetDirectory(0)
    return mc_hist


#def getTotalFromStack(stack):
#    hists = stack.GetHists()
#    tot   = hists.Last().Clone("total_"+stack.GetName())
#    tot.Reset()
#    tot.Merge(hists)
#    return tot


def addHists(histList, name=None, title=None):
    stack = getStackFromHists(histList)
    hist_tot = getStackTot(stack)
    if name:
        hist_tot.SetName(name)
    if title:
        hist_tot.SetTitle(title)
    return hist_tot

#def addHists(histlist):
#    stack = ROOT.THStack()
#    for h in histlist:
#        stack.Add(h)
#    return getTotalFromStack(stack)

def addTHns(thns, name=None, title=None):
    name = name if name else "total"
    h0 = thns[0].Clone(name)
    for h in thns[1:]:
        h0.Add(h)
    if title:
        hist_tot.SetTitle(name)
    return h0

def getHistMax(hist, retrunBinVal=False):
    """
        returns the x and y corresponding to maximum
        if retrunBinVal=True will the the value at the center of bin, otherwise the bin number.
    """
    nBinX = hist.GetNbinsX()
    histMax= max( [ (  hist.GetXaxis().GetBinCenter(x) if retrunBinVal else x, hist.GetBinContent(x) ) for x in range(1, nBinX+1)] , key= lambda f: f[1] )
    return histMax

def getHistMin(hist,onlyPos=False):
    nBinX = hist.GetNbinsX()
    binContents = [ (x, hist.GetBinContent(x) ) for x in range(1, nBinX+1)]
    if onlyPos:
        binContents=filter( lambda x: x[1]>0, binContents )
    ret = min( binContents , key= lambda f: f[1] ) if binContents else [0,0]
    return ret


def getHistsMax(hists):
    hist_maxs = [ (h,getHistMax(h)) for h in hists]
    return max( hist_maxs, key=lambda f: f[1][1])

def getHistsMin(hists):
    hist_mins = [ (h,getHistMin(h)) for h in hists]
    return min( hist_mins, key=lambda f: f[1][1])


def getRatio(h1,h2,eff=False):   
    r = h1.Clone()
    if eff:
        r.Divide(r,h2, 1.0, 1.0, "B")
    else:
        r.Divide(h2)
    return r


def getHistErrors(h):
    hnom = h.Clone()
    he   = h.Clone()
    for ib in range(h.GetNbinsX()+1):
        hnom.SetBinError(ib, 0)
    he.Divide(hnom)
    return he



def sortHists(hists, func = lambda x:x.Integral(), reverse=False):
    """
        Sort histograms. By default according to their integral.
    """
    return list(sorted(hists, key=func, reverse=reverse))



def getHistIntegralAndError(h, option=""): 
    """
    option can be '' or 'width' 
    width: multiplies bin content with bin width
    https://root.cern.ch/doc/master/classTH1.html#a72750b362fbf8fae7a4c26579ece18bb
    """
    from numpy import array
    s = array([0.0])
    v = h.IntegralAndError(1, h.GetNbinsX(),s, option )
    #print(v,s[0])
    try:
        from uncertainties import ufloat
        return ufloat(v,s[0])
    except:
        return (v,s[0])



def getFullWidthWrtMedian(hist, pr=0.68, med_prob = 0.5):
    """
       calculates full width of histogram as the difference of  quantiles at  med+pr/2.0 and med-pr/2.0
       
    """
    #import numpy as np
    #med_prob = 0.5
    arr_prob = np.array([med_prob, med_prob - pr/2, med_prob + pr/2])
    q = np.array([0.0]*len(arr_prob))
    
    y=hist.GetQuantiles(len(q), q, arr_prob)
    #return abs((abs(q[2]) - abs(q[1]))/2.0)
    return q[2]-q[1]


def getConfidenceInterval(hist, p=0.68, med_prob = 0.5, double_sided=True):
    """
       calculates confidence interval from the different in quantiles (Q)
        for double sided:    Q(med+p/2) - Q(med-p/2)
        for single sided:    Q(med+p) - Q(med)
       and returns their difference
    """
    #import numpy as np
    #med_prob = 0.5
    if double_sided:
        arr_prob = np.array([med_prob - p/2, med_prob + p/2])
    else:
        arr_prob = np.array([med_prob, med_prob + p])
    q = np.array([0.0]*len(arr_prob))
    
    y=hist.GetQuantiles(len(q), q, arr_prob)
    #return abs((abs(q[2]) - abs(q[1]))/2.0)
    #return abs(q[1]) - abs(q[0])
    return q[1]-q[0]



def integrate(hist, x0, x1):
    """
        Integrate with respect to x values (not bin number)
    """
    x0_bin = hist.GetXaxis().FindBin(x0)-1
    x1_bin = hist.GetXaxis().FindBin(x1)
    return hist.Integral(x0_bin, x1_bin)

def getBootstrappedQuantileError(hist, p=0.95, N=1000, uf=True):
    raise Exception("use getBootstrappedFullWidth")

def getBootstrappedFullWidth(hist, p=0.95, N=1000, half_width=False, uf=True):

    #import numpy as np
    new = hist.Clone()
    entries = hist.GetEntries()
    #entries = 1000
    if entries:
        vals = []
        for i in range(N):
            new.Reset()
            new.FillRandom(hist, int(entries))
            
            width = getFullWidthWrtMedian(new, p)
            if half_width:
                width = width/2.0 
            vals.append( width )
        
        new.Delete()
        res = (np.mean(vals), np.std(vals))
    else:
        res = (0,0)
    if uf:
        from uncertainties import ufloat
        return ufloat( *res  )
    else:
        return res


def getQuantile(hist, pr=0.5):

    arr_prob = np.array( [pr] if not isinstance(pr, (tuple,list)) else pr)
    q = np.array([0.0]*len(arr_prob))
    y=hist.GetQuantiles(len(q), q, arr_prob)
    
    #print(y,q)
    return q


def sigdig(x, n=2):
    """
    round number to n significant figures:
    TODO: make compatible with uncertainties.ufloat
    """
    if x:
        return round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
    else:
        return x


def round_uf(uf):
    from uncertainties import ufloat_fromstr
    return ufloat_fromstr("%s"%uf)




def getHistNearMax(h, per=0.9, reverse=False):
    """
    Try to find the x position where y reaches near 90% of maximum
    """
    from root_numpy import hist2array
    arr = hist2array(h) 
    arr_max = np.amax(arr)
    
    x_ys = list(enumerate(arr))
    if reverse:
        x_ys = x_ys[::-1]
    for x,y in x_ys:
        if y> per* arr_max:
            #print(x,y)
            break
            #return x
    return x


def getPunziFOM(seff, b, z=1.96):
    return seff/(z/2.0+math.sqrt(b))

def punziFOM(s,b, z=1.96, nSigtot=None):
        return (s/nSigTot)/(z/2.0 + math.sqrt(b)  ) if nSigTot else 0

def getPunziFOMfromHists(sig_hist, bkg_hist, z=1.96, nSigTot=None, per=0.99):
    """
        punzi fom based on eff(s)/(z/2 + sqrt(b))
        
        
        z: expected value for a corresponding C.L
        {
            90  : 1.645,
            95  : 1.96,  
            98  : 2.33,
            99.9: 3.291,
        }
        per: percentage of maximum as a suggested point for a cut
    """
    
    nSigTot = nSigTot if nSigTot else sig_hist.Integral()
    
    if not nSigTot:
        print("WARNING signal histogram is empty....")
    
    punziFOM = lambda s,b: (s/nSigTot)/(z/2.0 + math.sqrt(b)  ) if nSigTot else 0
    
    hfom_lt  = th1Func(sig_hist, lambda x, bc:  punziFOM( sig_hist.Integral(0,x) , bkg_hist.Integral(0,x)  ) )
    hfom_gt  = th1Func(sig_hist, lambda x, bc:  punziFOM( sig_hist.Integral(x,  sig_hist.GetNbinsX()+1 ) , bkg_hist.Integral(x,bkg_hist.GetNbinsX()+1 )))
    
    #c=ROOT.TCanvas()
    #hfom.Draw('hist')
    hfom_lt.SetTitle("FOM(<x)")
    hfom_lt.SetLineColor(ROOT.kViolet)

    hfom_gt.SetTitle("FOM(>x)")
    hfom_gt.SetLineColor(ROOT.kRed)
    
    
    if per:
        for h in [hfom_gt, hfom_lt]:
            reverse =  True if h==hfom_gt else False
            
            sug_x = getHistNearMax(h, per=per, reverse=reverse)
            sug = h.GetXaxis().GetBinCenter(sug_x  + (-1 if reverse else +1) )
            #sug = getQuantile(h, 0.90)[0]
            h.quant = sug
            h.line  = ROOT.TLine(sug, 0, sug, h.GetMaximum())
            h.line.SetLineColor(h.GetLineColor())
            h.line.SetLineStyle(h.GetLineStyle())
            h.line.SetLineWidth(h.GetLineWidth())
            
    return hfom_lt, hfom_gt



#def simpleFOM(s,b, bkg_fact=1):
#    return (s/sqrt(s+bkg_fact*b)) if s+b else 0

def getFOMfromHists(sig_hist, bkg_hist, fomFunc=punziFOM, fomArgs={'z':1.96, 'nSigTot':None}, per=0.99):
    """
        punzi fom based on eff(s)/(z/2 + sqrt(b))
        
        
        z: expected value for a corresponding C.L
        {
            90  : 1.645,
            95  : 1.96,  
            98  : 2.33,
            99.9: 3.291,
        }
        per: percentage of maximum as a suggested point for a cut
    """
    from functools import partial
    nSigTot = nSigTot if nSigTot else sig_hist.Integral()
    
    if not nSigTot:
        print("WARNING signal histogram is empty....")
    
    fomFuncWrapper = partial(fomFunc, **fomArgs)
    
    hfom_lt  = th1Func(sig_hist, lambda x, bc:  fomFuncWrapper( sig_hist.Integral(0,x) , bkg_hist.Integral(0,x)  ) )
    hfom_gt  = th1Func(sig_hist, lambda x, bc:  fomFuncWrapper( sig_hist.Integral(x,  sig_hist.GetNbinsX()+1 ) , bkg_hist.Integral(x,bkg_hist.GetNbinsX()+1 )))
    
    #c=ROOT.TCanvas()
    #hfom.Draw('hist')
    hfom_lt.SetTitle("FOM(<x)")
    hfom_lt.SetLineColor(ROOT.kViolet)

    hfom_gt.SetTitle("FOM(>x)")
    hfom_gt.SetLineColor(ROOT.kRed)
    
    
    if per:
        for h in [hfom_gt, hfom_lt]:
            reverse =  True if h==hfom_gt else False
            
            sug_x = getHistNearMax(h, per=per, reverse=reverse)
            sug = h.GetXaxis().GetBinCenter(sug_x  + (-1 if reverse else +1) )
            #sug = getQuantile(h, 0.90)[0]
            h.quant = sug
            h.line  = ROOT.TLine(sug, 0, sug, h.GetMaximum())
            h.line.SetLineColor(h.GetLineColor())
            h.line.SetLineStyle(h.GetLineStyle())
            h.line.SetLineWidth(h.GetLineWidth())
            
    return hfom_lt, hfom_gt


def makeHistFromList(lst, bins=None,name ="Histo", func=None):
    if not bins:
        bins = [len(lst),0,len(lst)]
    h = ROOT.TH1F(name,name,*bins)
    for ib,l in enumerate(lst,1):
        if func:
            l = func(l)
        if hasattr(l,"sigma"):
            h.SetBinContent(ib,l.val)
            h.SetBinError(ib,l.sigma)        
        if hasattr(l,"std_dev"):
            h.SetBinContent(ib,l.nominal_value)
            h.SetBinError(ib,l.std_dev)        
        else:
            h.SetBinContent(ib,l)
    return h


def makeHistFromDict(di , bins=None, name="Histo", bin_order=None,func=None):
    if bin_order:
        lst   = [ di.get(x,0) for x in bin_order ]
        labels = bin_order #[ x for x in bin_order if x in di]
    else:
        lst    = di.values()
        labels = di.keys()
    h = makeHistFromList(lst, bins, name, func)
    for ib, bin_label in enumerate(labels,1):
        h.GetXaxis().SetBinLabel( ib, str(bin_label))
    return h



def getDistributionFromList(l, nBins,  minVal=None, maxVal=None, name="dist", weighted=False):
    minVal = minVal if minVal else min(l)
    maxVal = maxVal if maxVal else max(l)
    h = ROOT.TH1D(name,name,nBins,minVal,maxVal)
    for v in l:
        try:
            v_ = getVal(v)
        except ValueError:
            v_ = v
        if weighted:
            sigma = getSigma(v)
            h.Fill( v_ , 1/sigma**2 )
        else:
            h.Fill( v_  )
    return h






def makeLegendFromHists( hists, name="Legend", loc=[], cols=1, opt='f', font=42):
    leg = ROOT.TLegend( *loc)
    leg.SetName(name)
    leg.SetFillColorAlpha(0,0.001)
    leg.SetBorderSize(0)
    leg.SetNColumns(cols)
    for h in hists:
        opt_ = getattr(h,"legOption", opt) 
        leg.AddEntry( h, "#font[%s]{%s}"%(font, h.GetTitle()), opt_)
    return leg 


#def makeLegend( data=None, mc_stack=None, sig_stack=None, leg_location=None , nBkgInLeg=None, legx=[0.75, 0.95 ], legy=[0.7, 0.87 ]):
#    mc_hists   = list( reversed( mc_stack.GetHists() ) ) if mc_stack else None
#    data_hists = [data] if data else []
#    sig_hists  = list( sig_stack.GetHists() ) if isinstance( sig_stack, ROOT.THStack) else [sig_stack]
#    hists  = mc_hists + sig_hists + data_hists
#    nhists = len(hists)
#    from math import ceil
#    nBkgInLeg = nBkgInLeg if nBkgInLeg else nhists
#    nhists_   = nBkgInLeg * int( ceil(1.*nhists/nBkgInLeg) )
#    legx=legx[:]
#    legs=[]
#    legy_dens = (legy[1]-legy[0]) / float(nBkgInLeg)
#    if mc_hists:
#        subBkgLists = [ hists[x:x+nBkgInLeg] for x in range(0, nhists , nBkgInLeg) ]
#        #print subBkgLists
#        nBkgLegs = len(subBkgLists)
#        for i , subBkgList in enumerate( subBkgLists ):
#            legloc = [legx[0], legy[1]- len(subBkgList)*legy_dens  ,legx[1],legy[1]]
#            bkgLeg = makeLegendFromHists( [h for h in subBkgList if h in mc_hists], name='leg_%s'%i , loc = legloc )# [legx[0], newLegY0 ,legx[1],legy[1]] )
#            #print "==========================================================================="
#            #print bkgLeg, subBkgList, "\n" , legloc #[legx[0], newLegY0  , legx[1],legy[1]]
#            #print "==========================================================================="
#            legs.append(bkgLeg)
#            #legx = [ 2*legx[0] -legx[1] , legx[0]  ]
#            dx   = legx[1]-legx[0]
#            legx = [ legx[0]-1.*dx , legx[1]-dx  ]
#    if sig_stack:
#        for sig in sig_hists:
#            bkgLeg.AddEntry(sig, sig.GetTitle(), 'l')
#    if data:
#        bkgLeg.AddEntry( data, data.GetTitle(), 'lp')
#    return legs


def makeLegend( data=None, mc_stack=None, sig_stack=None, leg_location=None , nBkgInLeg=None, legx=[0.75, 0.95 ], legy=[0.7, 0.87 ]):
    mc_hists   = list( reversed( mc_stack.GetHists() ) ) if mc_stack else None
    data_hists = [data] if data else []
    sig_hists  = list( sig_stack.GetHists() ) if isinstance( sig_stack, ROOT.THStack) else [sig_stack]
    #hists  = data_hists + mc_hists + sig_hists #+ data_hists
    hists  = mc_hists + data_hists + sig_hists #+ data_hists
    #if defaults.bigLeg:
    #    hists  = mc_hists + data_hists + sig_hists #+ data_hists
    #else:
    #    hists  = data_hists + mc_hists + sig_hists #+ data_hists
    nhists = len(hists)
    from math import ceil
    nBkgInLeg = nBkgInLeg if nBkgInLeg else nhists
    nhists_   = nBkgInLeg * int( ceil(1.*nhists/nBkgInLeg) )
    legx=legx[:]
    legs=[]
    legy_dens = (legy[1]-legy[0]) / float(nBkgInLeg)
    for h in mc_hists:
        h.legopt='f'
    for h in sig_hists:
        h.legopt='l'
    for h in data_hists:
        h.legopt='lpe'
    if mc_hists:
        subBkgLists = [ hists[x:x+nBkgInLeg] for x in range(0, nhists , nBkgInLeg) ]
        #print subBkgLists
        nBkgLegs = len(subBkgLists)
        for i , subBkgList in enumerate( subBkgLists ):
            legloc = [legx[0], legy[1]- len(subBkgList)*legy_dens  ,legx[1],legy[1]]
            #bkgLeg = makeLegendFromHists( [h for h in subBkgList if h in mc_hists], name='leg_%s'%i , loc = legloc )# [legx[0], newLegY0 ,legx[1],legy[1]] )
            bkgLeg = makeLegendFromHists( subBkgList, name='leg_%s'%i , loc = legloc )# [legx[0], newLegY0 ,legx[1],legy[1]] )
            #print "==========================================================================="
            #print bkgLeg, subBkgList, "\n" , legloc #[legx[0], newLegY0  , legx[1],legy[1]]
            #print "==========================================================================="
            legs.append(bkgLeg)
            #legx = [ 2*legx[0] -legx[1] , legx[0]  ]
            isLastCol = i == len(subBkgLists)-2
            dx   = legx[1]-legx[0]
            f_   = 0.08 if defaults.bigLeg else (isLastCol)*0.1
            legx = [ legx[0] - dx - f_ , legx[1] - dx   ]
    #if sig_stack:
    #    for sig in sig_hists:
    #        bkgLeg.AddEntry(sig, sig.GetTitle(), 'l')
    #if data:
    #    bkgLeg.AddEntry( data, data.GetTitle(), 'lp')
    return legs



def fixHistoStyle(h, x_title=None, y_title=None, min_max=None, **kwargs):
    """
        kwargs can be extra options, like kwargs={ 'lineColor':2, 'SetFillStyle':2}
    """
    if x_title:
        h.GetXaxis().SetTitle(x_title)
    if y_title:
        h.GetYaxis().SetTitle(y_title)
    if min_max:
        h.GetYaxis().SetRangeUser(*min_max)
    h.GetYaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(40)
    h.GetXaxis().SetTitleSize(40)
    h.GetXaxis().CenterTitle()
    h.GetYaxis().CenterTitle()
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(24)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(26)
    #h.GetXaxis().LabelsOption("d")
    h.GetXaxis().SetLabelOffset(0.01)
    h.GetXaxis().SetTitleOffset(1)
    h.GetYaxis().SetTitleOffset(1)
    for k,v in kwargs.items():
        root_set_attr_name = "Set"+k[0].upper()+k[1:]
        if hasattr(h,k):
            getattr(h,k)(v)
        elif hasattr(h, root_set_attr_name):
            getattr(h, root_set_attr_name)(v)

#def makeCanvasPads(    c1Name="canvas",  c1ww=defaults.canvas_width, c1wh=defaults.canvas_height,
#                       p1Name="pad1", p1M=defaults.pad1_loc , p1Gridx=False, p1Gridy=False,
#                       p2Name="pad2", p2M=defaults.pad2_loc, p2Gridx=False, p2Gridy=False,
#                       joinPads=True,
#                       func=None
#                    ):


#def drawCMSHeader( preliminary = "", lumi = 35.9, lxy = [0.16,0.915], rxy=[0.77,0.915], textR="%0.1f fb^{-1} (13 TeV)", cmsinside=True):
#    isPaper = preliminary.lower() in ['paper','']
#    latex = ROOT.TLatex()
#    latex.SetNDC()
#    latex.SetTextSize(0.04)
#    font=52
#    latex.SetTextFont(font)
#    #latexTextL = "#font[%s]{CMS %s}"%(font, preliminary)
#    #latexTextL = "CMS %s"%(preliminary)
#    cmstextsize = 0.08 if isPaper else 0.06
#    cmstext = "#font[61]{CMS}"
#    if not cmsinside:
#        latexTextL = cmstext
#        if preliminary:
#            latexTextL += "  #font[%s]{%s}"%(font,preliminary)
#        latex.DrawLatex(lxy[0],lxy[1],  latexTextL)
#    else:
#        textCMSlarge = ROOT.TLatex()
#        textCMSlarge.SetNDC()
#        textCMSlarge.SetTextSize(cmstextsize)
#        textCMSlarge.SetTextAlign(13)   
#        textCMSlarge.SetTextFont(42)
#        textCMSlarge.DrawLatex(0.20,0.85, cmstext)
#
#        if preliminary:
#            prelim = "#font[%s]{%s}"%(font,preliminary)
#            textPrelimlarge = ROOT.TLatex()
#            textPrelimlarge.SetNDC()
#            textPrelimlarge.SetTextSize(0.06*0.6)
#            textPrelimlarge.SetTextAlign(13)   
#            textPrelimlarge.SetTextFont(42)
#            textPrelimlarge.DrawLatex(0.21,0.78, prelim)
#
#       
#    if "%" in textR:
#        textR      = textR%lumi
#    latexTextR = "#font[42]{%s}"%(textR)
#    #latexTextR = "#font[%s]{%0.1f fb^{-1} (13 TeV)}"%(lumi)
#    latex.DrawLatex(rxy[0],rxy[1],  latexTextR)




def makeGraphFromDF(df, x=None, y=None, xe=None, ye=None, nPoints=None, name=None, title=None, xtitle=None, ytitle=None, markerColor=None, lineColor=None, markerStyle=20, lineStyle=None, verbose=False):
    if df.empty:
        print("df is empty")
        return ROOT.TGraphErrors()
    #import numpy as np
    arr     = np.array
    nPoints = nPoints if nPoints else len(df)
    zeros = arr([0]*nPoints)

    arr_x  = arr(getattr(df,x) )   if x  else arr(df.index)
    arr_y  = arr(getattr(df,y) )   if y  else zeros
    arr_xe = arr(getattr(df,xe))   if xe else zeros
    arr_ye = arr(getattr(df,ye))   if ye else zeros
    if verbose: print( arr_x, arr_y, arr_xe, arr_ye)
    if hasattr(arr_x[0], 'val') : # is a u_float!
        arr_x_ = arr([x_.val for x_ in arr_x])
        if not xe:
            arr_xe = arr([x_.sigma for x_ in arr_x])
        arr_x = arr_x_
    elif hasattr(arr_x[0], 'nominal_value') : # is a u_float!
        if verbose: print('x is ufloat')
        arr_x_ = arr([x_.nominal_value for x_ in arr_x])
        if not xe:
            arr_xe = arr([x_.std_dev for x_ in arr_x])
        arr_x = arr_x_
    if verbose: print(type( arr_y[0]), )
    if hasattr(arr_y[0], 'val') : # is a u_float!
        if verbose: print("y is u_float...")
        arr_y_ = arr([y_.val for y_ in arr_y])
        if not ye:
            if verbose: print("ye is u_float...")
            arr_ye = arr([y_.sigma for y_ in arr_y])
    elif hasattr(arr_y[0], 'nominal_value') : # is a ufloat!
        if verbose: print("y is ufloat...")
        arr_y_ = arr([y_.nominal_value for y_ in arr_y])
        if not ye:
            if verbose: print("ye is ufloat...")
            arr_ye = arr([y_.std_dev for y_ in arr_y])
        arr_y = arr_y_

    if verbose: print('n,x,y,xe,ye:',nPoints, arr_x, arr_y, arr_xe, arr_ye)
    graph = ROOT.TGraphErrors( nPoints, arr_x, arr_y, arr_xe, arr_ye )
    if verbose: graph.Print()


    if name:
        graph.SetName(name)
    if title:
        graph.SetTitle(title)
    if ytitle:
        graph.GetYaxis().SetTitle(ytitle)
    if xtitle:
        graph.GetXaxis().SetTitle(xtitle)
    if markerColor:
        graph.SetMarkerColor(markerColor)
    if lineColor:
        graph.SetLineColor(lineColor)
    if lineStyle:
        graph.SetLineStyle(lineStyle)
    if markerStyle:
        graph.SetMarkerStyle(markerStyle)

        
    return graph





########################################################################################################
###################     
########################################################################################################

def makeTGraphFromDF(df, x_col, y_col, title=None, x_title=None, y_title=None):
    """
        kinda dublicate of makeGraphFromDF but smarter.... need to combine the two 
    """
    from numpy import array
    nPoints = len(df)
    zeros   = array([0.0]*nPoints, dtype='float')
    x = df[x_col]
    y = df[y_col]
    #xe = array( x.apply( lambda v: getSigma(v)), dtype='float' )  if hasattr(x[0], 'sigma') or hasattr(x[0], 'std_dev') else zeros
    #ye = array( y.apply( lambda v: getSigma(v)), dtype='float' )  if hasattr(y[0], 'sigma') or hasattr(y[0], 'std_dev') else zeros
    xe = array( x.apply( lambda v: getSigma(v, strict=False, def_val=0)), dtype='float' )  
    ye = array( y.apply( lambda v: getSigma(v, strict=False, def_val=0)), dtype='float' )  
    withErrors = ( xe.any() or ye.any() ) # are all errors zero?
    #print(nPoints, ye,xe,withErrors)

    if withErrors:
        #print ( nPoints, array(x) , array(y.apply(lambda v: v.val)), xe, ye )
        #g = ROOT.TGraphErrors( nPoints, array(x,dtype='float') , array(y.apply(lambda v: getVal(v)), dtype='float'), xe, ye )
        g = ROOT.TGraphErrors( nPoints, array(x,dtype='float') , array(y.apply(lambda v: getVal(v, strict=False)), dtype='float'), xe, ye )
    else:
        g = ROOT.TGraphErrors( nPoints, array(x,dtype='float') , array(y,dtype='float'))
    g.GetXaxis().SetTitle(x_title if x_title else x_col)
    g.GetYaxis().SetTitle(y_title if y_title else y_col)
    if title:
        g.SetTitle(title)
    return g#,array(x), array(y.apply(lambda v: getVal(v))), xe,ye




def getDataMCRatios( data_hist, mc_hist ):
    import array as ar
    efill  = defaults.error_fill_style
    ecolor = 13#defaults.error_fill_color
     
    if isinstance( mc_hist, (tuple,list) ):
        mc_hist = getStackFromHists( mc_hist , "Stack" )

    if type(mc_hist) == ROOT.THStack :
        stack = mc_hist.Clone("stack")
        mc_hist = stack.GetHists().Last().Clone("mc_hist" )
        mc_hist.Reset()
        mc_hist.Merge( stack.GetHists() )
        

    unity = mc_hist.Clone( "IAmOne" )
    unity.SetLineColor(1)
    unity.SetLineWidth(1)
    unity.SetFillColor(0)
    nBins = unity.GetNbinsX()
    mc_noe = mc_hist.Clone( "mc_noerror" )
    #mc_noe.Sumw2(0)
    mc_noe.SetError(ar.array( "d",[0]*(nBins+1) ) ) 
    
    mc_e = mc_hist.Clone( "mc_error" )
    mc_e.Divide(mc_noe)
    mc_e.SetFillStyle(efill)
    #mc_e.SetLineColor( ecolor )
    mc_e.SetFillColor( ecolor -1)
    mc_e.SetMarkerSize(0)

    for ib in range( nBins+1 ):
        unity.SetBinContent(ib, 1)
        unity.SetBinError(ib, 0)

    data_ratio = data_hist.Clone( "data_ratio"  )
    data_ratio.Divide( mc_noe )

    mc_eb = mc_hist.Clone("mc_errorbar" )
    mc_eb.SetFillStyle( efill )
    mc_eb.SetLineColor( ecolor )
    mc_eb.SetMarkerSize(0)
    mc_eb.SetFillColor( ecolor )
    #mc_eb.SetFillColor(ROOT.kBlue-5)

    ROOT.gStyle.SetHatchesSpacing(0.01)
    ROOT.gStyle.SetHatchesLineWidth(1)

    ret = [data_ratio, mc_e, mc_eb , unity, mc_noe]
    for x in ret: x.SetDirectory(0)
    return ret



def drawNiceDataPlot( data_hist, mc_stack, sig_stack = None ,mc_total = None, options={} , saveDir = "./" , name = "plot", leg= None):
    """
           mc_total can be given in order to propegate errors fully, otherwise errors in mc_stack will be added in quad
    """

    canv = []
    ratios = [] 
    uq    = name+"_"
    #print uq

    if not mc_total:
        mc_total = mc_stack.GetHists()[0].Clone( "total"+"_"+uq)
        mc_total.Reset()
        mc_total.Merge( mc_stack.GetHists() )

    canv_hw=(800,800)
    canv  = makeCanvasPads() 
    #print canv
    canv[1].cd()
    setLogY = options.get('logy',1)
    canv[1].SetLogy( setLogY )

    ratios = getDataMCRatios( data_hist  , mc_total )
    data_ratio , mc_e, mc_eb, unity, mc_noe = ratios
    ymax = max( getHistMax( mc_noe )[1] , getHistMax( data_hist )[1] )
    #if sig_stack:
    #    ymin = min( [getHistMin( mc_stack.GetHists().First() )[1] , getHistMax( data_hist )[1] , ])
    #else:
    #    ymin = min( [getHistMin( mc_stack.GetHists().First() )[1] , getHistMax( data_hist )[1] , getHistMax( sig_stack.GetHists().First() )[1] ])
    ymin = options.get('ymin', defaults.ymin )
    extras = [mc_stack]
    #print '---------------', mc_stack
    #mc_stack.Print("all")
    #mc_eb.Draw("E2")

    # Recreating the stack here for some reason because ROOT segfaults if I use mc_stack ( no clue why! )
    stack = ROOT.THStack( mc_stack.GetTitle() + "2", mc_stack.GetName() )  
    for h in mc_stack.GetHists():
        h.SetDirectory(0)
        stack.Add(h)
        #print h.Draw("same")
    mc_stack = stack
    #
    ymax = options.get("ymax",ymax) 
    mc_stack.Draw("hist")
    ytitle = options.get( "ytitle", "Events")
    mc_stack.GetYaxis().SetTitle( ytitle )
    mc_stack.GetYaxis().SetTitleOffset(1.0)
    mc_stack.SetMaximum(ymax* ( 1.5 + 15*setLogY) )
    mc_stack.SetMinimum( ymin )
    mc_eb.Draw("E2same")
    #mc_e.Print("all")
    if sig_stack:
        dOptSig = "same hist"
        if isinstance( sig_stack, ROOT.THStack ):
            dOptSig += " nostack"
        sig_stack.Draw(dOptSig)
    #drawCMSHeader()
    data_hist.Draw("same")

    if leg:
        leg = [leg] if not type(leg) in [list, tuple] else leg
        for l in leg:
            l.Draw()

    ## draw ratio
    ytitle_r = options.get( "ytitle_r", "Data/#Sigma MC")
    #ytitle_size = options.get("ytitle_size", 0.12 )
    #ytitle_offset = options.get("ytitle_offset", 0.5)
    xtitle = options.get( "xtitle")
    canv[2].cd()
    unity.Draw() 
    unity.GetYaxis().SetTitle(       ytitle_r ) 
    #unity.GetYaxis().SetTitleSize(   ytitle_size )  
    #unity.GetYaxis().SetTitleOffset( ytitle_offset )
    #unity.GetYaxis().SetLabelSize( unity.GetYaxis().GetLabelSize()*2)
    if xtitle: unity.GetXaxis().SetTitle( xtitle )
    unity.SetNdivisions(505, "y")
    nBinsX = unity.GetNbinsX()
    xsize = canv_hw[0]/( nBinsX +1)/180. #180 scale is arbitrary (but emperical!)
    xsize = min([0.12, xsize])
    unity.GetXaxis().SetLabelSize( xsize )
    unity.GetXaxis().LabelsOption("v")
    mc_e.Draw("E2same")
    mc_e.Draw("E2same")
    data_ratio.Draw("E0p same")
    #data_ratio.SetMaximum(2)
    #data_ratio.SetMinimum(0)
    #degTools.saveCanvas( canv[0], saveDir , name)
    canv = list(canv)
    canv[1],canv[2] = adjustRatioPadTitleSizes( canv[1], canv[2] )
    canv = tuple(canv)

    defaults.CMS_lumi.CMS_lumi(canv[1], defaults.lumi, defaults.energy , defaults.iPosX) 
    canv[1].Update()
    canv[2].Update()
    if saveDir:
        if not os.path.isdir( saveDir ):
            os.makedirs( saveDir ) 
        save_path = os.path.join( saveDir, name )
        canv[0].SaveAs( save_path +".root")
        canv[0].SaveAs( save_path +".pdf")
        canv[0].SaveAs( save_path +".png")

        #canv[0].SetGrayscale(1)
        #canv[0].SaveAs( save_path +"_bw.pdf")
        #canv[0].SaveAs( save_path +"_bw.png")
        
    return canv, ratios, mc_stack



def getCanvPrims( canv , doClone=True):
    if doClone:
        ret= [canv.GetPrimitive( x.GetName() ).Clone() for x in canv.GetListOfPrimitives() ]
    else:
        ret= [canv.GetPrimitive( x.GetName() ) for x in canv.GetListOfPrimitives() ]
    for x in ret:
        if hasattr(x,"SetDirectory"): x.SetDirectory(0)
    return ret

def getCanvFirstHist( canv , doClone=False):
    prims = getCanvPrims(canv, doClone=doClone)
    hist_types = [ROOT.TH1D, ROOT.TH2D, ROOT.TH1F, ROOT.TH2F ]
    for p in prims:
        if any( [isinstance( p, htype) for htype in hist_types ] ):
            return p
    
    

def adjustRatioPadTitleSizes(pad1,pad2):
    h1=getCanvFirstHist( pad1, doClone=False)
    h2=getCanvFirstHist( pad2, doClone=False)
    pads_ratio = pad1.GetAbsHNDC() / pad2.GetAbsHNDC()
    #print pads_ratio 
    #print h1.GetYaxis().GetTitleOffset() , h2.GetYaxis().GetTitleOffset()

    h1.GetXaxis().SetTitleSize( getattr( defaults, "xtitle_size", 0.06 ) )
    h1.GetXaxis().SetLabelSize( getattr( defaults, "xtitle_size", 0.06 ) )
    h1.GetXaxis().SetTitleOffset( getattr( defaults, "xtitle_offset", 0.001 ) )

    h1.GetYaxis().SetTitleSize(   getattr( defaults, "ytitle_size", 0.06 ) )
    h1.GetYaxis().SetLabelSize(   getattr( defaults, "ytitle_size", 0.06 ) )
    h1.GetYaxis().SetTitleOffset( getattr( defaults, "ytitle_offset", 1.25 ) )

    xtitle_r_factor = getattr( defaults, "xtitle_r_factor", 0.8)
    h2.GetXaxis().SetTitleSize(  h1.GetXaxis().GetTitleSize()      * pads_ratio )
    h2.GetXaxis().SetLabelSize(  h1.GetXaxis().GetLabelSize()  * pads_ratio )
    h2.GetXaxis().SetTitleOffset(  h1.GetXaxis().GetTitleOffset()   )
    
    h2.GetYaxis().SetTitleSize(  h1.GetYaxis().GetTitleSize()      * pads_ratio )
    h2.GetYaxis().SetTitleOffset(  h1.GetYaxis().GetTitleOffset()  / pads_ratio * xtitle_r_factor)
    h2.GetYaxis().SetLabelSize(  h1.GetYaxis().GetLabelSize()  * pads_ratio * xtitle_r_factor)
    h2.GetYaxis().SetLabelOffset(0.01)
     
    h2.GetYaxis().SetRangeUser(*defaults.ratio_range)

    #print h1.GetXaxis().GetLabelOffset() , h2.GetXaxis().GetLabelOffset()
    print(h1.GetYaxis().GetLabelSize(), h2.GetYaxis().GetLabelSize())
    #print  h2.GetXaxis().GetTitleSize()

    return pad1, pad2

default_widths = {'y_width':500, 'x_width':500, 'y_ratio_width':200}
def makeRatioCanv( name="Canvas" , **widths):
        __doc__ = """ 
                default_widths = {default_widths}
        """.format(default_widths=default_widths)

        widths_ =deepcopy( default_widths )
        widths_.update(widths)
        widths = widths_
        widths['y_width'] += widths['y_ratio_width']
        scaleFacRatioPad = widths['y_width']/float( widths['y_ratio_width'] )
        y_border = widths['y_ratio_width']/float( widths['y_width'] )

        c1 = ROOT.TCanvas(str(uuid.uuid4()).replace('-','_'), "name",200,10, widths['x_width'], widths['y_width']     ) 
        c1.Divide(1,2,0,0)
        
        topPad = c1.cd(1)
        topPad.SetBottomMargin(0)
        topPad.SetLeftMargin(0.15)
        topPad.SetTopMargin(0.07)
        topPad.SetRightMargin(0.05)
        topPad.SetPad(topPad.GetX1(), y_border, topPad.GetX2(), topPad.GetY2())
        bottomPad = c1.cd(2)
        bottomPad.SetTopMargin(0)
        bottomPad.SetRightMargin(0.05)
        bottomPad.SetLeftMargin(0.15)
        bottomPad.SetBottomMargin(scaleFacRatioPad*0.13)
        bottomPad.SetPad(bottomPad.GetX1(), bottomPad.GetY1(), bottomPad.GetX2(), y_border)
        return (c1,topPad, bottomPad)




def th2Func2(hist, func = lambda x,y,bc: bc):
    """
        Returns a new histogram after applying func(xbin,ybin,bincontent) to each bin of hist
    """
    newhist = hist.Clone()
    newhist.Reset()
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    for x in range(nx):
        for y in range(ny):
            bc = hist.GetBinContent(x+1, y+1 )
            newbc = func(x,y,bc)
            newhist.SetBinContent(x+1, y+1, newbc)
    return newhist



def getTH2MaxBinContent(hist):
    bcs = getTH2FbinContent(hist)
    return max( itertools.chain( *[ y.values() for x,y in bcs.items() ] ) )
    


def getTH2DwithVarBins( c, var,  cutString = "(1)", weight = "weight"  , xbins=[0,2], ybins=[0,3], name = "testhist"):
    from array import array
    from PythonTools.NavidTools import uniqueName
    htmp = name +"_"+uniqueHash()
    print ( len(xbins)-1, array('d', xbins), len(ybins)-1, array('d', ybins) )
    h = ROOT.TH2D(htmp, htmp, len(xbins)-1, array('d', xbins), len(ybins)-1, array('d', ybins) )
    c.Draw(var+">>%s"%htmp, weight+"*("+cutString+")", 'goff')
    return h






def getHistBins(hist ):
    xbins = [hist.GetXaxis().GetBinLowEdge(ix+1) for ix in range(hist.GetNbinsX() +1 )  ]
    ybins = [hist.GetYaxis().GetBinLowEdge(iy+1) for iy in range(hist.GetNbinsY() +1 )  ]
    return xbins,ybins




def getTH2FbinContent(hist , legFunc= lambda xtitle,ytitle : (xtitle,ytitle), getError=False):
    """
       legFunc can be used to change the xtitle and ytitle in the output dictionary 
    """
    nbinsx = hist.GetNbinsX()
    nbinsy = hist.GetNbinsY()
    cont = {}
    for x in range(1,nbinsx+1):
        xbin = int( hist.GetXaxis().GetBinCenter(x) )
        #cont[xbin]={}
        for y in range(1,nbinsy+1):
            ybin = int( hist.GetYaxis().GetBinCenter(y) )
            bincontent = hist.GetBinContent(x,y)
            if bincontent:
                xtitle,ytitle = legFunc(xbin,ybin)
                if not cont.has_key(xtitle):
                    cont[xtitle]={}
                if getError:
                    binerror = hist.GetBinError(x,y)
                    ret      = u_float( bincontent, binerror )
                else:   
                    ret      = bincontent
                cont[xtitle][ytitle]=ret
    return cont







def getTH1FbinContent(hist, labeled=False, ordered = False, errors = True):
    '''
        returns bincontent of th1f with errors
        TODO: add option for over/under flow bins
    '''
    from collections import OrderedDict
    if type(hist) in [ROOT.THStack]:
        stack = hist.Clone()
        #hist  = getTotalFromStack(stack)
        hist  = getStackTot(stack)
         
    if errors:
        bin_values = [ u_float( hist.GetBinContent(ib) , hist.GetBinError(ib) )  for ib in range(1, hist.GetNbinsX() +1 ) ]
    else:
        bin_values = [ hist.GetBinContent(ib)  for ib in range(1, hist.GetNbinsX() +1 ) ]
    
    if not labeled:
        bin_labels = [ ib for ib in range(1, hist.GetNbinsX() +1 ) ]
        return bin_values
    else:
        bin_labels = [ hist.GetXaxis().GetBinLabel(ib) for ib in range(1, hist.GetNbinsX() +1 ) ]
        labels_values =  zip( bin_labels, bin_values) 
        d = OrderedDict if ordered else dict
        return d( labels_values ) 



#########################################################
###   efficiencies using RootTools.Sample
##########################################################


def getEfficiency(samples,samp, plot, cutInst_pass, cutInst_tot ,ret = False ):

    str_pass = cutInst_pass.fullName
    str_tot  = cutInst_tot.fullName

    try:
        h_pass = samples[samp]['cuts'][str_pass][plot]
        h_tot  = samples[samp]['cuts'][str_tot][plot]
    except KeyError:
        print ("!!!!!!!!!!!!!!!!!!!!!")
        print ("Plot key not for pass or tot not found.")
        print ("pass: samples[{samp}]['cuts'][{str_pass}][{plot}]".format(samp=samp, str_pass = str_pass, plot=plot))
        print ("tot:  samples[{samp}]['cuts'][{str_tot}][{plot}]".format(samp=samp, str_tot = str_tot, plot=plot))
        return False    

    #g_efficiency    =   ROOT.TGraphAsymmErrors()
    #g_efficiency.Divide(h_pass,h_tot,"cl=0.683 b(1,1) mode")

    h_eff = ROOT.TEfficiency(h_pass, h_tot)

    eff_name = 'EFF_%s_WRT_%s'%(str_pass,str_tot)
    eff_plot_name = plot + "_" + eff_name     

    #decorHist( samples[samp], cutInst_pass, h_eff, plots[plot]['decor'] ) 
    #decorHist( samples[samp], cutInst_pass, h_eff, {} ) 
    h_eff.SetName(samples[samp].name+"_"+eff_plot_name)
    h_eff.SetMarkerStyle(0)
    #h_eff.SetLineColor( samples[samp]['tree'].GetLineColor() )
    h_eff.SetLineColor( sample_colors[samp] ) 
    
    if samp in samples.bkgList():
        h_eff.SetLineWidth(2)
        h_eff.SetLineStyle(3)

    h_eff.SetTitle("{TITLE};{X};{Y}".format(TITLE=samples[samp].name+"_"+eff_plot_name,  X= h_pass.GetXaxis().GetTitle()  , Y= "#frac{%s}{%s}"%(str_pass, str_tot)  ))
    
    if not samples[samp]['cuts'].has_key(eff_name):
        samples[samp]['cuts'][eff_name] = {}
    samples[samp]['cuts'][eff_name][plot] = h_eff

    samples[samp]['plots'][eff_plot_name] = h_eff
    
    if ret:
        return h_eff







def getSampEfficiency2(samp, plot, cut_common="(1)", weight_num="(1)", weight_den="(1)", draw=False, 
                      total_title='total',
                      passed_title='passed',
                      title=None,
                      plot_name=None,
                     ):

    #exec(open('/afs/desy.de/user/n/nrad/analysis/commontools/PythonTools/ROOTHelpers.py').read())
    #from PythonTools import getAndSetEList
    getAndSetEList(samp.chain, cut_common, retrieve=True)
    h_den = getHistoFromSample(samp, var=plot.var, binning=plot.bins, cutString=str(cut_common) , weight=weight_den ,  addOverFlowBin=plot.overflow, combineWithSampleWeight=False )
    h_num = getHistoFromSample(samp, var=plot.var, binning=plot.bins, cutString=str(cut_common) , weight=weight_num ,  addOverFlowBin=plot.overflow, combineWithSampleWeight=False )
    samp.chain.SetEventList(0) ## reset eventlist
    h_den.Sumw2()
    h_num.Sumw2()
    
    total_title = total_title if total_title else weight_den
    passed_title = passed_title if passed_title else weight_num
    h_den.SetTitle(total_title)
    h_num.SetTitle(passed_title)
    
    h_den.SetLineColor(ROOT.kRed)
    h_num.SetLineColor(ROOT.kBlue)
    h_den.SetLineWidth(2)
    h_num.SetLineWidth(2)
    h_den.SetMarkerSize(1)
    h_num.SetMarkerSize(1)
    h_den.SetMarkerColor(ROOT.kRed)
    h_num.SetMarkerColor(ROOT.kBlue)
    trig_eff = ROOT.TEfficiency(h_num, h_den)
    #trig_eff.Draw()
    
    for h in [h_den, h_num]:
        h.GetXaxis().SetTitle(plot.xTitle)
        h.GetYaxis().SetTitle("Events")

    
    r=h_num.Clone()    
    r.SetLineColor(ROOT.kBlack)
    r.SetMarkerColor(ROOT.kBlack)
    r.Divide(h_den)
    #r.Draw()
    r.SetMinimum(0.2)
    r.SetMaximum(1.1)
    r.GetYaxis().SetTitle("Efficiency")
    #rr = r.Clone()
    
    rr = r.Clone()
    rr.SetLineColor(0)
    rr.SetMarkerColor(0)
    #rr.Reset()
    #rr.drawOption = 'hist'
    #rr.SetLineColor(0)
    #rr.SetMarkerColor(0)
    #r = trig_eff
    
    ret = {'ratio':r, 'eff':trig_eff, 'num':h_num, 'den':h_den}
    if draw:
        stack = getStackFromHists([h_den, h_num])
        stack.drawOption='nostack'
        if title:
            stack.SetTitle(title)
        r.SetTitle('')
        canvs = drawRatioPlot([stack], [r, trig_eff], 
                               widths={'x_width':1000, 'y_width':400, 'y_ratio_width':400}
                             )
        canvs['canv'].Draw()
        canvs['pads'][1].SetGridy(1)
        leg=canvs['pads'][0].BuildLegend()
        fixLegend(leg)
        canvs['stack']=stack
        ret.update(**canvs)
        ret['junk'] = [rr]
        if isinstance(draw, str):
            plot_name = plot_name if plot_name else p.name
            saveCanvas(canvs['canv'], draw, plot_name, formats=['png', 'pdf'], extraFormats=[])
    return ret




def getHistoEfficiency(h_den, h_num,
                       draw=False, 
                       total_title='total',
                       passed_title='passed',
                       title=None,
                       plot_name=None,
                       y_title="number of events",
                     ):
    #exec(open('/afs/desy.de/user/n/nrad/analysis/commontools/PythonTools/ROOTHelpers.py').read())
    #from PythonTools import getAndSetEList
    h_den.Sumw2()
    h_num.Sumw2()
    
    total_title  = total_title if total_title else "Total"
    passed_title = passed_title if passed_title else "Passed"
    h_den.SetTitle(total_title)
    h_num.SetTitle(passed_title)
    
    h_den.SetLineColor(ROOT.kRed)
    h_num.SetLineColor(ROOT.kBlue)
    h_den.SetLineWidth(2)
    h_num.SetLineWidth(2)
    h_den.SetMarkerSize(1)
    h_num.SetMarkerSize(1)
    h_den.SetMarkerColor(ROOT.kRed)
    h_num.SetMarkerColor(ROOT.kBlue)
    trig_eff = ROOT.TEfficiency(h_num, h_den)
    

    
    #r=h_num.Clone()    
    #r.SetLineColor(ROOT.kBlack)
    #r.SetMarkerColor(ROOT.kBlack)
    #r.Divide(h_den)
    r = getRatio(h_num, h_den, eff=True)

    #r.Draw()
    r.SetMinimum(0.2)
    r.SetMaximum(1.1)
    r.GetYaxis().SetTitle("Efficiency")
    #rr = r.Clone()
    
    #rr = r.Clone()
    #rr.SetLineColor(0)
    #rr.SetMarkerColor(0)
    
    ret = {'ratio':r, 'eff':trig_eff, 'num':h_num, 'den':h_den}
    if draw:
        stack = getStackFromHists([h_den, h_num])
        stack.drawOption='nostack'
        if title:
            stack.SetTitle(title)
        r.SetTitle('')
        
        canvs = drawRatioPlot([stack], [r, trig_eff], 
                               widths={'x_width':1000, 'y_width':400, 'y_ratio_width':400}
                             )
        #canvs['canv'].Draw()
        stack.GetYaxis().SetTitle(y_title)
        canvs['pads'][1].SetGridy(1)
        leg=canvs['pads'][0].BuildLegend()
        fixLegend(leg)
        canvs['stack']=stack
        ret.update(**canvs)
        #ret['junk'] = [rr]
        if isinstance(draw, str):
            plot_name = plot_name if plot_name else "eff"
            canvs['canv'].Modified()
            canvs['canv'].Update()
            canvs['canv'].cd()
            ret['title'] = drawLatex(title, x=0.55, y=0.97, font=42, align=22 )
            saveCanvas(canvs['canv'], draw, plot_name, formats=['png', 'pdf'], extraFormats=[])
            print('------------------- title:', title)
    return ret



def getSampEfficiency(samp, plot, cut_common="(1)", weight_num="(1)", weight_den="(1)", draw=False, 
                      total_title='total',
                      passed_title='passed',
                      title=None,
                      plot_name=None,
                     ):
    #exec(open('/afs/desy.de/user/n/nrad/analysis/commontools/PythonTools/ROOTHelpers.py').read())
    #from PythonTools import getAndSetEList
    getAndSetEList(samp.chain, cut_common, retrieve=True)
    h_den = getHistoFromSample(samp, var=plot.var, binning=plot.bins, cutString=str(cut_common) , weight=weight_den ,  addOverFlowBin=plot.overflow, combineWithSampleWeight=False )
    h_num = getHistoFromSample(samp, var=plot.var, binning=plot.bins, cutString=str(cut_common) , weight=weight_num ,  addOverFlowBin=plot.overflow, combineWithSampleWeight=False )
    samp.chain.SetEventList(0) ## reset eventlist
    h_den.Sumw2()
    h_num.Sumw2()
    
    total_title  = total_title if total_title else weight_den
    passed_title = passed_title if passed_title else weight_num
    plot_name    = plot_name if plot_name else p.name

    for h in [h_den, h_num]:
        h.GetXaxis().SetTitle(plot.xTitle)
        h.GetYaxis().SetTitle("Events")
    ret = getHistoEfficiency(h_den, h_num, total_title=total_title, passed_title=passed_title, title=title, draw=draw, plot_name=plot_name)
    
    return ret


def getSampsEfficiency(samps, plot, cut_common="(1)", weight_num="(1)", weight_den="(1)", draw=False, 
                      total_title='total',
                      passed_title='passed',
                      title=None,
                      plot_name=None,
                      combineWithSampleWeight=False,
                     ):
    """
        gets the efficiency for all samples
        combineWithSampleWeight should be False if samples include 
    """
    
    samps.getElists(cut_common)
    #getAndSetEList(samp.chain, cut_common, retrieve=True)
    histos_den = samps.getHistos([plot], cut=str(cut_common), weight=weight_den, rdf=False, combineWithSampleWeight=combineWithSampleWeight)
    histos_num = samps.getHistos([plot], cut=str(cut_common), weight=weight_num, rdf=False, combineWithSampleWeight=combineWithSampleWeight)

    hden = getSigBkgData(histos_den[plot.name])
    hnum = getSigBkgData(histos_num[plot.name])
    #hden_sig, hden_bkg, hden_data = getSigBkgData(histos_den[plot.name])
    #hnum_sig, hnum_bkg, hnum_data = getSigBkgData(histos_num[plot.name])
   
    hden_sig  = hden['sig']
    hden_bkg  = hden['bkg']
    hden_data = hden['data']

    hnum_sig  = hnum['sig']
    hnum_bkg  = hnum['bkg']
    hnum_data = hnum['data']
    
 
    hden_mc = addHists([hden_sig+hden_bkg], title="MC_den")
    hnum_mc = addHists([hnum_sig+hnum_bkg], title="MC_num")


    for h in [hden_mc, hnum_mc, hden_data, hnum_data, hden_sig, hnum_sig]:
        h.GetXaxis().SetTitle(plot.xTitle)
        h.GetYaxis().SetTitle("Events")
        h.SetFillColor(0) 
    
    total_title  = total_title if total_title else weight_den
    passed_title = passed_title if passed_title else weight_num
    plot_name    = plot_name if plot_name else p.name

    ret_mc    = getHistoEfficiency(hden_mc,   hnum_mc,   draw=draw, total_title=total_title, passed_title=passed_title, title=f"{title} (MC)",   plot_name=f"MC_{plot_name}")
    ret_data  = getHistoEfficiency(hden_data, hnum_data, draw=draw, total_title=total_title, passed_title=passed_title, title=f"{title} (Data)", plot_name=f"Data_{plot_name}")
    ret_sig   = getHistoEfficiency(hden_sig, hnum_sig, draw=draw, total_title=total_title, passed_title=passed_title,   title=f"{title} (Taupair)", plot_name=f"taupair_{plot_name}")
    
    #ret_mc['histos'] = histos_den
    #ret_data['histos'] = histos
    
    return {'data':ret_data, 'mc':ret_mc, 'histos_den':histos_den, 'histos_num':histos_num, 'taupair':ret_sig}


def getSigBkgData(hists, foms=False):
    """
        hists is a list of lists
        sig = hists[0][0]
        bkgs = hists[0][1:]
        data = hists[1][:]
    """
    h_sig = hists[0][0]
    h_bkgs = hists[0][1:]
    h_bkg  = addTHns(h_bkgs)
    h_sig.SetTitle('Signal')
    h_bkg.SetTitle('Background')
    h_mc_tot = addTHns([h_bkg,h_sig])
    ret = {'sig':h_sig, 'bkg':h_bkg, 'mc_tot':h_mc_tot}

    if len(hists)>1:
        h_datas = hists[1][:]
        h_data = addTHns(h_datas)
        h_ratio  = getRatio(h_data, h_mc_tot)
        h_ratio.SetTitle("data/mc")
        h_ratio.SetName("ratio")
        ret.update(data=h_data, ratio=h_ratio) 
    #h_bkg.SetFillColor(0)  
 
    if foms:
        fom    = th2Func(h_sig, lambda x,y,bc: fomFunc(bc, h_bkg.GetBinContent(x,y)))
        fom100 = th2Func(h_sig, lambda x,y,bc: fomFunc(bc, h_bkg.GetBinContent(x,y), 100))
        fom.SetTitle("FOM")
        fom100.SetTitle("FOM (100xBKG)")
        
        ret.update(fom=fom, fom100=fom100)
    

 
    return ret


#########################################################
###   Optimization Tools (FOM, etc)
##########################################################




def getTHn(samp, keys_bins):
    from PythonTools.SparseTools import SparseUp
    
    keys_bins = np.array(keys_bins)
    var_names = keys_bins.T[0]
    binnings  = keys_bins.T[1]
    nbins = np.vstack( binnings ).T[0].astype('int')
    #xmin  = np.vstack( binnings ).T[1]
    #xmax  = np.vstack( binnings ).T[2]
    
    
    #samp = samps.tau_3pi
    sup = SparseUp(samp.name, keys_bins )
    sup.get_arrays(samp.upr)
    print(samp.name)
    arr = sup.get_arrays( samp.upr, keys=['hie'])
    trig = arr['hie']==1
    sup.get_arrays(samp.upr, idx=trig)
    #sup.get_arrays(samp.upr, idx=trig)
    if 'exp' in samp.name or 'data' in samp.name:
        #arr = sup.get_arrays( samp.upr, keys=['__experiment__', 'psnm_32', 'psnm_27'])
        #exp = arr['__experiment__']
        #psnm32 = arr['psnm_32']
        #psnm27 = arr['psnm_27']
        #first  = np.logical_and(np.logical_or( exp==7, exp==8), psnm27)
        #second = np.logical_and(exp==10, psnm32)
        ## (((__experiment__==7||__experiment__==8)&&psnm_27)||(__experiment__==10&&psnm_32))
        #trig = np.logical_or(first,second)
        #arr = sup.get_arrays( samp.upr, keys=['hie'])
        #trig = arr['hie']==1
        thn = sup.get_thn() ## NO TRIG SELECTION FOR DATA!!!
        
        #print(f"{samp.name}: I have to ignore the trigger selelction")
    else:
        #return
        #sup.get_arrays
        if "hie" in samp.weightString:
            weight_string = samp.weightString.replace("* (hie)","")
        else:
            weight_string = samp.weightString
        thn = sup.get_thn(weight=eval(weight_string))
    return thn




def getTHnToBeFixed(samp, keys_bins, trigger_string=None):
    """

    """
    from PythonTools.SparseTools import SparseUp

    keys_bins = np.array(keys_bins)
    var_names = keys_bins.T[0]
    binnings  = keys_bins.T[1]
    nbins = np.vstack( binnings ).T[0].astype('int')
    #xmin  = np.vstack( binnings ).T[1]
    #xmax  = np.vstack( binnings ).T[2]


    #samp = samps.tau_3pi
    sup = SparseUp(samp.name, keys_bins )
    sup.get_arrays(samp.upr)
    print(samp.name)

    trigger_vars = getVariablesFromString(trigger_string)
    selection_arr = sup.get_arrays( samp.upr, keys=trigger_vars)
    trig = selection_arr['hie']==1
    sup.get_arrays(samp.upr, idx=trig)
    #sup.get_arrays(samp.upr, idx=trig)
    if 'exp' in samp.name or 'data' in samp.name:
        #arr = sup.get_arrays( samp.upr, keys=['__experiment__', 'psnm_32', 'psnm_27'])
        #exp = arr['__experiment__']
        #psnm32 = arr['psnm_32']
        #psnm27 = arr['psnm_27']
        #first  = np.logical_and(np.logical_or( exp==7, exp==8), psnm27)
        #second = np.logical_and(exp==10, psnm32)
        ## (((__experiment__==7||__experiment__==8)&&psnm_27)||(__experiment__==10&&psnm_32))
        #trig = np.logical_or(first,second)
        #arr = sup.get_arrays( samp.upr, keys=['hie'])
        #trig = arr['hie']==1
        thn = sup.get_thn() ## NO TRIG SELECTION FOR DATA!!!

        #print(f"{samp.name}: I have to ignore the trigger selelction")
    else:
        #return
        #sup.get_arrays
        if "hie" in samp.weightString:
            weight_string = samp.weightString.replace("* (hie)","")
        else:
            weight_string = samp.weightString
        thn = sup.get_thn(weight=eval(weight_string))
    return thn




def fomFunc(s,b, bkg_fact=1): 
    return s/math.sqrt(s + bkg_fact*b) if s else 0

def purityFunc(s, b):
    return s/(s+b) if s else 0


#def simpleFOM(s,b, bkg_fact=1):
#    return (s/sqrt(s+bkg_fact*b)) if s+b else 0

#bkg_fact=1
#h_sig_arr[10]=0
#h_bkg_arr[10]=0

def simpleFOM(s,b, bkg_fact=1):
    s = np.array(s,dtype='float64')
    b = np.array(b,dtype='float64')
    fom = np.divide(s, np.sqrt(s+bkg_fact*b), out=np.zeros_like(s), where= (s+b)!=0)
    return fom

def purity(s,b, bkg_fact=1):
    s = np.array(s, dtype='float64')
    b = np.array(b, dtype='float64')
    pur = np.divide(s, s+bkg_fact*b, out=np.zeros_like(s), where= (s+b)!=0)
    return pur

def getCutStringFromDict(di):
    """ gets a cut string from a dict of lower and upper value
        di= {
                'var1': (minVal, maxVal),
                'var2': (None, 1),
            }

    """
    cuts = []
    for k,v in di.items():
        cuts_ =[]
        lower_cut = v[0] #if  not v[0] is None  else False 
        upper_cut = v[1] #if  not v[1] is None  else False 
        
        #print(k, lower_cut, upper_cut, (lower_cut is False) ) 
        cuts_ += [f"{k}>=%s"%round(lower_cut,2)] if not lower_cut is None else []
        cuts_ += [f"{k}<=%s"%round(upper_cut,2)] if not upper_cut is None else []
        #print(cuts_)
        cuts.append(" && ".join(cuts_))
    cut = " && ".join(cuts)
    return cut


nSigAll = 769508.625 # n(tau pairs) for 8.8 /fb
def report(s_hist, b_hist, title="", n_sig_tot=nSigAll, verbose=True):
    
    s = s_hist.Integral()
    b = b_hist.Integral()
    fom = fomFunc(s, b)
    fom100 = fomFunc(s, b, 100)
    pur = purityFunc(s, b)
    eff = s/n_sig_tot
    if verbose: print(f"""{title}
    nSig:           {round(s,2)}
    nBkg:           {round(b,2)}
    Purity:         {round(pur*100,2)}%
    Signal Eff:     {round(eff*100,2)}%
    FOM:            {round(fom,2)}
    FOM(bkg x 100): {round(fom100,2)}
    
    """)
    ret = { 'nSig':     round(s,2),
            'nBkg':     round(b,2),
            'Purity':   round(pur*100,2),
            'eff':      round(eff*100,2),
            'fom':      round(fom,2),
            'fom100':   round(fom100,2),
          
          }
    return ret

def getCutPerformance(samps, cut="",pvar=None, nProc=1):
    """
         
    """
    #pvar   = plot_variables[pname]
    histos = samps.getHistos([pvar], cut=cut, rdf=False, nProc=nProc)
    hists = getSigBkgData(histos[pvar.name])
    h_sig = hists['sig']
    h_bkg = hists['bkg']
    return report(h_sig, h_bkg, title="", verbose=False)

#simpleFOM(h_sig_arr, h_bkg_arr)
maximizeFunc = lambda fom, pur, eff: fom*pur*pur

def findBestUpperLowerCut(h_sig_arr, h_bkg_arr, edges=[], maximize_func=maximizeFunc, bkg_fact=1, verbose=False):
    """
        calculates the best lower and upper cut which maximizes the maximize_func 
    """


    n_sig_tot  = h_sig_arr.sum()
    if not n_sig_tot:
        print('there is basically no signal here...')
    
    
    sig_cumsum = np.cumsum(h_sig_arr)
    bkg_cumsum = np.cumsum(h_bkg_arr)
    
    sig_cumsum_rev = np.flip( np.cumsum(np.flip(h_sig_arr) ) )
    bkg_cumsum_rev = np.flip( np.cumsum(np.flip(h_bkg_arr) ) )
    
    fom_upper_cut = simpleFOM(sig_cumsum, bkg_cumsum, bkg_fact=bkg_fact)
    fom_lower_cut = simpleFOM(sig_cumsum_rev, bkg_cumsum_rev, bkg_fact=bkg_fact)
    
    pur_upper_cut = purity(sig_cumsum, bkg_cumsum)
    pur_lower_cut = purity(sig_cumsum_rev, bkg_cumsum_rev)
    
    eff_upper_cut = sig_cumsum/n_sig_tot
    eff_lower_cut = sig_cumsum_rev/n_sig_tot
    
    upper_cut_to_maximize = maximize_func(fom_upper_cut, pur_upper_cut, eff_upper_cut)
    lower_cut_to_maximize = maximize_func(fom_lower_cut, pur_lower_cut, eff_lower_cut)
    
    best_upper_cut_bin = np.argmax(upper_cut_to_maximize) 
    best_lower_cut_bin = np.argmax(lower_cut_to_maximize)
    
    best_upper_cut = edges[best_upper_cut_bin] if sum(edges) else 999
    best_lower_cut = edges[best_lower_cut_bin] if sum(edges) else 999
    

    n_sig = round(h_sig_arr[best_lower_cut_bin:best_upper_cut_bin].sum(), 3)
    n_bkg = round(h_bkg_arr[best_lower_cut_bin:best_upper_cut_bin].sum(), 3)

    summary = {
                'best_upper_cut'    :best_upper_cut,
                'best_upper_cut_bin':best_upper_cut_bin,
                'best_lower_cut'    :best_lower_cut,
                'best_lower_cut_bin':best_lower_cut_bin,
                'n_sig': n_sig,
                'n_bkg': n_bkg,
                'eff'  : n_sig/n_sig_tot,
                'fom': simpleFOM(n_sig, n_bkg, bkg_fact=bkg_fact),
                'pur': purity(n_sig, n_bkg, bkg_fact=1.0),
                
    }
    
    
    if verbose:
        from pprint import pprint
        pprint(summary)
        print(f"""
          lower cut val: {round(best_lower_cut,3)},  ibin: {best_lower_cut_bin}, 
          upper cut val: {round(best_upper_cut,3)},  ibin: {best_upper_cut_bin},
         """)
    return summary




def findBestCuts(sig_thn, bkg_thn, var_names=None, maximize_func=lambda f,p,e:f, bkg_fact=1.0, reset_axes=True):
    from root_numpy import hist2array
    """
        loops over var_names and finds the best lower and upper cuts which maximizes maximize_func
        WARNING: sig_thn and bkg_thn need to have the same order in their axis!! FIXME: THERE IS NO CHECK FOR THIS!!! 

    """
    #thn = thns[0][0]
    nDim = sig_thn.GetNdimensions()
    
    
    axis_idx_dict = {sig_thn.GetAxis(iax).GetTitle():iax for iax in range(nDim)}
    
    if reset_axes:
        for iax in axis_idx_dict.values():
            for thn in [sig_thn, bkg_thn]:
                thn.GetAxis(iax).SetRange() ## Reset All Axes first!
            
    
    var_names = var_names if not var_names is None else list(axis_idx_dict)
    for v in var_names:
        if v not in axis_idx_dict:
            raise ValueError("The requested variable is not in the given THN")

    #for ivar, (vname, binning) in enumerate(var_names[:]):   
    results = {}
    cuts = {}
    for vname in var_names:
        #print(ivar, vname, binning)


        axis_idx = axis_idx_dict[vname]
        axis_sig = sig_thn.GetAxis(axis_idx)        
        axis_bkg = bkg_thn.GetAxis(axis_idx)

        sig_thn.GetAxis(axis_idx).SetRange() # make sure the axis is not limited before optimizing
        bkg_thn.GetAxis(axis_idx).SetRange() # make sure the axis is not limited before optimizing
        sig_thn_arr, edges = hist2array( sig_thn.Projection(axis_idx), return_edges=True )
        bkg_thn_arr, edges = hist2array( bkg_thn.Projection(axis_idx), return_edges=True )


        #ret = findBestUpperLowerCut(sig_thn_arr, bkg_thn_arr, edges[0], maximize_func=lambda f,p,e:f, bkg_fact=bkg_fact)
        ret = findBestUpperLowerCut(sig_thn_arr, bkg_thn_arr, edges[0], maximize_func=maximize_func, bkg_fact=bkg_fact)
        best_lower_cut_bin = int(ret['best_lower_cut_bin'])
        best_upper_cut_bin = int(ret['best_upper_cut_bin'])
        best_lower_cut = round(ret['best_lower_cut'],3)
        best_upper_cut = round(ret['best_upper_cut'],3)

        results[vname] = ret

        print(f"""\
        {vname}, iax: {axis_idx}, shape: {sig_thn_arr.shape},  
         best values : [{best_lower_cut}, {best_upper_cut}]
         bins  :       [{best_lower_cut_bin}, {best_upper_cut_bin}]
         eff:          {round(ret['eff'],2)}
        """)
        cuts[vname] = (best_lower_cut, best_upper_cut)
        sig_thn.GetAxis(axis_idx).SetRange(best_lower_cut_bin, best_upper_cut_bin) # limit to the best range before moving to the next variable
        bkg_thn.GetAxis(axis_idx).SetRange(best_lower_cut_bin, best_upper_cut_bin) # limit to the best range before moving to the next variable

        
    return results, cuts





#findBestUpperLowerCut(h_sig_arr, h_bkg_arr, edges=edges[0], maximize_func=lambda f,p,e: f*p, bkg_fact=100, verbose=True)







