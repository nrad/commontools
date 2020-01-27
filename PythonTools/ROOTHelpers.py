import ROOT
import itertools
import time
import uuid
import os

from PythonTools.NavidTools import uniqueName, hashString, hashObj


###################################
#      
###################################


def getPlotFromChain(c, var, binning, cutString = "(1)", weight = "weight", binningIsExplicit=False ,addOverFlowBin='',variableBinning=(False, 1) , name=None):
  """
          
        From HEPHYPythonTools! author unknown.
  """
  hname_tmp = uniqueName("h_tmp")

  if binningIsExplicit:
    h = ROOT.TH1D(hname_tmp, hname_tmp, len(binning)-1, array('d', binning))
#    h.SetBins(len(binning), array('d', binning))
  else:
    if len(binning)==6:
      h = ROOT.TH2D(hname_tmp, hname_tmp, *binning)
    else:
      h = ROOT.TH1D(hname_tmp, hname_tmp, *binning)

  c.Draw(var+">>%s"%hname_tmp, weight+"*("+cutString+")", 'goff')

  if variableBinning[0]:
    h.Sumw2()
    h.Scale(variableBinning[1],"width")

  res = h.Clone(name) if name else h.Clone()
  h.Delete()
  del h

  if addOverFlowBin.lower() == "upper" or addOverFlowBin.lower() == "both":
    nbins = res.GetNbinsX()
#    print "Adding", res.GetBinContent(nbins + 1), res.GetBinError(nbins + 1)
    res.SetBinContent(nbins , res.GetBinContent(nbins) + res.GetBinContent(nbins + 1))
    res.SetBinError(nbins , sqrt(res.GetBinError(nbins)**2 + res.GetBinError(nbins + 1)**2))
  if addOverFlowBin.lower() == "lower" or addOverFlowBin.lower() == "both":
    res.SetBinContent(1 , res.GetBinContent(0) + res.GetBinContent(1))
    res.SetBinError(1 , sqrt(res.GetBinError(0)**2 + res.GetBinError(1)**2))
  return res


def getEList(chain, cut, eListName=None, resetEList=True, keepAsChainAttr=False):
  if resetEList:
    chain.SetEventList(0) 
  #eListName_temp = uniqueName("eList_tmp")
  eListName_temp = getEListUniqueName(chain,cut,eListName)
  chain.Draw('>>%s'%eListName_temp, cut)
  #eListTemp = ROOT.gROOT.Get(eListName_temp)
  eListTemp = getattr(ROOT, eListName_temp)
  eList = eListTemp.Clone(eListName) if eListName else eListTemp.Clone()
  del eListTemp
  return eList.Clone()


def getAndSetEList(chain, cut, eListName=None, resetEList=True, retrieve=False):
    if retrieve:
        if not resetEList:
            raise Exception("Not safe to retrieve the EList if it is not being reset")
        eList = retrieveEList(chain,cut,eListName)
    else:
        eList = getEList(chain, cut, eListName=eListName, resetEList=resetEList) 

    chain.SetEventList(eList)
    return eList.Clone()

def getEListUniqueName(chain,cut,eListName=None):
    """
        Produces unique name for eList
        TODO: at each run the unique name is different because of the python object ID, maybe fix this?
    """
    eListName = eListName + "_" if eListName else "eList_"
    eListUniqueName = eListName + hashObj( [chain,cut] )
    return eListUniqueName

def retrieveEList(chain, cut, eListName=None):
    """
        tries to get EList from chain according the eListUniqueName, if not there, then creates eList
        TODO: keep a dictionary to check the checksum of eList itself after it is retrieved
    """
    eListUniqueName = getEListUniqueName(chain,cut,eListName) 
    eList = getattr(chain, eListUniqueName, False)
    if eList is False:
        print("Creating and Storing eList: %s"%eListUniqueName )
        eList = getEList(chain,cut,eListName,resetEList=True,keepAsChainAttr=True)
        setattr(chain, eListUniqueName, eList)
    return eList.Clone()


def getObjFromFile(fname, hname):
  f = ROOT.TFile(fname)
  assert not f.IsZombie()
  f.cd()
  htmp = f.Get(hname)
  if not htmp:  return htmp
  ROOT.gDirectory.cd('PyROOT:/')
  res = htmp.Clone()
  f.Close()
  return res

def th2Func(hist, func = lambda x,y,bc: bc, ignoreZeros=True):
    """
        Returns a new histogram after applying func(xbin,ybin,bincontent) to each bin of hist.
        if ignoreZeros is true, the new hist will not be filled if round(bincontent,10)==0
        Useful to perform operations on multiple Histograms:
        new_hist = th2Func(h1, func= lambda x,y,bc : h1.GetBinContent(x,y) + 0.5 h2.GetBinContent(x,y) )
        BUT YOU should make sure they  have the same binnings otherwise the result would be nonsense.
    """
    newhist = hist.Clone()
    newhist.Reset()
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    for x in range(nx+1):
        for y in range(ny+1):
            bc = hist.GetBinContent(x, y )
            newbc = func(x,y,bc)
            if ignoreZeros and round(newbc,10)==0:
               continue 
            newhist.SetBinContent(x, y, newbc)
    return newhist


def getTH2FbinContent(hist , legFunc= lambda x,y : (x,y), getError=False, ignoreZeros=True, binContentFunc=None):
    """
       returns a bin content as a dicitonary of dictionaries.
       legFunc can be used to change the xtitle and ytitle in the output dictionary 
       TODO: Return a numpy array
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
            #if bincontent or not ignoreZeros:
            if binContentFunc:
                bincontent=binContentFunc(bincontent)
            if bincontent or not ignoreZeros:
                xtitle,ytitle = legFunc(x,y)
                #if not cont.has_key(xtitle):
                if not xtitle in cont: 
                    cont[xtitle]={}
                if getError:
                    binerror = hist.GetBinError(x,y)
                    ret      = u_float( bincontent, binerror )
                else:
                    ret      = bincontent
                cont[xtitle][ytitle]=ret
    return cont


def getDistributionFromList(l, nBins,  minVal=None, maxVal=None, name="dist"):
    minVal = minVal if minVal else min(l)
    maxVal = maxVal if maxVal else max(l)
    h = ROOT.TH1D(name,name,nBins,minVal,maxVal)
    for v in l:
        v_ = getattr(v,"val",v)
        h.Fill( v_ )
    return h


def getTH2dist(th2d,*binning):
    """
        Returns a histogram with the distribution of the th2d
    """
    bc_dict = getTH2FbinContent(th2d)
    l = list(itertools.chain( * [[z for i,z in y.items()] for x,y in bc_dict.items()] ))
    binning = binning if binning else (100, min(l), max(l))
    h = ROOT.TH1D("dist","dist",*binning)
    for x in l:
        h.Fill(x)
    return h




def saveCanvas(canv,dir="./",name="",formats=["pdf","png"], extraFormats=["C","root"] , make_dir=True):
    if "$" in dir: 
        dir = os.path.expandvars(dir)
        if "$" in dir:
            raise Exception("Unresolved environmental variables! %s"%dir)
    if not os.path.isdir(dir) and make_dir: 
        makeDir(dir)
    if type(formats)!=type([]):
        formats = [formats]
    for form in formats:
        canv.SaveAs(dir+"/%s.%s"%(name,form) )
    if extraFormats:
        extraDir = dir+"/extras/"
        if not os.path.isdir(extraDir): mkdir_p(extraDir)
        for form in extraFormats:
            canv.SaveAs(extraDir+"/%s.%s"%(name,form) )


def drawLatex( txt = '', x=0.4, y=0.9,  font=22, size=0.04 , align=11, angle=0, setNDC=True):
    latex = ROOT.TLatex()
    if setNDC:
        latex.SetNDC()
    latex.SetTextSize(size)
    latex.SetTextFont(font)
    latex.SetTextAlign(align)
    latex.SetTextAngle(angle)
    if txt: latex.DrawLatex(x,y,  txt)
    return latex


def splitTree(tree, nSplit, branchesToKeep=None):
    nEntries = tree.GetEntries()
    
    if branchesToKeep:
        tree.SetBranchStatus("*",0) #deactivate branches
        for b in branchesToKeep:
            tree.SetBranchStatus(b,1)

    chunkFirstLast = getChunkFirstLast(nEntries, nSplit)
    print(f"Splitting tree ({nEntries} Entries) to {len(chunkFirstLast)}")
    for iTree, (first,last) in enumerate(chunkFirstLast):
        new = tree.CloneTree(0)
        new.SetName( "%s_%s"%(tree.GetName(),iTree) )
        new.SetTitle( "%s_%s"%(tree.GetTitle(),iTree) )
        for iEvt in range(first,last):
            tree.GetEntry(iEvt)
            new.Fill()
        yield new

    
def CanvasPartition(canv, nx, ny, lMargin=0.07, rMargin=0.05,
                      bMargin=0.07, tMargin=0.07):
    """
        Pythonized version of https://root.cern/doc/v610/canvas2_8C.html
    """
    
    # Setup Pad layout:
    vSpacing = 0.0
    vStep  = (1.- bMargin - tMargin - (ny-1) * vSpacing) / ny
    hSpacing = 0.0
    hStep  = (1.- lMargin - rMargin - (nx-1) * hSpacing) / nx;
    pads = []
    for i in range(nx):
        if i==0:
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr-hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        elif i == nx-1:
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = rMargin / (hposr-hposl);
        else:
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = 0.0;
        for j in range(ny):
            if j==0: 
               vposd = 0.0;
               vposu = bMargin + vStep;
               vfactor = vposu-vposd;
               vmard = bMargin / vfactor;
               vmaru = 0.0;
            elif j == ny-1:
               vposd = vposu + vSpacing;
               vposu = vposd + vStep + tMargin;
               vfactor = vposu-vposd;
               vmard = 0.0;
               vmaru = tMargin / (vposu-vposd);
            else:
               vposd = vposu + vSpacing;
               vposu = vposd + vStep;
               vfactor = vposu-vposd;
               vmard = 0.0;
               vmaru = 0.0;
            canv.cd(0);
            #char name[16];
            name = "pad_%i_%i"%(i,j)
            pad = ROOT.gROOT.FindObject(name)
            if pad:
                try:
                    pad.Delete()
                except:
                    print("failed to del canvas")
                    pass
                del(pad)
                
            pad = ROOT.TPad(name,"",hposl,vposd,hposr,vposu);
            pad.SetLeftMargin(hmarl);
            pad.SetRightMargin(hmarr);
            pad.SetBottomMargin(vmard);
            pad.SetTopMargin(vmaru);
            pad.SetFrameBorderMode(0);
            pad.SetBorderMode(0);
            pad.SetBorderSize(0);
            pad.Draw()
            pads.append(pad)
    return canv, pads


def removeBranchesFromTree(chain, branchesToKeep=[]):
    nEvents = chain.GetEntries()
    chain.SetBranchStatus("*",0)
    for br in branchesToKeep:
        chain.SetBranchStatus(br,1)
    new = chain.CloneTree(-1)
    return new

###################################
##
##   Non Root Helpers
##
###################################


def get_basename (f):
    return os.path.basename(f)
def get_filename (f):
    return os.path.splitext(os.path.basename(f))[0]
def get_ext (f):
    return os.path.splitext(os.path.basename(f))[1]


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        #if exc.errno == errno.EEXIST and os.path.isdir(path):
        if os.path.isdir(path):
            pass
        else:
            raise

def makeDir(path):
    if "." in path[-5:]:
        path = path.replace(os.path.basename(path),"")
        print (path)
    if os.path.isdir(path):
        return
    else:
        #mkdir_p(path)
        os.makedirs(path)



def getChunkFirstLast(nEntries,nSplit):
    """
    Creates list of first and last for splitting nEntries into nSplit chunks.
    to be used as  [first,last) , i.e range(first,last)
    follows the method of CloneTree as: https://root.cern.ch/root/html/tutorials/tree/copytree3.C.html
    """
    import math
    nPerTree = math.floor( nEntries/nSplit ) 
    chunkFirstLast = [ ( (i-1)*nPerTree,(i )*nPerTree  if i<nSplit else nEntries+1 )  for i in range(1,nSplit+1)  ]
    return chunkFirstLast


