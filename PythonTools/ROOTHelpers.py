import ROOT
import itertools
import functools
import time
import uuid
import os
import math
import glob
import shutil
import getpass


from PythonTools.NavidTools import uniqueName, hashString, hashObj
ROOT.TH1.SetDefaultSumw2()

user = getpass.getuser()
user_init = user[0]

###################################
#      
###################################


def joinStrings(*args, delim="&&", container="(%s)"):
    return (" %s "%delim).join([container%x for x in args])

joinCuts    = functools.partial(joinStrings, delim="&&")
joinWeights = functools.partial(joinStrings, delim="*")


def getVariablesFromString(var):
    return list(set([x for x in re.findall("\w+",cut) if not x[0].isdigit() ]))

###################################
#      
###################################

def getPlotFromChain(c, var, binning, cutString = "(1)", weight = "(1)", binningIsExplicit=False ,addOverFlowBin='',variableBinning=False , name=None, nEvents=None):
    """
                From HEPHYPythonTools!
                author unknown (probably Robert Schoefbeck)
    """
    hname_tmp = uniqueName("h_tmp")

    ndims = var.count(":")+1
    if ndims>3:
        raise NotImplemented("Trying to plot variable %s which has %s dimensions... but now only up to 3D is implemented... look into ThnSparse maybe?"%(var, ndims))
    TH = getattr(ROOT, "TH%sD"%ndims)

    if binningIsExplicit:
        from array import array
        if ndims==1:
            h = TH(hname_tmp, hname_tmp, 
                   len(binning)-1, array('d', binning))
        elif ndims==2:
            h = TH(hname_tmp, hname_tmp, 
                   len(binning[1])-1, array('d', binning[1]),
                   len(binning[0])-1, array('d', binning[0]))
        elif ndims==1:
            h = TH(hname_tmp, hname_tmp, 
                   len(binning[2])-1, array('d', binning[2]),
                   len(binning[1])-1, array('d', binning[1]),
                   len(binning[0])-1, array('d', binning[0]))
#        h.SetBins(len(binning), array('d', binning))
    else:
        if len(binning) in [3,6,9]:
            h = TH(hname_tmp, hname_tmp, *binning)
        else:
            raise ValueError("Can't make sense of the bins! %s"%binning)
    nEventsArgs = (nEvents,0) if nEvents else ()
    #c.Draw(var+">>%s"%hname_tmp, weight+"*("+cutString+")", 'goff', *nEventsArgs)
    c.Draw(var+">>%s"%hname_tmp, "(%s)*(%s)"%(cutString, weight), 'goff', *nEventsArgs)

    if variableBinning != False:
        h.Sumw2()
        h.Scale(variableBinning,"width")

    res = h.Clone(name) if name else h.Clone()
    h.Delete()
    del h

    if addOverFlowBin.lower() == "upper" or addOverFlowBin.lower() == "both":
        nbins = res.GetNbinsX()
#        print "Adding", res.GetBinContent(nbins + 1), res.GetBinError(nbins + 1)
        res.SetBinContent(nbins , res.GetBinContent(nbins) + res.GetBinContent(nbins + 1))
        res.SetBinError(nbins , math.sqrt(res.GetBinError(nbins)**2 + res.GetBinError(nbins + 1)**2))
    if addOverFlowBin.lower() == "lower" or addOverFlowBin.lower() == "both":
        res.SetBinContent(1 , res.GetBinContent(0) + res.GetBinContent(1))
        res.SetBinError(1 , math.sqrt(res.GetBinError(0)**2 + res.GetBinError(1)**2))
    return res


def getPlotFromChainOLD(c, var, binning, cutString = "(1)", weight = "(1)", binningIsExplicit=False ,addOverFlowBin='',variableBinning=(False, 1) , name=None, nEvents=None):
    """
                From HEPHYPythonTools!
                author unknown (probably Robert Schoefbeck)
    """
    hname_tmp = uniqueName("h_tmp")
    
    ndims = var.count(":")+1
    if ndims>3:
        raise NotImplemented("Trying to plot variable %s which has %s dimensions... but now only up to 3D is implemented... look into ThnSparse maybe?"%(var, ndims))
    TH = getattr(ROOT, "TH%s"%ndims)

    if binningIsExplicit:
        if ndims==1:
            h = TH(hname_tmp, hname_tmp, len(binning)-1, array('d', binning))
        if ndims==2:
            h
#        h.SetBins(len(binning), array('d', binning))
    else:
        if len(binning)==6:
            h = TH(hname_tmp, hname_tmp, *binning)
        elif len(binning)==3:
            h = TH(hname_tmp, hname_tmp, *binning)
        elif len(binning)==9:
            h = TH(hname_tmp, hname_tmp, *binning)
        else:
            raise ValueError("Can't make sense of the bins! %s"%binning)
    nEventsArgs = (nEvents,0) if nEvents else () 
    #c.Draw(var+">>%s"%hname_tmp, weight+"*("+cutString+")", 'goff', *nEventsArgs)
    c.Draw(var+">>%s"%hname_tmp, "(%s)*(%s)"%(cutString, weight), 'goff', *nEventsArgs)

    if variableBinning[0]:
        h.Sumw2()
        h.Scale(variableBinning[1],"width")

    res = h.Clone(name) if name else h.Clone()
    h.Delete()
    del h

    if addOverFlowBin.lower() == "upper" or addOverFlowBin.lower() == "both":
        nbins = res.GetNbinsX()
#        print "Adding", res.GetBinContent(nbins + 1), res.GetBinError(nbins + 1)
        res.SetBinContent(nbins , res.GetBinContent(nbins) + res.GetBinContent(nbins + 1))
        res.SetBinError(nbins , math.sqrt(res.GetBinError(nbins)**2 + res.GetBinError(nbins + 1)**2))
    if addOverFlowBin.lower() == "lower" or addOverFlowBin.lower() == "both":
        res.SetBinContent(1 , res.GetBinContent(0) + res.GetBinContent(1))
        res.SetBinError(1 , math.sqrt(res.GetBinError(0)**2 + res.GetBinError(1)**2))
    return res

def getChain(files, tree_name='tau3x1'):
    """
        gets TChain named "tree_name" from from files
        files can be a list of files, or str with file patern 
    """
    chain = ROOT.TChain(tree_name)
    
    if isinstance(files,str):
        files = glob.glob(files)
    if not len(files):
        print("Warning no files empty! %s"%files)
    for f in files:
        chain.Add(f)
    return chain




opt_dict = {"*":"x", "+":"plus", "/":"over", "-":"minus", "(":"op", ")":"cp", ".": "p", ":":"cln", ",":"com", " ":"", "||":"or"}
def defineVariableStringForRDF(v, prefix="defined_", opt_dict=opt_dict, ):
    vname = v
    redefined = False
    for opt,optname in opt_dict.items():
        if opt in v:
            vname = vname.replace(opt, f"__{optname}__")
            redefined = True
    if redefined:
        return f"{prefix}{vname}"


def defineVariablesForRDF(rdf, variable_list, prefix="defined_", opt_dict=opt_dict, ):
    new_vars_dict = {}
    for v in variable_list:
        redef = defineVariableStringForRDF(v)
        if redef:
            try:
                #print(f"define new variable, {redef}, {v}")
                rdf = rdf.Define(redef,v)
            except TypeError as exp:
                print ("Tried to Define Variable For RDF:", v, redef)
                raise exp
            new_vars_dict[v]=redef 
        else:
            new_vars_dict[v]=v
    return rdf, new_vars_dict





def getHistoFromRDF(rdf, var, binning, cut=None, weight=None, name=None, title=None): #, xtitle=None, ytitle=None):
    """
        Helper function to get histograms from RDataFrames
        'var' can be either variable string or [xvar,yvar]
             if there are formulas in the var (or weight), I will try to Define them for the RDF.
        'cut' can be applied for better performance (especially if you have multiple histograms to make) 
            apply Filter to RDF before passing it on to this function.
    """
    allvars = []

    if cut:
        try:
            rdf = rdf.Filter(str(cut))
        except Exception as exp:
            print(f"FAILED TO FILTER RDF WITH THE FOLLOWING SELECITON: \n {cut}")
            raise exp
        #allvars.append(cut)
    if weight:
        allvars.append(str(weight))

    if isinstance(var, (list,tuple)):
        nHistoDimension = len(var)
        var = tuple(str(x) for x in var)
        histoFuncString = "Histo2D" ## Howabout dimensional Histos?
        if nHistoDimension != 2:
            raise NotImplementedError("O.K... only dimensions up to Histo2D are implemented... fix it!")
        assert len(binning)==6, f"Binning doesn't match for Histo2D: {binning}"
    else:
        nHistoDimension = 1 
        var = (str(var),)
        histoFuncString = "Histo1D"
        assert len(binning)==3, f"Binning doesn't match for Histo1D: {binning}"

    allvars.extend(var) 
    rdf, varsdict = defineVariablesForRDF(rdf, allvars)
    newvars = tuple(new for old,new in varsdict.items() if old in var)


    name = name if name else uniqueName("Hist")
    title = title if title else name
    #histo2d_args = ( (name,title) +binning ,) + tuple( var_.split(":") )
    histo_args = ( (name,title) +binning ,) + newvars

    if weight:
        histo_args = histo_args + (varsdict[weight],)
    try:
        h = getattr(rdf, histoFuncString)(*histo_args)
    except TypeError as exp:
        print("FAILED MAKING RDF HISTOGRAM WITH ARGS: \n", histo_args)
        raise exp
    return h



def getRDFHistoFromSample(sample, *args, **kwargs): #, xtitle=None, ytitle=None):
    """
    wrapper for getHistoFromRDF
    Sample is expected to be a RootTools Sample
    args = cut, weight
    """
    if not hasattr(sample, "combineWithSampleWeight"):
        raise Exception(f"sample ({sample}) doesn't seem to be an instance of RootTools.Sample") 
    cut = kwargs.get("cut",None)
    weight = kwargs.get("weight",None)
    cut = str(cut) if cut else cut
    weight = str(weight) if weight else weight
    kwargs['cut']    = sample.combineWithSampleSelection( cut )   if  kwargs.pop("combineWithSampleSelection", True) else cut 
    kwargs["weight"] = sample.combineWithSampleWeight( weight )   if  kwargs.pop("combineWithSampleWeight", True) else weight
    kwargs["title"]  = kwargs.get("title", sample.name)
    #print('getting histo from rdf', sample.name, args, kwargs)
    histo = getHistoFromRDF(sample.rdf, *args, **kwargs)
    histo.cut    = kwargs['cut']
    histo.weight = kwargs['weight']
    #print('                  ....... done')
    return histo
     
def getHistoFromSample(sample, *args, **kwargs): #, xtitle=None, ytitle=None):
    """
    wrapper for getPlotFromChain
    Sample is expected to be a RootTools Sample
    args = cut, weight
    """
    if not hasattr(sample, "combineWithSampleWeight"):
        raise Exception(f"sample ({sample}) doesn't seem to be an instance of RootTools.Sample") 
    cut = kwargs.get("cutString",None)
    weight = kwargs.get("weight",None)
    cut    = str(cut)    if cut else cut
    weight = str(weight) if weight else weight

    kwargs['cutString']    = sample.combineWithSampleSelection( cut )   if  kwargs.pop("combineWithSampleSelection", True) else cut 
    kwargs["weight"]       = sample.combineWithSampleWeight( weight )   if  kwargs.pop("combineWithSampleWeight", True) else weight
  
    #if 'nFractEvents' in kwargs:
    #    kwargs['nEvents'] = int(sample.chain.GetEntries()*kwargs['nFractEvents'])
    #    kwargs["weight"] = f"({kwargs['weight']}*{kwargs['nFractEvents']})"
    histo = getPlotFromChain(sample.chain, *args, **kwargs)
    histo.SetName(sample.name)
    histo.SetTitle(sample.texName) 
    histo.cut = kwargs['cutString']
    histo.weight = kwargs['weight']
        

    return histo
     


def getGraphsFromContours(contours):
    """
        h.Draw("CONT LIST")
        contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
        gs=getGraphsFromContours(contours)

    """
    gs = []
    for igl in range(contours.GetSize()):
        gl = contours.At(igl)
        for ig in range(gl.GetSize()):
            g = gl.At(ig)
            gs.append(g)
            print(g, g.Integral())
    gs = sorted(gs, key=lambda x: x.Integral(), reverse=True)
    return gs





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
    #eListName = eListName + "_" if eListName else "eList_"
    eListName = eListName + "_" if eListName else "%s_"%chain.GetName()
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


#
def th2Func(hist, func = lambda x,y,bc: bc, ignoreZeros=True, name=None):
    """
        Returns a new histogram after applying func(xbin,ybin,bincontent) to each bin of hist.
        if ignoreZeros is true, the new hist will not be filled if round(bincontent,10)==0
        Useful to perform operations on multiple Histograms:
        new_hist = th2Func(h1, func= lambda x,y,bc : h1.GetBinContent(x,y) + 0.5 h2.GetBinContent(x,y) )
        BUT YOU should make sure they  have the same binnings otherwise the result would be nonsense.
    """
    newhist = hist.Clone()
    if name:
        newhist.SetTitle(name)
        newhist.SetName(name)
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


def th1Func(hist, func = lambda x,y,bc: bc, ignoreZeros=True):
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
    for x in range(nx+1):
        bc = hist.GetBinContent(x)
        newbc = func(x,bc)
        if ignoreZeros and round(newbc,10)==0:
           continue 
        newhist.SetBinContent(x,newbc)
    return newhist

def th1FuncError(hist, func = lambda x,y,bc: bc, ignoreZeros=True):
    """
        same as th1Func but includes errors
        Returns a new histogram after applying func(xbin,ybin,bincontent) to each bin of hist.
        if ignoreZeros is true, the new hist will not be filled if round(bincontent,10)==0
        Useful to perform operations on multiple Histograms:
        new_hist = th2Func(h1, func= lambda x,y,bc : h1.GetBinContent(x,y) + 0.5 h2.GetBinContent(x,y) )
        BUT YOU should make sure they  have the same binnings otherwise the result would be nonsense.
    """
    from uncertainties import ufloat
    newhist = hist.Clone()
    newhist.Reset()
    nx = hist.GetNbinsX()
    for x in range(nx+1):
        bc = ufloat(hist.GetBinContent(x),hist.GetBinError(x))
        newbc = func(x,bc)
        if ignoreZeros and round(getVal(newbc, strict=False),10)==0 and round(getSigma(newbc, strict=False),10)==0:
            continue 
        newhist.SetBinContent(x,getVal(newbc, strict=False) )
        newhist.SetBinError(x, getSigma(newbc, strict=False, def_val=0))
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

def getTF1Params(f, errors=True, di=True):

    """
        returns a dictionary or list with the parameter values of the function
    """

    nPars = f.GetNpar()
    if errors:
        from uncertainties import ufloat
        vals =  [ufloat(f.GetParameter(i), f.GetParError(i)) for i in range(nPars)]
    else:
        vals = [f.GetParameter(i) for i in range(nPars)]
    if di:
        names = [f.GetParName(i) for i in range(nPars)]
        vals  = dict(zip(names, vals))
    return vals



def copyIndexPHP(directory):
    ''' Copy index.php to directory
    '''
    index_php_to = os.path.join( directory, 'index.php' )
    os.makedirs(directory, exist_ok=True) 
    index_php = f'/afs/desy.de/user/{user_init}/{user}/www/index.php'
    if not os.path.isfile(index_php):
        raise Exception("unable to find the original PHP file...what have you done with it?!! %s"%index_php)
    try:
        shutil.copyfile( index_php, index_php_to )
    except Exception as exp:
        print("Failed to copy index.php... looked for it in:\n%s \ntried to copy it to:%s\n"%(index_php,index_php_to))
        raise exp

#def saveCanvas(canv,dir="./",name="",formats=["pdf","png"], extraFormats=["C","root"] , make_dir=True):
def saveCanvas(canv, dir="./",name="", formats=["png"], extraFormats=[] , make_dir=True, showWebLink=False, copyPHP=True):
    if "$" in dir: 
        dir = os.path.expandvars(dir)
        if "$" in dir:
            raise Exception("Unresolved environmental variables! %s"%dir)
    if not os.path.isdir(dir) and make_dir:
        makeDir(dir)
    copyIndexPHP(dir)
    if type(formats)!=type([]):
        formats = [formats]
    for form in formats:
        canv.SaveAs(dir+"/%s.%s"%(name,form) )
    if extraFormats:
        extraDir = dir+"/extras/"
        if not os.path.isdir(extraDir): mkdir_p(extraDir)
        for form in extraFormats:
            canv.SaveAs(extraDir+"/%s.%s"%(name,form) )
    if showWebLink and 'www' in dir:
        print( f'https://www.desy.de/~{user}/' + dir.split('www')[-1] + ("/index.php" if copyPHP else "") )

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

def makeCanvasPads(    c1Name="canvas",  c1ww=800, c1wh=650,
                       p1Name="pad1", p1M=(0.0, 0.4, 1, 1) , p1Gridx=False, p1Gridy=False,
                       p2Name="pad2", p2M=(0.0, 0,1,0.4), p2Gridx=False, p2Gridy=False,
                       joinPads=True,
                       func=None
                    ):
  c = ROOT.TCanvas(c1Name,c1Name,c1ww,c1wh)

  pad1 = ROOT.TPad(p1Name, p1Name, *p1M)
  pad1.SetBottomMargin(0)  # joins upper and lower plot
  if p1Gridx: pad1.SetGridx()
  if p1Gridy: pad1.SetGridy()

  # Lower ratio plot is pad2
  c.cd()  # returns to main canvas before defining pad2
  pad2 = ROOT.TPad(p2Name, p2Name, *p2M)

  if joinPads: pad2.SetTopMargin(0)  # joins upper and lower plot
  pad2.SetBottomMargin(0.3)
  if p2Gridx: pad2.SetGridx()
  if p2Gridy: pad2.SetGridy()
  if func:
    func(pad1,pad2)
  pad1.Draw()
  pad2.Draw()
  return c, pad1, pad2


def removeBranchesFromTree(chain, branchesToKeep=[]):
    nEvents = chain.GetEntries()
    chain.SetBranchStatus("*",0)
    for br in branchesToKeep:
        chain.SetBranchStatus(br,1)
    new = chain.CloneTree(-1)
    return new

###






    
def adjustRatioStyle2(r, plot=None, ref=None, size_factor = 1):
    print(r,plot,ref, size_factor)
    r.GetXaxis().SetLabelSize(0.12  if not ref else size_factor*ref.GetXaxis().GetLabelSize() )
    r.GetYaxis().SetLabelSize(0.12  if not ref else size_factor*ref.GetYaxis().GetLabelSize() )
    r.GetXaxis().SetTitleSize(0.1   if not ref else size_factor*ref.GetXaxis().GetTitleSize() )
    r.GetYaxis().SetTitleSize(0.1   if not ref else size_factor*ref.GetYaxis().GetTitleSize() )
    r.GetXaxis().SetTitleOffset(1   if not ref else ref.GetXaxis().GetTitleOffset() )
    r.GetYaxis().SetTitleOffset(0.5 if not ref else ref.GetYaxis().GetTitleOffset()/size_factor )
    r.GetYaxis().SetNdivisions(5,4,2)
    
    if plot:
        #r.GetYaxis().SetTitle("new/old   ")
        r.GetXaxis().SetTitle( plot.xTitle )
    elif ref:
        r.GetXaxis().SetTitle( ref.GetXaxis().GetTitle() )

def adjustHistAxis(h):
    
    for ax in [h.GetXaxis(), h.GetYaxis()]:
        ax.SetLabelFont(43)
        ax.SetTitleFont(43)
        ax.SetLabelSize(26)
        ax.SetTitleSize(26)
    #h.GetXaxis().SetTitleOffset(3)
    #h.GetYaxis().SetTitleOffset(1.5)

    #if plot:
    #    #r.GetYaxis().SetTitle("new/old   ")
    #    r.GetXaxis().SetTitle( plot.xTitle )
    #elif ref:
    #    r.GetXaxis().SetTitle( ref.GetXaxis().GetTitle() )

        
        
def createCanvasPads(cname="drawHistos", wtopx=200, wtopy=10, xwidth=500, ywidth=500, ratiowidth=200):
    c1 = ROOT.TCanvas(str(uuid.uuid4()).replace('-','_'), cname,wtopx,wtopy, xwidth, ywidth)
    y_border = ratiowidth/(ywidth+ratiowidth)
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
    #bottomPad.SetBottomMargin(0.3)
    bottomPad.SetBottomMargin(ywidth/ratiowidth*0.16)
    bottomPad.SetPad(bottomPad.GetX1(), bottomPad.GetY1(), bottomPad.GetX2(), y_border)
    return c1, topPad, bottomPad

def drawRatioPlot(  topPadObjects, 
                    bottomPadObjects=None, 
                    widths={}, 
                    plot=None, 
                    doDraw=False, 
                    y_title_offset=2, 
                    x_title_offset=4.2, 
                    y_label_offset=0.02,
                    fix_overlapping_labels=False,
                ):
    """
    
    
    """
    ratio = True if bottomPadObjects else False
    bottomPadObjects = bottomPadObjects if isinstance(bottomPadObjects, (list,tuple)) else [bottomPadObjects]
    topPadObjects = topPadObjects if isinstance(topPadObjects, (list,tuple)) else [topPadObjects]
    
    #y_title_offset = 2
    default_widths = {'x_width':1000, 'y_width':400, 'y_ratio_width':400}
    default_widths.update( widths )
    #y_border = default_widths['y_ratio_width']/float( default_widths['y_width'] )
    #default_widths['y_width'] += default_widths['y_ratio_width']

    if ratio is not None:
        scaleFacRatioPad = default_widths['y_width']/float( default_widths['y_ratio_width'] )
        default_widths['y_width'] += default_widths['y_ratio_width']
        #y_border         = default_widths['y_ratio_width']/float( default_widths['y_width'] )

    widths = default_widths
    #c1 = ROOT.TCanvas(str(uuid.uuid4()).replace('-','_'), "drawHistos",200,10, default_widths['x_width'], default_widths['y_width'])
    c1, topPad, bottomPad = createCanvasPads('canvname', xwidth=widths['x_width'], ywidth=widths['y_width'] , ratiowidth=widths['y_ratio_width'] )
    #c1 = ROOT.TCanvas(str(uuid.uuid4()).replace('-','_'), "drawHistos", default_widths['x_width'], default_widths['y_width'])
    #c1.Divide(1,2,0,0)
    #topPad = c1.cd(1)
    #topPad.SetBottomMargin(0)
    #topPad.SetLeftMargin(0.15)
    #topPad.SetTopMargin(0.07)
    #topPad.SetRightMargin(0.05)
    #topPad.SetPad(topPad.GetX1(), y_border, topPad.GetX2(), topPad.GetY2())
    #bottomPad = c1.cd(2)
    #bottomPad.SetTopMargin(0.005)
    #bottomPad.SetRightMargin(0.05)
    #bottomPad.SetLeftMargin(0.15)
    #bottomPad.SetBottomMargin(scaleFacRatioPad*0.13)
    #bottomPad.SetPad(bottomPad.GetX1(), bottomPad.GetY1(), bottomPad.GetX2(), y_border)

    if plot:
        topPad.SetLogy( plot.logY )
    
    c1.cd()
    topPad.cd()
    topRefHist = topPadObjects[0]
    #topRefHist.Reset()
    bottomRefHist = bottomPadObjects[0]
    #bottomRefHist.Reset()
    
    
    
    
    dOpt = ""
    for s in topPadObjects:
        s.Draw(dOpt + getattr(s,"drawOption","") )
        dOpt = "same"


        
    topPad.Update()
    topPad.Draw()
    c1.Update()

    adjustHistAxis(topRefHist)
    topRefHist.GetYaxis().SetTitleOffset(y_title_offset)
    topRefHist.GetYaxis().SetMaxDigits(4)

    if ratio:
        bottomPad.cd()
       
        bottomRefHist.Draw()
        bottomPad.Update()
        #adjustRatioStyle(bottomPadObjects[0], plot, topPadObjects[0], scaleFacRatioPad)
        #adjustRatioStyle(1,2,3,4)
        adjustHistAxis(bottomRefHist)
        if plot and plot.xTitle:
            bottomRefHist.GetXaxis().SetTitle(plot.xTitle)
        
        if fix_overlapping_labels:
            topRefHist.GetYaxis().ChangeLabel(1,-1,0.) # remove the first axis label in Y 
            bottomRefHist.GetYaxis().ChangeLabel(-1,-1,0.) # remove the last axis label in Y
            bottomRefHist.GetYaxis().ChangeLabel(1,-1,0.) # remove the last axis label in Y
        
        #bottomRefHist.GetYaxis().SetNdivisions( int(np.ceil(bottomRefHist.GetMaximum())) ,0,2, False)
        #bottomRefHist.GetYaxis().SetNdivisions( 5 ,0,2, False)
        
        bottomRefHist.GetYaxis().SetLabelOffset(y_label_offset)
        bottomRefHist.GetXaxis().SetTitleOffset(x_title_offset)
        bottomRefHist.GetYaxis().SetTitleOffset(y_title_offset)
        #bottomRefHist.GetYaxis().SetNdivisions( int(np.ceil(bottomRefHist.GetMaximum())+1) ,2,0, False) 
        bottomRefHist.GetYaxis().SetNdivisions( 5,3,0) 
        
        dOpt = ""
        #bottomRefHist.Draw()
        for s in bottomPadObjects:
            s.Draw(dOpt + getattr(s,"drawOption","") )
            dOpt = "same"
            
        #bottomPadObjects[0].GetYaxis().SetTitleOffset(0.5)
        #bottomPadObjects[0].GetYaxis().SetTitleSize(0.08)
        
        
        bottomPad.Draw()
    if doDraw:
        c1.Draw()

    ret = {'canv':c1, 'pads': [topPad, bottomPad]}
    return ret


def convertUtoX(x, pad=None):
    pad = pad if pad else ROOT.gPad
    return (x-ROOT.gPad.GetX1())/(ROOT.gPad.GetX2()-ROOT.gPad.GetX1())
def convertVtoY(v, pad=None):
    pad = pad if pad else ROOT.gPad
    return (v-ROOT.gPad.GetY1())/(ROOT.gPad.GetY2()-ROOT.gPad.GetY1())

def drawBelle2Info(lumi=4.7, label="Phase 3 (Preliminary)", y=0.87, x=0.19, size=27, size_diff=2, y_diff=0.1, x_diff=0, font=43):
    tlatex = ROOT.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(font)  
    tlatex.SetTextSize(size)
    #tlatex.SetTextSize(0.060)
    #y = 0.87
    #x = 0.19
    #size_diff = 2
    if label:
        tlatex.DrawLatex(x, y, r"#bf{Belle II} %s"%label)
    tlatex.SetTextSize(size-size_diff)
    #tlatex.SetTextSize(0.050)
    if lumi:
        #tlatex.DrawLatex(x, y-y_diff, "#int#it{L}dt = %s pb^{-1} "%lumi)
        tlatex.DrawLatex(x-x_diff, y-y_diff, "#lower[0.00]{#int}#it{L}dt = %s fb^{-1} "%lumi)
    return tlatex
    


def getCanvPos(canvas, margin = 0.01):
    x1  = round( margin +  ( canvas.GetUxmin()-canvas.GetX1() ) / (canvas.GetX2()-canvas.GetX1()) , 3)
    x2  = round(-margin +  ( canvas.GetUxmax()-canvas.GetX1() ) / (canvas.GetX2()-canvas.GetX1()) , 3)
    y1  = round( margin +  ( canvas.GetUymin()-canvas.GetY1() ) / (canvas.GetY2()-canvas.GetY1()) , 3)
    y2  = round(-margin +  ( canvas.GetUymax()-canvas.GetY1() ) / (canvas.GetY2()-canvas.GetY1()) , 3)
    return (x1,y1,x2,y2)



def getDistributionFromList(l, nBins,  minVal=None, maxVal=None, name="dist", weighted=False):
    """
        Fill histogram with the values in list l
    """
    minVal = minVal if minVal else min(l)
    maxVal = maxVal if maxVal else max(l)
    h = ROOT.TH1D(name,name,nBins,minVal,maxVal)
    for v in l:
        v_ = getattr(v,"val",v)
        if weighted:
            sigma = v.sigma
            h.Fill( v_ , 1/sigma**2 )
        else:
            h.Fill( v_  )
    return h



def getFitError(h, cl=0.68, color=None, style=3544):
    """
    Has to be done right after the fit!!, i.e
    h.Fit("gaus")
    he = getFitError(h, cl=0.96)
    """
    he = h.Clone()
    if hasattr(he,'Reset'):
        he.Reset()
    ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(he, cl)
    he.SetFillColor( color if color else h.GetLineColor() )
    he.SetFillStyle( style )
    he.SetMarkerStyle(22)
    he.SetMarkerSize(0)
    return he



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
    os.makedirs(path, exist_ok=True)
    #if "." in path[-5:]:
    #    path = path.replace(os.path.basename(path),"")
    #    print (path)
    #if os.path.isdir(path):
    #    return
    #else:
    #    #mkdir_p(path)
    #    os.makedirs(path)



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

