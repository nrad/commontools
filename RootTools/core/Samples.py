"""

    Samples Class for easir handling of multiple sample files.
    nice features:
        samples.stack.getHistos( .... )
    somewhat compatible with RDataFrame
    
    TO DO:
        remove the outside RootTools dependances:
            param module
            PythonTools.ROOTHelper dependencies

"""



#import core.standard as RootTools
from RootTools.core.Sample import check_equal_
from RootTools.core.Sample import Sample 
from RootTools.plot.Stack import Stack
from RootTools.plot.Color import Color
from RootTools.plot import styles


import itertools
import param
#import NavidTools.NavidTools as nt
import glob
from copy import deepcopy
import os
from collections import OrderedDict




class SampleParam(param.Parameterized):
    reco_tag = param.String()
    mc_tag = param.String()
    proc_tag = param.String()
    postfix  = param.String(default="")
    #exp    = param.Integer()
    isData = param.Boolean(default=False)
    isBkg = param.Boolean(default=False)
    isSignal = param.Boolean(default=False)
    color    = param.Integer(default=1)
    style    = param.Integer(default=0)
    name   = param.String()
    sample_name = param.String()
    data_tag    = param.String()
    xsec   = param.Number()
    lumi   = param.Number()
    files_pattern = param.String()
    base_dir = param.String()
    n_files = param.Integer()
    tree_name = param.String(default="tau3x1")
    weight_string = param.String(default="(1)")
    selection_string = param.String(default="(1)")
    lumi = param.Number(default=1.0)
    texName = param.String() 
    
    def _get_files(self):
        import glob
        #self._files_pattern = self.files_pattern.format( mc_tag=self.mc_tag, name=self.name)
        
        self._files_pattern = os.path.join( self.base_dir, self.files_pattern).format( ** self._dict() )
        self._files = glob.glob( self._files_pattern)
        if not self._files:
            print("WARNING: no files found for sample %s in %s"%(self.name, self._files_pattern))
                    
        return self._files
        #self._files = glob.glob( 

    def _get_df(self):
        from root_pandas import read_root
        self.df = read_root( self._files, key=self.tree_name) 


    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.sample_name = self.sample_name if self.sample_name else self.name
        
        self._get_files()

    def _dict(self):
        return dict( self.get_param_values() )



def makeSampParams(sample_set , **kwargs):
    samp_params = {}
    for samp, samp_info in sample_set.items():
        samp_info = deepcopy(samp_info)
        names = samp_info.pop('names')
        samp_info.update(kwargs)
        #print(samp_info)
        weight_dict = samp_info.pop("weight_dict") if "weight_dict" in samp_info else {}
        lumi_dict   = samp_info.pop("lumi_dict") if "lumi_dict" in samp_info else {}
        for name in names:
            samp_info_ = deepcopy(samp_info)
            if weight_dict:
                weight_string = weight_dict[name]
                samp_info_.update({'weight_string':weight_string})
            if lumi_dict and name in lumi_dict:
                lumi = lumi_dict[name]
                samp_info_.update({'lumi':lumi})
            samp_params[name] = SampleParam( name=name, **samp_info_ )
    return samp_params      

def makeRootToolsSampFromParams(samp_param ):

    #from RootTools.core.Sample import Sample
    samp = Sample.fromFiles(samp_param.name, 
                               samp_param._files,  
                               weightString=samp_param.weight_string,
                               selectionString=samp_param.selection_string,
                               treeName=samp_param.tree_name,
                               texName=samp_param.texName,
                                )
    return samp

def makeSampsAll( samp_params , samp_styles={}, strict=True):
    samps_all = {}
    for sname, samp_param in samp_params.items():
        style = samp_styles[sname] if sname in samp_styles else None
        if not samp_param._files and not strict:
            print(f"O.K.... no files found for {sname}...and apparently you didn't wanna be too strict, so will skip this.") 
            continue
        samps_all[sname] = makeRootToolsSampFromParams( samp_param )
        for tag in ['isData', 'isBkg','isSignal']:
            setattr( samps_all[sname], tag, getattr(samp_param,tag, False) )
    return samps_all


def getDataBkgSigTag(samp):
    tags = ['isData', 'isBkg','isSignal']
    ret  = {tag :  getattr(samp, tag, False) for tag in tags }
    if (sum(ret.values()) <= 1):
        ret['isBkg']=True
    return ret

 

def makeCombinedSamps( samps_all, samps_to_combine, samp_styles):
    samps = {}
    samps_to_combine = deepcopy(samps_to_combine)
    for sname, slist in samps_to_combine.items():
        tags  = {}
        for s in slist:
            #print (sname,s, tags)                
            if s not in samps_all:
                print("WARNING: Failed combining sample %s in %s \n %s"%(s,sname, samps_all.keys() ))
                slist.pop(slist.index(s))
                continue
            sample_tags = getDataBkgSigTag( samps_all[s] )
            if tags:
                if not tags == sample_tags:
                    raise Exception("Inconsistant samples are being combined together for %s: %s vs %s"%(sname, tags,sample_tags))
            else:
                tags.update( sample_tags )
        if not slist:
            print("WARNING: Didn't find matching sample for %s"%sname)
            continue
        samps[sname] = Sample.combine( sname, [samps_all[s] for s in slist] )
       
        samps[sname].addWeightString(check_equal_([samps_all[s].weightString for s in slist]))
        for tag, val in tags.items():
            #print (sname,tag,val)
            setattr(samps[sname], tag, val)

        if sname in samp_styles:
            samps[sname].style = samp_styles[sname]
        else:
            print(" sample %s has no style! :("%(sname) )
    return samps








class Samples():
    #import RootTools.core.standard as RootTools
    def __init__(self, name , sample_set, settings=None, samps_to_combine=[], samp_styles=[], tree_name='tau3x1', strict=True):

        self.name = name
        self._sample_set = sample_set
        self._settings   = settings
        self._samps_to_combine = samps_to_combine
        self._samp_styles      = samp_styles
        


        self._samp_params = makeSampParams( self._sample_set, **self._settings)
        #self._samps_all   = makeSampsAll( self._samp_params, tree_name)
        self._samps_all   = makeSampsAll( self._samp_params, self._samp_styles, strict=strict)

        self._sig_list_all  = [sname for sname,s  in self._samps_all.items() if getattr(s, 'isSignal', False) ]
        self._bkg_list_all  = [sname for sname,s  in self._samps_all.items() if getattr(s, 'isBkg', False) ]
        self._data_list_all = [sname for sname,s  in self._samps_all.items() if getattr(s, 'isData', False) ]

        self._samps       = makeCombinedSamps( self._samps_all, self._samps_to_combine, self._samp_styles) if samps_to_combine else self._samps_all

        self._sig_list   = [sname for sname,s  in self._samps.items() if getattr(s, 'isSignal', False) ]
        self._bkg_list   = [sname for sname,s  in self._samps.items() if getattr(s, 'isBkg', False) ]
        self._data_list  = [sname for sname,s  in self._samps.items() if getattr(s, 'isData', False) ]

        
        for sname, s in self._samps.items():
            setattr(self, sname,s)


    def _addSample(self, **kwargs):
        if isinstance(kwargs, dict):
            sample_info = SampleParam(**kwargs)
        elif isinstance(kwargs, SampleParam):
            sample_info = kwargs
        else:
            raise Exception(f"arguments are not compatible with a SampleParam instance ({kwarg})")

        sample = makeRootToolsSampFromParams( sample_info )
        sname = sample_info.name
        
        #color = getattr(sample_info,"color") if getattr(sample_info,"color") else Color()
        color = Color(getattr(sample_info,"color")) #if getattr(sample_info,"color") else Color()
        style={'markerStyle':0}
        if sample_info.isBkg:
            self._bkg_list_all.append(sname)
            self._bkg_list.append(sname)
            style.update(**{'lineColor':color, 'lineWidth':2 })
        elif sample_info.isSignal:
            self._sig_list_all.append(sname)
            self._sig_list.append(sname)
            style.update(**{'lineColor':color, 'lineStyle':3, 'lineWidth':2 })
        elif sample_info.isData:
            self._data_list_all.append(sname)
            self._data_list.append(sname)
            style.update(**{'lineColor':1, 'lineWidth':2, 'markerStyle':20, 'markerSize':1 })
        else:
            raise Exception("Couldnt find any of the flags isBkg, isSignal or isData.... but maybe I should continue?")
        lineStyle = getattr(sample_info, 'style')
        if lineStyle:
            style.update(**{'lineStyle':lineStyle})
        print(kwargs)
        print("STYLE", style)
        sample.style = styles.styler(**style)
        self._samps_all[sname] = sample
        self._samps[sname] = sample
        setattr(self, sname, sample)

    def _getStack(self, sig_list=None, bkg_list=None, data_list=None ):
        sig_list   = self._sig_list  if sig_list  == None else sig_list
        bkg_list   = self._bkg_list  if bkg_list  == None else bkg_list
        data_list  = self._data_list if data_list == None else data_list
        stack_list = [ [ self._samps[s] for s in sig_list] + [ self._samps[s] for s in bkg_list] ]+\
                     [ [ self._samps[s] for s in data_list]]
        #print( stack_list )
        stack = Stack( * [ s for s in stack_list if s ] )
        #self._stack = stack
        #setattr(self, '_stack',  stack )
        return stack


    @property
    def stack(self):
        if not hasattr(self,'_stack'):
            #self._getStack()
            setattr(self, '_stack',  self._getStack() )
        return self._stack

    def getStackHistosRDF(self, *args,**kwargs):
        from PythonTools.ROOTHelpers import getRDFHistoFromSample
        return self.stack.applyFunc( lambda x: getRDFHistoFromSample(x, *args, **kwargs) )

    def getStackHistos(self, *args,**kwargs):
        from PythonTools.ROOTHelpers import getHistoFromSample
        #def getPlotFromChain(c, var, binning, cutString = "(1)", weight = "weight", binningIsExplicit=False ,addOverFlowBin='',variableBinning=(False, 1) , name=None, nEvents=None)
        return self.stack.applyFunc( lambda x: getHistoFromSample(x, *args, **kwargs) )


    def getElists(self, cut, resetEList=True, retrieve=False):
        #def getAndSetEList(chain, cut, eListName=None, resetEList=True, retrieve=False):
        from PythonTools.ROOTHelpers import getAndSetEList
        self.stack.applyFunc(lambda s: getAndSetEList(s.chain, cut, resetEList=resetEList, retrieve=retrieve) )

    def getHistos(self, plots, cut=None, weight=None, rdf=True, nEvents=None):
        rhistos = {}
        #Register the plots in the RDF (Fast step)
        histos = {}
        if nEvents:
            print(f"-------- WARNINGL: Will only plot {nEvents} of the events..THE NORMALIZATINS WILL BE WRONG! This option should only be used for tests!")
            if rdf: 
               raise Exception("nEvents option does not work with RDFs")

        if rdf:
            for p in plots:
                rhistos[p.name] = self.getStackHistosRDF(p.var, tuple(p.bins), cut=str(cut) if cut else cut, weight=weight )
            print("Registered Histograms in the RDFs (Lazy action)")
            #actually get the plots in the next step... should only take "time" for the first plot per sample
            print(f"Now getting them")
            for p in plots:
                histos[p.name]  = self.stack.applyFuncToStack(rhistos[p.name], lambda x: x.GetValue().Clone() )
                self.stack.applyFuncToStackByIndex(histos[p.name], lambda i,j,s: [ 
                                                                       setattr(histos[p.name][i][j], "cut",    rhistos[p.name][i][j].cut   ),
                                                                       setattr(histos[p.name][i][j], "weight", rhistos[p.name][i][j].weight) ,
                                                                       #setattr(histos[p.name][i][j], "ptr", rhistos[p.name][i][j]) ,
                                                                              ])
        else:
            for p in plots:
                print(f"Getting Histogram: {p.name}")
                histos[p.name] = self.getStackHistos(p.var, p.bins, cutString=str(cut) if cut else cut, weight=str(weight) if weight else weight, nEvents=nEvents )
                #self.applyFuncToStacByIndex(histos[p.name], lambda i,j,h: histos[p.name][i][j].SetTitle( self.stack[i][j].texName) )
        for p in plots:
            print(f"Plotting: {p.name}")
            #histos[p.name]  = self.stack.applyFuncToStack(rhistos[p.name], lambda x: x.GetValue().Clone() )
            self.stack.applyFuncToStackByIndex(histos[p.name], lambda i,j,s: self.stack[i][j].style(s) if hasattr(self.stack[i][j], 'style') else None )
            self.stack.applyFuncToStackByIndex(histos[p.name], lambda i,j,s: s.SetTitle(self.stack[i][j].texName) )
        return histos




