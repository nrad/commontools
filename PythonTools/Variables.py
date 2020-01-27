
from PythonTools.RecursiveFormatter import RecursiveFormatter
from copy import deepcopy
import re 


def selectKeysFromDict( di, keys=[] ):
    return {k:di[k] for k in keys}

def varsToFormat( varstring ):
    return re.findall( r'{(.*?)}', varstring ) 

def replaceItemInList( lst, this, withthat):
    if not this in lst:
        return lst
    idx = lst.index(this) 
    return lst[0:idx] + ( withthat if isinstance(withthat, (tuple,list) ) else [withthat] )  + lst[idx+1:]

def expandItemsInList( lst, mapdict):
    """
        replaces items in list if if they match to a key in mapdict (semi-iteratively)
        newlist = expandItemsInList( ['var1','var2'] , {'var1':['var4','var5'], 'var4':'var6' } )
    """
    newlist = []
    i=0
    while (lst != newlist):
        i+=1
        if i==20:
            raise Exception("Got stuck in a loop somehow!")
        lst = newlist if newlist else lst
        #print 'cutlist', cutlist
        #print 'cutlist_', newlist
        for c in lst:
            #print '   ', c
            for k,v in mapdict.items():
                if "{%s}"%k==c:
                    #print k
                    newlist = replaceItemInList(lst, "{%s}"%k , getattr(v,'list',v)   )
        if not newlist:
            return lst
    return newlist



def formatVar(  var, vars_dict ):
    vars_dict_ = {k: v if not isinstance(v,(tuple,list)) else combineCutStrings(v) for k,v in vars_dict.items() }
    var_ = var if not isinstance(var, (tuple,list)) else combineCutStrings(var)
    formatted = RecursiveFormatter().format( var_, **vars_dict_ )
    return formatted


def isVarInstance(var):
    isVar = isinstance(var, Variable) or getattr(var, 'isVariable', False)
    return isVar

def findRelevantVarsToFormat( var , vars_dict ):
    #isVarInst = isinstance(var, Variable)
    isVarInst = isVarInstance(var)
    varname   = var.name   if isVarInst else var
    varstring = var.string if isVarInst else var
    relevant_vars = varsToFormat( varstring )
    if not relevant_vars:
        return relevant_vars
    for v in relevant_vars:
        if v not in vars_dict:
            raise Exception("Variable %s (%s) depends on definition of (%s) but it was not found in the vars_dict given!"%(varname,varstring, v))
        relevant_vars.extend( [v_ for v_ in findRelevantVarsToFormat( vars_dict[v], vars_dict ) if v_ not in relevant_vars] )
    return relevant_vars 


def combineCutStrings(cutStringList):
  #print cutStringList
  return "%s"%( " && ".join([ "%s"%c for c in cutStringList if c not in ["(1)", "((1))"] ]) )


def getFlowFromList(cutList):
    flow = []
    for i,l in enumerate(cutList):
        flow.append(cutList[0:i+1] )
    return flow

class Variable(object):
    """
        For use within the Variables class only
    """
    def __init__(self, name, string, latex=None, relevant_vars_dict={} ):
        self.name   = name
        self.isVariable = True
        if isinstance(string, (tuple,list)):
            self.list = string
            self.string = string #combineCutStrings(string)
            self.isRegion = True
        else:
            self.list   = [string]
            self.string = str(string)
            self.isRegion = False
        self.latex  = latex if latex else ''
        self.relevant_vars_dict = relevant_vars_dict 

    def format(self, **kwargs):
        newargs = deepcopy( self.relevant_vars_dict )
        newargs.update( kwargs )
        return formatVar( self.string, newargs )
    @property
    def formatted(self):
        return formatVar( self.string, self.relevant_vars_dict )
        #return formatVar( self.expandedVar, self.relevant_vars_dict )
    def __str__(self):
        #return self.string
        return self.formatted
    def __repr__(self):
        return "<%s.%s %s: %s >"%(self.__module__, self.__class__.__name__, self.name, self.string)

    def replaceCut(self, this, withthat, name=None, relevant_vars_dict={}):
        newlist = replaceItemInList( self.list, this, withthat )
        return Variable( name   = name if name else self.name, 
                         string = newlist,
                         relevant_vars_dict = relevant_vars_dict if relevant_vars_dict else self.relevant_vars_dict ,
                          )

    def expandList(self, relevant_regions_dict={}):
        relevant_regions_dict = relevant_regions_dict if relevant_regions_dict else self.relevant_regions_dict
        cutlist = self.list[:]
        newlist = [] 
        newlist = expandItemsInList( cutlist, relevant_regions_dict )
        return newlist

    @property
    def relevant_regions(self):
        return [ k for k,v in self.relevant_vars_dict.items() if (getattr(v, 'isRegion', False) or isinstance(v, (tuple,list))) ]

    @property
    def relevant_regions_dict(self):
        return selectKeysFromDict( self.relevant_vars_dict, self.relevant_regions ) 

    @property
    def expandedList(self):
        return self.expandList( self.relevant_regions_dict )

    @property
    def expandedVar(self):
        return Variable( self.name, self.expandedList, relevant_vars_dict=self.relevant_vars_dict)

    def expandVar(self, relevant_vars_dict = {} ):
        raise Exception("Not implemented yet.. issues with loop getting stuck in expandItemsInList")
        return False

    
    def clone(self):
        return Variable( self.name, newlist, relevant_vars_dict=self.relevant_vars_dict)


class Variables():
    """
    vars_dict = { 
                    'var1': 'LepIndex_{someopt}' ,
                    'var2': 'GenPart_pt{var1}',            # treated as "cut"
                   'sr'   : ['{var1}>200', '{var1}<300'],  # treated as a "region"
                }
    
    variables = Variables( vars_dict, someopt=123)
    or
    variables = Variables( vars_dict, **{'someopt':123})
    
    
        variables.var1.formatted
        variables.var1.format( someopt='SOMETHINGELSE_{newopt}' , newopt='var3' )
    
    definitions can be changed globally later on by:
        variables._update( lep='mu')
        or 
        variables._update( **{'lep'='mu','lt':'def'})

    for regions it might be usefull to do:
  
    variables.sr.list  
    
    """

    def __init__( self, vars_dict = None , verbose = False, **kwargs):
        self.verbose = verbose
        self.vars_dict = {}
        self._update(vars_dict, **kwargs)
        self.vars_dict_orig = deepcopy( vars_dict )

    def _update(self, vars_dict = {} , **kwargs):
        self.vars_dict.update( self._make_vars_dict( vars_dict )  )
        kwarg_vars_dict = self._make_vars_dict( kwargs )
        self.vars_dict.update( kwarg_vars_dict )

        for k,v in self.vars_dict.items():
            var = Variable( name=k, string=v , relevant_vars_dict = self.vars_dict )
            if var.isRegion:
                self.vars_dict[k] = var.string
            #vstring = vstring if not isinstance(vstring, (tuple,list)) else combineCutStrings(vstring)
            setattr(self, k, var )

    def _make_vars_dict(self, vars_dict):
        di = {}
        for k,v in vars_dict.items():
            vargs = {'name':k}
            if isVarInstance(v):
                vargs['string'] = v.string
            elif not isinstance(v,dict):
                vargs['string'] = v
            elif not v.get('string'): 
                if v.has_key('var'):
                    v['string'] = v.pop('var')
            else:
                vargs.update(v)
            #di[k] = Variable(**vargs)
            vstring = vargs['string']
            #vstring = vstring if not isinstance(vstring, (tuple,list)) else combineCutStrings(vstring)
            di[k] = vstring
        #di = deepcopy( vars_dict )
        return di

    def _addFormatFunc(self, var):
        def func():
            return formatVar( var, self.vars_dict )
        func.__name__ = "format_func_%s"%var.name
        setattr( var, "format", func)


if __name__ == '__main__':
    #print (__doc__)
    vars_dict_ = {
                    
    
                     'lepIndex': 'Index{lepCol}_{lep}{lt}',
                     'lepIndex1': '{lepIndex}[0]',
                     'lepIndex2': 'Max$(Alt$(Index{lepCol}_{lep}{lt}[1],-999))',
                     'lepPt': '{lepCol}_pt[{lepIndex1}]',
                     'lepPhi': '{lepCol}_phi[{lepIndex1}]',
                     'lepIndex_tight': 'Index{lepCol}_{lep}{tightWP}',
                     'lepIndex_tight1': '{lepIndex_tight}[0]',
                     'lepIndex_tight_lep': 'Index{lepCol}_lep{tightWP}',
                     'lepIndex_veto': 'Index{lepCol}_lep{lt}',
                     'lepMT': '{lepCol}_mt[{lepIndex1}]',
                     'lepIndex_lep': 'Index{lepCol}_lep{lt}',
                     'lepIndex_lep1': '{lepIndex_lep}[0]',

    }
    
    kwargs_ = {
                     'lepCol':'LepGood',
                     'lep'   :'lep',
                     'lt'    :'_lowpt',
                     'tightWP':'TIGHT',
                    'looseWP':'LOOSE',
             }
    
    vs= Variables( vars_dict_ , **kwargs_ )
    
