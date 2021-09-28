''' A stack of samples (not plots).
Must be a list of lists.
'''
# Standard imports
import uuid
import array

# RootTools
from RootTools.core.Sample import Sample
from RootTools.plot.Plot import Plot
from RootTools.plot.Immutable import Immutable
from RootTools.plot.Binning import Binning

class Stack ( list ):
        
    def __init__(self, *stackList):

        # change [[...], X, [...] ...]  to [[...], [X], [...], ...]
        stackList = [ s if type(s)==type([]) else [s] for s in stackList]

        # Check the input. LBYL.
        for s in stackList:
            if not type(s)==type([]) or not all(isinstance(p, Sample) for p in s):
                raise ValueError("Stack should be a list of lists of Sample instances. Got %r."%( stackList ) )

        # Make Immutable
        # stackList = map(lambda l: map(lambda s:Immutable(s), l), stackList)

        super(Stack, self).__init__( stackList )

    @property
    def samples(self):
        ''' Get all unique samples for this stack
        '''
        return list(set(sum(self,[])))

    def make_histos(self, plot):
        '''Make histograms for plot for this stack. Structure is list of lists of histos parallel to the stack object
        '''
        res = []
        for i, l in enumerate(self):
            histos = [] 
            for j, s in enumerate(l):
                if isinstance(plot.binning, Binning):
                    if plot.binning.binning_is_explicit:
                        # explicit binning with thresholds
                        histo = plot.histo_class(\
                            "_".join([plot.name, s.name, str(uuid.uuid4()).replace('-','_')]), 
                            "_".join([plot.name, s.name]), 
                             len(plot.binning.binning)-1, array.array('d', plot.binning.binning) )
                    else:
                        # default case (but using Binning class)
                        histo = plot.histo_class(\
                            "_".join([plot.name, s.name, str(uuid.uuid4()).replace('-','_')]), 
                            "_".join([plot.name, s.name]), 
                             *plot.binning.binning )
                elif type(plot.binning)==type([]) or type(plot.binning)==type(()): 
                    # default case: Binning is specified as [n, x0, x1]
                    histo = plot.histo_class(\
                        "_".join([plot.name, s.name, str(uuid.uuid4()).replace('-','_')]), 
                        "_".join([plot.name, s.name]), 
                         *plot.binning )
                else:
                    raise ValueError( "Don't know what to do with binning of plot %s: %r"%( plot.name, plot.binning ) )

                histo.Reset()

                # Default sumW2
                try: histo.Sumw2()
                except: pass

                # Exectute style function on histo
                if hasattr(s, "style"):
                    s.style(histo)

                histos.append(histo)
            res.append(histos)
            
        return res 

    def getSampleIndicesInStack(self, sample):
        ''' Find the indices of a sample in the stack
        '''
        indices = []
        for i, l in enumerate(self):
            for j, s in enumerate(l):
                if s==sample: indices.append((i,j))
        return indices

    def applyFuncInParal(self, func=None, nProc=16):
        from PythonTools.NavidTools import runFuncInParal
        flat =  Stack.flatten(self)
        ret_flat  =  runFuncInParal(func, flat, nProc=nProc, verbose=False)
        samps_idx = list( map( self.getSampleIndicesInStack, self.flatten(self) )) 

        ret = []
        for iret, [(i,j)] in enumerate(samps_idx):
            if len(ret) <= i:
                ret.append([])
            ret[i].append( ret_flat[iret] )
 
        return ret

    def applyFunc(self, func=None):
        """
            Applies func to all samples in stack and returns the outputs in the same structure as the stack
        """
        return self.applyFuncToStack(self, func = func)

    @staticmethod
    def applyFuncToStack(stackTypeList, func = None):
        """
            stack method:
            stackTypeList must be list of lists.
            Applies func to each element in stackTypeList and returns the outputs in the same structure as the stackTypeList
            
        """
        ret = []
        for i,l in enumerate(stackTypeList):
            ret.append([])
            for j,s in enumerate(l):
                ret[-1].append( func(s) )
        return ret

    @staticmethod
    def applyFuncToStackByIndex(stackTypeList, func = None):
        """
            stack method:
            stackTypeList must be list of lists.
            Applies func to each element in stackTypeList and returns the outputs in the same structure as the stackTypeList
            func should be of the form:
            func = lambda i,j, s: ....
            where stackTypeList[i][j]=s
 
        """
        ret = []
        for i,l in enumerate(stackTypeList):
            ret.append([])
            for j,s in enumerate(l):
                ret[-1].append( func(i,j,s) )
        return ret

    #def getHistos(self, var, binning, cut=None, weight=None, name=None, title=None):
    #    from PythonTools.ROOTHelpers import getHistoFromRDF
    #    #                         
    #    getHistoFromRDF(
    #    return self.applyFunc( lambda x: getHistoFromRDF(x.rdf, *args, **kwargs) )


    @staticmethod
    def flatten(stackTypeList):
        import itertools
        return list(itertools.chain(*stackTypeList))

    def getYieldsFromDraw(self, selectionString = None, weightString = None, returnInfo = True, returnDataFrame=False):
        ret = self.applyFunc( lambda s: s.getYieldFromDraw( selectionString, weightString, returnInfo=returnInfo if not returnDataFrame else True) )
        if returnDataFrame:
            import pandas
            df = pandas.DataFrame( self.flatten(ret) )
            df.set_index('name', inplace=True)
            ret = df
        return ret



