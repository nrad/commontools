''' Sample class.
    Implements definition and handling of the TChain.
'''

# Standard imports
import ROOT
import uuid
import os
import random
from array import array
from math import sqrt
import subprocess

# Logging
import logging
logger = logging.getLogger(__name__)

# RootTools imports
import RootTools.core.helpers as helpers
import RootTools.plot.Plot as Plot
from   RootTools.core.SampleBase import SampleBase

# new_name method for sample counting
@helpers.static_vars( sample_counter = 0 )
def new_name():
    result = "Sample_"+str( new_name.sample_counter )
    new_name.sample_counter += 1
    return result

def check_equal_(vals):
    if not len(set(vals)) == 1:
        raise ValueError( "Sample combine check failed on: %r"%vals )
    else:
        return vals[0]

class Sample ( SampleBase ): # 'object' argument will disappear in Python 3

    def __init__(self, 
            name, 
            treeName , 
            files = [], 
            normalization = None, 
            xSection = -1,
            selectionString = None, 
            weightString = None,
            isData = False, 
            color = 0, 
            texName = None):
        ''' Handling of sample. Uses a TChain to handle root files with flat trees.
            'name': Name of the sample, 
            'treeName': name of the TTree in the input files
            'normalization': can be set in order to later calculate weights, 
            'xSection': cross section of the sample
            e.g. to total number of events befor all cuts or the sum of NLO gen weights
            'selectionString': sample specific string based selection (can be list of strings)
            'weightString': sample specific string based weight (can be list of strings)
            'isData': Whether the sample is real data or not (simulation)
            'color': ROOT color to be used in plot scripts
            'texName': ROOT TeX string to be used in legends etc.
        '''
        
        super(Sample, self).__init__( name=name, files=files, normalization=normalization, xSection=xSection, isData=isData, color=color, texName=texName)

        self.treeName = treeName
        self._chain = None
       
        self.__selectionStrings = [] 
        self.setSelectionString( selectionString )

        self.__weightStrings = [] 
        self.setWeightString( weightString )

        # Other samples. Add friend elements (friend, treeName)
        self.friends = []
             
        logger.debug("Created new sample %s with %i files, treeName %s,  selectionStrings %r and weightStrings %r.", 
            name, len(self.files), treeName, self.__selectionStrings, self.__weightStrings)

    def setSelectionString(self, selectionString):
        if type(selectionString)==type(""):
            self.__selectionStrings = [ selectionString ]
        elif type(selectionString)==type([]): 
            self.__selectionStrings =  selectionString 
        elif selectionString is None:
            self.__selectionStrings = []
        else:
            raise ValueError( "Don't know what to do with selectionString %r"%selectionString )
        logger.debug("Sample now has selectionString: %s", self.selectionString)
        self.clear()

    def addSelectionString(self, selectionString):
        if type(selectionString)==type(""):
            self.__selectionStrings += [ selectionString ] 
            self.clear()
        elif type(selectionString)==type([]): 
            self.__selectionStrings +=  selectionString 
            self.clear()
        #elif (selectionString is None ) or selectionString in [ [], "(1)", "" ] :
        elif (selectionString is None ) or selectionString == [] :
            pass 
        else:
            raise ValueError( "Don't know what to do with selectionString %r"%selectionString )

    def setWeightString(self, weightString):
        if type(weightString)==type(""):
            self.__weightStrings = [ weightString ]
        elif type(weightString)==type([]): 
            self.__weightStrings =  weightString 
        elif weightString is None:
            self.__weightStrings = []
        else:
            raise ValueError( "Don't know what to do with weightString %r"%weightString )
        logger.debug("Sample now has weightString: %s", self.weightString)
        self.clear()

    def addWeightString(self, weightString):
        if type(weightString)==type(""):
            self.__weightStrings += [ weightString ] 
            self.clear()
        elif type(weightString)==type([]): 
            self.__weightStrings +=  weightString 
            self.clear()
        elif (weightString is None ) or weightString in [ [], "(1)", "" ] :
        #elif (weightString is None ) or weightString == [] :
            pass 
        else:
            raise ValueError( "Don't know what to do with weightString %r"%weightString )

    @property
    def selectionString(self):
        return self.__selectionStrings if type(self.__selectionStrings)==type("") else helpers.combineStrings(self.__selectionStrings, stringOperator = "&&") 

    @property
    def weightString(self):
        return self.__weightStrings if type(self.__weightStrings)==type("") else helpers.combineStrings(self.__weightStrings, stringOperator = "*") 

    @classmethod
    def combine(cls, name, samples, texName = None, maxN = None, color = 0):
        '''Make new sample from a list of samples.
           Adds normalizations if neither is None
        '''
        if not (type(samples) in [type([]), type(())]) or len(samples)<1:
            raise ValueError( "Need non-empty list of samples. Got %r"% samples)

        normalizations = [s.normalization for s in samples]
        if None not in normalizations:
            normalization = sum(normalizations)
        else:
            normalization = None

        files = sum([s.files for s in samples], [])
        maxN = maxN if maxN is not None and maxN>0 else None
        files = files[:maxN]

        return cls(name = name, \
                   treeName = check_equal_([s.treeName for s in samples]),
                   xSection = check_equal_([s.xSection for s in samples]),
                   normalization = normalization,
                   files = files,
                   selectionString = check_equal_([s.selectionString for s in samples]),
                   isData = check_equal_([s.isData for s in samples]),
                   color = color, 
                   texName = texName
            )
 
    @classmethod
    def fromFiles(cls, name, files, 
        treeName = "Events", normalization = None, xSection = -1, 
        selectionString = None, weightString = None, 
        isData = False, color = 0, texName = None, maxN = None):
        '''Load sample from files or list of files. If the name is "", enumerate the sample
        '''

        # Work with files and list of files
        files = [files] if type(files)==type("") else files

        # If no name, enumerate them.
        if not name: name = new_name()

        # restrict files 
        maxN = maxN if maxN is not None and maxN>0 else None
        files = files[:maxN]

        sample =  cls(name = name, treeName = treeName, files = files, normalization = normalization, xSection = xSection,\
                selectionString = selectionString, weightString = weightString,
                isData = isData, color=color, texName = texName)

        logger.debug("Loaded sample %s from %i files.", name, len(files))
        return sample


    @classmethod
    def fromDirectory(cls, name, directory, treeName = "Events", normalization = None, xSection = -1, \
                selectionString = None, weightString = None,
                isData = False, color = 0, texName = None, maxN = None):
        '''Load sample from directory or list of directories. If the name is "", enumerate the sample
        '''
        # Work with directories and list of directories
        directories = [directory] if type(directory)==type("") else directory 

        # Automatically read from dpm if the directories indicate so
        if all( d.startswith('/dpm/') for d in directories ):
            return Sample.fromDPMDirectory( name=name, directory=directory, treeName=treeName, normalization=normalization, xSection=xSection,
                                            selectionString=selectionString, weightString=weightString, isData=isData, color=color, texName=texName, maxN=maxN) 

        # If no name, enumerate them.
        if not name: name = new_name()

        # find all files
        files = [] 
        for d in directories:
            fileNames = [ os.path.join(d, f) for f in os.listdir(d) if f.endswith('.root') ]
            if len(fileNames) == 0:
                raise helpers.EmptySampleError( "No root files found in directory %s." %d )
            files.extend( fileNames )
        if not treeName: 
            treeName = "Events"
            logger.debug("Argument 'treeName' not provided, using 'Events'.") 

        # restrict files 
        maxN = maxN if maxN is not None and maxN>0 else None
        files = files[:maxN]

        sample =  cls(name = name, treeName = treeName, files = files, normalization = normalization, xSection = xSection,\
            selectionString = selectionString, weightString = weightString,
            isData = isData, color=color, texName = texName)
        logger.debug("Loaded sample %s from %i files.", name, len(files))
        return sample
    

    def split(self, n, nSub = None, clear = True, shuffle = False):
        ''' Split sample into n sub-samples
        '''
        
        if n==1: return self

        if not n>=1:
            raise ValueError( "Cannot split into: '%r'" % n )

        files = self.files
        if shuffle: random.shuffle( files ) 
        chunks = helpers.partition( files, min(n , len(files) ) ) 

        if clear: self.clear() # Kill yourself.

        splitSamps = [Sample(
            name            = self.name + "_%i" % n_sample,
            treeName        = self.treeName,
            files           = chunks[n_sample],
            xSection        = self.xSection,
            normalization   = self.normalization,
            selectionString = self.selectionString,
            weightString    = self.weightString,
            isData          = self.isData,
            color           = self.color,
            texName         = self.texName) for n_sample in range(len(chunks))]

        if hasattr(self, 'json'):
            for s in splitSamps:
                s.json = self.json

        if nSub == None:
            return splitSamps 
        else:
            if nSub<len(chunks):
                return splitSamps[nSub]
            else:
                return None
        
    # Handle loading of chain -> load it when first used 
    @property
    def chain(self):
        if not self._chain: 
            logger.debug("First request of attribute 'chain' for sample %s. Calling __loadChain", self.name)
            self.__loadChain()
        return self._chain

    # Handle loading of rdf -> load it when first used
    @property
    def rdf(self):
        if not hasattr(self, "_rdf"):
            logger.debug("First request of attribute 'rdf' for sample %s. Getting the RDataFrame from chain", self.name)
            self._rdf = ROOT.RDataFrame(self.chain)
        return self._rdf

    # "Private" method that loads the chain from self.files
    def __loadChain(self):
        ''' Load the TChain. Private.
        '''
        if len(self.files) == 0:
            raise helpers.EmptySampleError("Sample {name} has no input files! Can not load.".format(name = self.name) )
        else:
            self._chain = ROOT.TChain(self.treeName)
            counter = 0
            for f in self.files:
                logger.debug("Now adding file %s to sample '%s'", f, self.name)
                try:
                    if helpers.checkRootFile(f, checkForObjects=[self.treeName]):
                        self._chain.Add(f)
                        counter+=1
                    else:
                        logger.error( "Check of root file failed. Skipping. File: %s", f )
                except IOError as e:
                    logger.error( "Could not load file %s", f )
                    #raise e
            if counter==0:
                raise helpers.EmptySampleError( "No root files for sample %s." %self.name ) 
            logger.debug( "Loaded %i files for sample '%s'.", counter, self.name )

        # Add friends
        if hasattr( self, 'friends'):  # Catch cases where cached samples have no default value for friends attribute
            for friend_sample, friend_treeName in self.friends:
                self.chain.AddFriend(friend_sample.chain, friend_treeName)

    # branch information
    @property
    def leaves( self ):
        ''' Get the leaves in the chain
        '''
        if hasattr( self, "__leaves" ):
            return self.__leaves
        else:
            self.__leaves = [ {'name':s.GetName(), 'type':s.GetTypeName()} for s in  self.chain.GetListOfLeaves() ]
            return self.__leaves 

    def clear(self): 
        ''' Really (in the ROOT namespace) delete the chain
        '''
        if self._chain:

            self._chain.IsA().Destructor( self._chain )
            logger.debug("Called TChain Destructor for sample '%s'.", self.name)

            self._chain = None

        if hasattr(self, "__leaves"):
            del self.__leaves

        return

    def sortFiles( self, sample, filename_modifier = None):
        ''' Remake chain from files sorted wrt. to another sample (e.g. for friend trees)
        '''

        # Check if file lists are identical
        filenames       = map(os.path.basename, self.files)
        other_filenames = [ f if filename_modifier is None else filename_modifier(f) for f in map(os.path.basename, sample.files) ]

        # Check if we have the same number of files
        if len(filenames)!=len(other_filenames):
            raise RuntimeError( "Can not sort files of sample %s according to sample %s because lengths are different: %i != %i", self.name, sample.name, len(self.files), len(sample.files) ) 

        new_filelist = []
        for f in other_filenames:
            # find position of file from other sample
            try:
                index = filenames.index(f)
            except ValueError:
                logger.error("Can not file %s from sample %s in files of sample %s", f, sample.name, self.name)
                raise

            new_filelist.append( self.files[index]  )

        # Destroy 
        self.clear()
        # Recreate files
        self.files = new_filelist
        return self

    def addFriend( self, other_sample, treeName, sortFiles = False):
        ''' Friend a chain from another sample.
        '''
        if sortFiles:
            other_sample.sortFiles( self )

        # Add Chains 
        self.friends.append( (other_sample, treeName) )

    def treeReader(self, *args, **kwargs):
        ''' Return a Reader class for the sample
        '''
        from RootTools.core.TreeReader import TreeReader
        logger.debug("Creating TreeReader object for sample '%s'.", self.name)
        return TreeReader( self, *args, **kwargs )

    # Below some helper functions to get useful 

    def combineWithSampleSelection(self, selectionString):
        if selectionString is None: return self.selectionString
        if not type(selectionString)==type(""): raise ValueError( "Need 'None' or string for selectionString, got %s" % selectionString )
        if self.__selectionStrings:
            logger.debug("For Sample %s: Combining selectionString %s with sample selectionString %s", \
                self.name, selectionString, self.selectionString )
            return helpers.combineStrings( [selectionString]+self.__selectionStrings, stringOperator = "&&")
        else:
            logger.debug("For Sample %s: Return selectionString %s because sample has no selectionString", \
                self.name, selectionString )
            return selectionString

    def combineWithSampleWeight(self, weightString):
        if weightString is None: return self.weightString
        if not type(weightString)==type(""): raise ValueError( "Need 'None' or string for weightString, got %s" % weightString )
        if self.__weightStrings:
            logger.debug("For Sample %s: Combining weightString %s with sample weightString %s", \
                self.name, weightString, self.weightString )
            return helpers.combineStrings( [weightString]+self.__weightStrings, stringOperator = "*")
        else:
            logger.debug("For Sample %s: Return weightString %s because sample has no weightString", \
                self.name, weightString )
            return weightString

    def getEventList(self, selectionString=None):
        ''' Get a TEventList from a selectionString (combined with self.selectionString, if exists).
        '''

        selectionString_ = self.combineWithSampleSelection( selectionString )

        tmp=str(uuid.uuid4())
        logger.debug( "Making event list for sample %s and selectionString %s", self.name, selectionString_ )
        self.chain.Draw('>>'+tmp, selectionString_ if selectionString_ else "(1)")
        elistTMP_t = ROOT.gDirectory.Get(tmp)

        return elistTMP_t


    def getChainHash(self):
        # get and order chain files
        # make sure all files are active? (nGetEntries?)
        # hash filenames+filesizes+nGetEntries
        pass

    def getSelectionStringHash(self, selectionString=None):
        selectionString_ = self.combineWithSampleSelection( selectionString )
        #hash(selectionString_)    

    def setEventList(self, selectionString=None):
        eListHash    = self.getSampleHash() + getSelectionStringHash(selectionString_)
        eListName = eListHash
        ## getEventList
        ##    if eList exists (sample attr or stored?):
        ##         getattr(chain, eListName)
        ##    else:
        ##          getEList
        ##
        pass


    def getYieldFromDraw(self, selectionString = None, weightString = None, split = 1, returnInfo = True, ignoreSampleSelectionString = False, ignoreSampleWeightString = False):
        ''' Get yield from self.chain according to a selectionString and a weightString
        ''' 

        if split > 1:
            results = [ subsample.getYieldFromDraw( selectionString = selectionString, weightString = weightString) for subsample in self.split( n = split, shuffle = True ) ]
            return {'val':sum( [r['val'] for r in results], ), 'sigma':sqrt( sum( [r['sigma']**2 for r in results], 0 ) ) }
        elif split == 1:
            selectionString_ = self.combineWithSampleSelection( selectionString )  if not ignoreSampleSelectionString  else selectionString
            weightString_    = self.combineWithSampleWeight( weightString )        if not ignoreSampleWeightString     else weightString

            tmp=str(uuid.uuid4())
            h = ROOT.TH1D(tmp, tmp, 1,0,2)
            h.Sumw2()
            #weight = weightString if weightString else "1"
            logger.debug( "getYieldFromDraw for sample %s with chain %r", self.name, self.chain )
            self.chain.Draw("1>>"+tmp, "("+weightString_+")*("+selectionString_+")", 'goff')
            res = h.GetBinContent(1)
            resErr = h.GetBinError(1)
            del h

            info = {'cut':selectionString_, 'weight':weightString_, 'name':self.name}
            ret = {'val': res, 'sigma':resErr}
            if returnInfo:
                ret.update(info)
            try:
                from PythonTools.u_float import u_float
                ret['u_float'] = u_float( ret['val'], ret['sigma'] ) 
            except:
                pass
            return ret
        else:
            raise ValueError( "Can't split into %r. Need positive integer." % split )

    def get1DHistoFromDraw(self, variableString, binning, selectionString = None, weightString = None, binningIsExplicit = False, addOverFlowBin = None, isProfile = False):
        ''' Get TH1D/TProfile1D from draw command using selectionString, weight. If binningIsExplicit is true, 
            the binning argument (a list) is translated into variable bin widths. 
            addOverFlowBin can be 'upper', 'lower', 'both' and will add 
            the corresponding overflow bin to the last bin of a 1D histogram.
            isProfile can be True (default) or the TProfile build option (e.g. a string 's' ), see
            https://root.cern.ch/doc/master/classTProfile.html#a1ff9340284c73ce8762ab6e7dc0e6725'''

        selectionString_ = self.combineWithSampleSelection( selectionString )
        weightString_    = self.combineWithSampleWeight( weightString )

        tmp=str(uuid.uuid4())
        if binningIsExplicit:
            binningArgs = (len(binning)-1, array('d', binning))
        else:
            binningArgs = binning

        if isProfile:
            if type(isProfile) == type(""):
                res = ROOT.TProfile(tmp, tmp, *( binningArgs + (isProfile,)) )
            else:
                res = ROOT.TProfile(tmp, tmp, *binningArgs)
        else:
                res = ROOT.TH1D(tmp, tmp, *binningArgs)

        #weight = weightString if weightString else "1"

        self.chain.Draw(variableString+">>"+tmp, "("+weightString_+")*("+selectionString_+")", 'goff')
       
        Plot.addOverFlowBin1D( res, addOverFlowBin )
 
        return res

    def get2DHistoFromDraw(self, variableString, binning, selectionString = None, weightString = None, binningIsExplicit = False, isProfile = False):
        ''' Get TH2D/TProfile2D from draw command using selectionString, weight. If binningIsExplicit is true, 
            the binning argument (a tuple of two lists) is translated into variable bin widths. 
            isProfile can be True (default) or the TProfile build option (e.g. a string 's' ), see
            https://root.cern.ch/doc/master/classTProfile.html#a1ff9340284c73ce8762ab6e7dc0e6725
        '''

        selectionString_ = self.combineWithSampleSelection( selectionString )
        weightString_    = self.combineWithSampleWeight( weightString )

        tmp=str(uuid.uuid4())
        if binningIsExplicit:
            if not len(binning)==2 and type(binning)==type(()):
                raise ValueError( "Need a tuple with two lists corresponding to variable bin thresholds for x and y axis. Got % s"% binning )
            binningArgs = (len(binning[0])-1, array('d', binning[0]), len(binning[1])-1, array('d', binning[1]))
        else:
            if not len(binning)==6:
                raise ValueError( "Need binning in standard 2D form: [nBinsx,xLow,xHigh,nBinsy,yLow,yHigh]. Got %s" % binning )
            binningArgs = binning

        if isProfile:
            if type(isProfile) == type(""):
                res = ROOT.TProfile2D(tmp, tmp, *( binningArgs + (isProfile,)) )
            else:
                res = ROOT.TProfile2D(tmp, tmp, *binningArgs)
        else:
                res = ROOT.TH2D(tmp, tmp, *binningArgs)

        self.chain.Draw(variableString+">>"+tmp, "("+weightString_+")*("+selectionString_+")", 'goff')

        return res

    def get3DHistoFromDraw(self, variableString, binning, selectionString = None, weightString = None, binningIsExplicit = False, isProfile = False):
        ''' Get TH3D/TProfile3D from draw command using selectionString, weight. If binningIsExplicit is true, 
            the binning argument (a tuple of two lists) is translated into variable bin widths. 
            isProfile can be True (default) or the TProfile build option (e.g. a string 's' ), see
            https://root.cern.ch/doc/master/classTProfile.html#a1ff9340284c73ce8762ab6e7dc0e6725
        '''

        selectionString_ = self.combineWithSampleSelection( selectionString )
        weightString_    = self.combineWithSampleWeight( weightString )

        tmp=str(uuid.uuid4())
        if binningIsExplicit:
            if not len(binning)==3 and type(binning)==type(()):
                raise ValueError( "Need a tuple with three lists corresponding to variable bin thresholds for x, y and z axis. Got % s"% binning )
            binningArgs = (len(binning[0])-1, array('d', binning[0]), len(binning[1])-1, array('d', binning[1]),  len(binning[2])-1, array('d', binning[2]))
        else:
            if not len(binning)==9:
                raise ValueError( "Need binning in standard 3D form: [nBinsx,xLow,xHigh,nBinsy,yLow,yHigh,nBinsz,zLow,zHigh]. Got %s" % binning )
            binningArgs = binning

        if isProfile:
            logger.warning( "Not sure TTree::Draw into TProfile3D is implemented in ROOT." )
            if type(isProfile) == type(""):
                res = ROOT.TProfile3D(tmp, tmp, *( binningArgs + (isProfile,)) )
            else:
                res = ROOT.TProfile3D(tmp, tmp, *binningArgs)
        else:
                res = ROOT.TH3D(tmp, tmp, *binningArgs)

        self.chain.Draw(variableString+">>"+tmp, "("+weightString_+")*("+selectionString_+")", 'goff')

        return res
