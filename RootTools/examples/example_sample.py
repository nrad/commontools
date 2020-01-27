''' A simple example that demonstrates how to create Sample instances from data files.
'''
# Standard imports
import sys
import os

#import logging
import ROOT

# RootTools
from RootTools.core.standard import *

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, None)

filename =  "example_data/file_0.root"

# simplest way
s0 = Sample.fromFiles("", filename )
s0.chain

# works as well with lists:
s0_2 = Sample.fromFiles("", [filename, filename], selectionString = ["met_pt>200", "Jet_pt[0]>100"])
s0_2.chain

afsRootToolsExamples =  os.path.expandvars( "$CMSSW_BASE/src/RootTools/examples/example_data" ) 

# from files
s1 = Sample.fromFiles("TTZToQQ", files = [os.path.join( afsRootToolsExamples, "file_0.root" )], treeName = "Events")
s1.chain

## from CMG output
#td = os.path.join( afsRootToolsExamples, "DoubleMuon_Run2016D-PromptReco-v2" )
#s2 = Sample.fromCMGOutput("DoubleMuon_Run2016D", baseDirectory = td, treeName = "tree")
#s2.chain

## from a directory with root files
#s3=Sample.fromDirectory(name="same_as_s2", directory = os.path.join( afsRootToolsExamples, "sample_data" ), treeName = "Events")
#s3.chain
