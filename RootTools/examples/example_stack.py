''' Little demonstration how to define stacks.
Note: A stack is made from Samples (not from histos)
'''

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

# RootTools
from RootTools.core.standard import *

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#make a sample
s0 = Sample.fromFiles("", "example_data/file_0.root" )
s1 = Sample.fromFiles("", "example_data/file_1.root" , selectionString = 'met_pt>100')

stack = Stack( [s0, s1], [ s0 ] )
