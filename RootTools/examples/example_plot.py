''' Little demonstration how to plot.
Note: A stack is made from Samples (not from histos)
'''

#Standard imports
import ROOT
import os
from math import sqrt, cos

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

argParser.add_argument('--plotPath', 
      action='store',
      nargs='?',
      default='.',
      help="where to write the plots"
)

# RootTools
from RootTools.core.standard import *

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

#make samples. Samples are statisticall compatible.
s0 = Sample.fromFiles("", "example_data/file_0.root" )
# Add a cut to the sample instead of the plot. This is useful when requiring filters on data or triggers etc.
s1 = Sample.fromFiles("", "example_data/file_1.root" , selectionString = 'met_pt>20')
s2 = Sample.fromFiles("", "example_data/file_2.root" , selectionString = 'met_pt>40')

# styles are functions to be executed on the histogram
s0.style = styles.lineStyle( color = ROOT.kBlue )
s1.style = styles.lineStyle( color = ROOT.kRed )
s1.style = styles.fillStyle( color = ROOT.kGreen )

# scaling a sample
# Let's scale s2 up by a factor of 2 and compare it to s0+s1
s2.scale = 2 

# Define the stack.
# The whole point is that a stack is a stack of samples, not plots. Plots take stacks as arguments. 
stack = Stack( [ s0, s1], [ s2 ] )

# Let's use a trivial weight. All functions will 
plot_weight   = lambda event, sample : 1

# Two selection strings
selectionString = "nJet>0"
selectionString_2 = "nJet>1"

# Variables to be read from the tree
read_variables = []
#Sequence to be executed
sequence = []

plot1 = Plot(\
    name = "met_coarseBinning", # Name is not needed. If not provided, variable.name is used as filename instead.
    stack = stack,
    # met_pt is in the chain
    attribute = TreeVariable.fromString( "met_pt/F" ), 
    binning = [40,0,500], 
    selectionString = selectionString,
    weight = plot_weight 
)

def make_met_sqrtpt2( event, sample ):
    event.met_sqrtpt2 = sqrt(event.met_pt**2)

sequence.append( make_met_sqrtpt2 )
plot2 = Plot(\
    name = "met", # Name is not needed. If not provided, variable.name is used as filename instead.
    stack = stack,
    # met_pt is in the chain
    attribute = "met_sqrtpt2", 
    binning = [10,0,500], 
    selectionString = selectionString,
    weight = plot_weight 
)

read_variables +=  ["Jet[pt/F]", "met_pt/F"]
plot3 = Plot(\
    name = "met_plus_jet0Pt",
    stack = stack,
    # Here, I create a new variable that was not in the chain. 
    attribute = lambda event, sample: event.met_pt + event.Jet_pt[0], 
    binning = [20,0,500], 
    selectionString = selectionString,
    weight = plot_weight
)


plotting.fill([plot1, plot2, plot3], read_variables = read_variables, sequence = sequence)

if not os.path.exists( args.plotPath ): os.makedirs( args.plotPath )
for plot in [plot1, plot2, plot3]:
    plotting.draw(plot, plot_directory = args.plotPath, 
        ratio = {}, # Add a default ratio from the first two elements in the stack 
        logX = False, logY = True, 
# inside each stack member, sort wrt. histo.Integral. Useful when the stack looks like [ [mc1,mc2,mc3,...], [data]]
        sorting = True,
# Adjust y Range such that all bins are within canvas and that nothing is shadowed by the legend. Alternatively (yLow, yHigh)
        yRange = "auto",
# Make legend. Alternatively provide (x0,y0,x1,y1) 
        legend = "auto"  
        )
# Configuration of the ratio
# ratio = {'num':1, 'den':0, 'logY':False, 'style':None, 'texY': 'Data / MC', 'yRange': (0.5, 1.5)}
