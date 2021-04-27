"""

use functions like this which take chain, and event as arguments
and set the SOMEVARIABLE as the attribute of event:

def someFunc(chain, event, *args, **kwarg):
    <so some physics>
    p4 = getP4(chain)
    event.SOMEVARIABLE = 42*p4.M()
then:

sequence = [ someFunc, someOtherFunc, ...]
new_vars = [ SOMEVARIABLE, SOMEOTHERVARIABLE ] 

filter_module(infile, outfile, sequence, new_vars, cut, drop_branches, overwrite_output, tree_name=args.tree_name, n_max=n_max    )


then hopefull you will get your outfile which contains the new branches SOMEVARIABLE and SOMEOTHERVARIABLE

"""

#import sys
#sys.path.append("/afs/desy.de/user/n/nrad/analysis/commontools/FilterTools/") #need to do this in a nicer way



from filter_module import filter_module
from Helpers import *

import numpy as np



#track_names = ['track_tag_%s', 'pion1_signal_%s', 'pion2_signal_%s', 'lepton_signal_%s' ]

###
###
###
###
###

def isSignal(chain, event, *args, **kwargs):
    event.IS_SIGNAL = 1.0
    event.IS_BKG    = 0.0

def isBKG(chain, event, *args, **kwargs):
    event.IS_SIGNAL = 0.0
    event.IS_BKG    = 1.0

def mTauMinMax(chain, event, *args, **kwargs):
    p4s = getTrackP4s(chain, mode='3x1')
    P1 = p4s[0]
    P2 = sumP4s(p4s[1:])
    m_min, m_max = calcMTau(P1, P2)

    a, b, c = calcMTau(P1, P2, ret='coef')
    disc    = b*b-4*a*c

    event.mmin = m_min
    event.mmax = m_max
    event.disc = disc
    event.coef_1 = a
    event.coef_2 = b
    event.coef_3 = c

def invMelectron(chain, event, *args, **kwargs):
    p4s = get3prongP4s(chain, p4_vars=p4_vars, mass=electron_mass )
    p12 = p4s[0] + p4s[1]
    p23 = p4s[1] + p4s[2]
    p13 = p4s[0] + p4s[2]

    event.invM12 = p12.M()
    event.invM23 = p23.M()
    event.invM13 = p13.M()
    


get3prongInvMs
###
###
###
###
###

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("--tree_name", default="tau3x1")
    parser.add_argument("--cut", default='')
    parser.add_argument("--n", type=int, default=-1)
    parser.add_argument("--signal",        action='store_true')
    parser.add_argument("--overwrite",     action='store_true')
    parser.add_argument("--drop_branches", action='store_true')

    args = parser.parse_args()

    infile  = args.infile
    outfile = args.outfile

    overwrite_output = args.overwrite
    drop_branches    = args.drop_branches
    cut   = args.cut
    n_max = args.n

    #new_vars = ['mmin', 'mmax', 'IS_SIGNAL', 'IS_BKG', 'disc', 'coef_1', 'coef_2', 'coef_3', 'bdt']
    #sequence = [mTauMinMax] 

    new_vars = ['invM12', 'invM23', 'invM13']
    sequence = [invMelectron]

    if args.signal:
        sequence.append(isSignal)
    else:
        sequence.append(isBKG)


    filter_module(infile, outfile, sequence, new_vars, cut, drop_branches, overwrite_output, tree_name=args.tree_name, n_max=n_max    )

