import ROOT
from math import *
import os, sys, time
import itertools
import os
import glob
import functools 
import numpy as np


print("ROOT Version", ROOT.gROOT.GetVersion())




cuts = {
                #'pipi_signal'  :  "(tauPlusMCMode==3 && tauMinusMCMode==3)",  
                'pipi_signal'  :  '(tauPlusMCMode==3 || tauMinusMCMode==3) && (tauPlusMCMode==5 || tauMinusMCMode==5)',  
                'pipi_pi0_veto'  :  '(dmID_1prong0==211 && dmID_1prong1==211) && (nPi0s_1prong0==0 && nPi0s_1prong1==0)', 
                'no_cut' : '(1)', 
                'ten_percent' : '__event__%10==1', 
       }


def filter_module(infile, outfile, sequence, new_vars, cut='', drop_branches=False, overwrite_output=False, tree_name='tau3x1', n_max=-1, *args, **kwargs):
    """
            sequence should be a list of eventFuncs which takes the "event" structure as input
            attributes of event corresponding to the new variables can be set, as in:

            

            def eventFunc(chain, event, *args, **kwargs)
                event.lepton_pt = chain.lep_px**2 + chain.lep_py**@

    """


    print('processing:\n input:', infile, '\n output:', outfile)

    startM = time.clock()
    startC = time.time()
    
    assert outfile.endswith(".root")
    dirname = os.path.dirname(outfile)
    if dirname:
        os.makedirs(dirname, exist_ok=True)

    if not overwrite_output:
        if os.path.isfile(outfile) and os.path.getsize(outfile)>500:
            raise Exception("Output file already exists: \n %s"%outfile)
    
    #new_vars        = []
    new_vector_vars = []   # [ ('vname', nmax), (...) ] 

    strstr = "\
    struct newvars_t {\n\
        %s\n\
        %s\n\
    \n}"%( '\n'.join(['Double_t {v};'.format(v=v) for v in new_vars]),
           '\n'.join(['Double_t {v}[{nMax}];'.format(v=v, nMax=nMax) for v, nMax in new_vector_vars])
        )
    print(strstr)

    files = glob.glob(infile)
    input_tree = ROOT.TChain(tree_name)
    for f in files:
        input_tree.Add(f)

    #print(f"Input files added. nfiles={len(files)}. first input file is {files[:1]}")
    g = ROOT.TFile(outfile,"recreate")
    
    if drop_branches:
        raise NotImplementedError("Sorry... please fix me? thanks!")
        from branch_lists import vars_to_keep
        keepbranches = vars_to_keep 
        input_tree.SetBranchStatus("*", 0)
        for br in keepbranches:
            print (br)
            input_tree.SetBranchStatus(br, 1)
        new_tree = input_tree.CloneTree(0)
        input_tree.SetBranchStatus("*", 1)
    else:
        new_tree = input_tree.CloneTree(0)
    
    
    print ('adding new branches')
    ROOT.gROOT.ProcessLine(strstr)
    event = ROOT.newvars_t()
    for v in new_vars:
        new_tree.Branch(v, ROOT.AddressOf(event, v), f"{v}/D")
    for v, nMax in new_vector_vars:
        new_tree.Branch(v, ROOT.AddressOf(event, v), f"{v}[{nMax}]/D")
    
    n_total = input_tree.GetEntries()
    print('nevents', n_total)
   


    cut = cuts.get(cut, cut) 
 
    if cut:
        input_tree.Draw(">>eList", cut)
        eList = ROOT.eList
        n_selected = eList.GetN()
        print('selected events', n_selected)
    else:
        n_selected = n_total
    
    if n_max > 0 and n_selected > n_max:
        n_selected = n_max
        print("\n WARNING: WILL ONLY RUN OVER %s EVENTS \n"%n_selected)
    
    for irow in range(n_selected):
        i = irow if not cut else eList.GetEntry(irow)
    
        if(i%10000==0): print( i)
        input_tree.GetEntry(i)
   
        for eventFunc in sequence:
            eventFunc(input_tree, event,*args, **kwargs) 
    
        new_tree.Fill()
     
    g.cd()
    new_tree.Write()
    g.Close()
    
    endM = time.clock()
    endC = time.time()
    
    print ("-"*20)
    print ("machine:", (endM-startM), "wall clock:", (endC-startC))




if __name__ == '__main__':


    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("--tree_name", default="tau3x1")
    parser.add_argument("--cut", default='')
    parser.add_argument("--n", type=int, default=-1)
    parser.add_argument("--overwrite", default=False, action='store_true')
    parser.add_argument("--drop_branches", default=False, action='store_true')

    args = parser.parse_args()

    infile  = args.infile
    outfile = args.outfile

    overwrite_output = args.overwrite
    drop_branches    = args.drop_branches
    cut   = args.cut
    n_max = args.n

    new_vars = []
    sequence = []

    #filter_module(infile, outfile, sequence, new_vars, cut, drop_branches, overwrite_output, tree_name=args.tree_name, n_max=n_max    )
