import ROOT
from math import *
import os, sys, time
import itertools
import os
import glob
import functools 
import numpy as np


print("ROOT Version", ROOT.gROOT.GetVersion())


cut_phase3_3pi = 'thrust>0.9 && thrust<0.99 && nPi0s_3prong==0 && nPhotons_3prong == 0 && nPi0s_1prong<=1 && nPhotons_1prong == 0 && track1_3prong_EoverP<0.8 && track2_3prong_EoverP<0.8 && track3_3prong_EoverP<0.8 && track1_3prong_EoverP > 0 && track2_3prong_EoverP > 0 && track3_3prong_EoverP > 0 && track_1prong_EoverP > 0 && visibleEnergyOfEventCMS<10.2 && visibleEnergyOfEventCMS>2.5'

cuts = {
                #'pipi_signal'  :  "(tauPlusMCMode==3 && tauMinusMCMode==3)",  
                'phase3_3pi' : cut_phase3_3pi,
                'phase3_3pi_mmin': cut_phase3_3pi + "&& Mmin>=1.5 && Mmin <= 2.2",
                'pipi_signal'    :  '(tauPlusMCMode==3 || tauMinusMCMode==3) && (tauPlusMCMode==5 || tauMinusMCMode==5)',  
                'pipi_pi0_veto'  :  '(dmID_1prong0==211 && dmID_1prong1==211) && (nPi0s_1prong0==0 && nPi0s_1prong1==0)', 
                'no_cut' : '(1)', 
                'ten_percent' : '__event__%10==1', 
       }

new_sel = 'track1_3prong_pt>=0.6 && track2_3prong_pt>=0.2 && track3_3prong_pt>=0.1 && thrust>=0.87 && thrust<=0.97 && visibleEnergyOfEventCMS>=2.5 && visibleEnergyOfEventCMS<=9 && missingMomentumOfEventCMS_theta>=0.15 && missingMomentumOfEventCMS_theta<=2.9 && missingMomentumOfEventCMS>=0.05 && missingMomentumOfEventCMS<=3.5 &&  nPi0s_3prong<=0 && nPi0s_1prong<=1 && nPhotons_3prong<=0 && nPhotons_1prong<=0 && tau_3prong_invM_2e_min>=0.7 && tau_3prong_invM_2e_min<=1.5 && tau_3prong_invM_2e_max>=0.8 && tau_3prong_invM_2e_min<=1.5'

cuts = {
                'phase3_3pi' : 'thrust>0.9 && thrust<0.99 && nPi0s_3prong==0 && nPhotons_3prong == 0 && nPi0s_1prong<=1 && nPhotons_1prong == 0 && track1_3prong_EoverP<0.8 && track2_3prong_EoverP<0.8 && track3_3prong_EoverP<0.8 && track1_3prong_EoverP > 0 && track2_3prong_EoverP > 0 && track3_3prong_EoverP > 0 && track_1prong_EoverP > 0 && visibleEnergyOfEventCMS<10.2 && visibleEnergyOfEventCMS>2.5',
                'track_EoP' : ' track1_3prong_EoverP > 0 && track2_3prong_EoverP > 0 && track3_3prong_EoverP > 0 && track_1prong_EoverP > 0',
                'ten_percent' : '(__event__%10==0)',
                'no_cut': '(1.0)',
                'phase3_3pi_noNeutralSelection' : 'thrust>0.9 && thrust<0.99  && track1_3prong_EoverP<0.8 && track2_3prong_EoverP<0.8 && track3_3prong_EoverP<0.8 && track1_3prong_EoverP > 0 && track2_3prong_EoverP > 0 && track3_3prong_EoverP > 0 && track_1prong_EoverP > 0 && visibleEnergyOfEventCMS<10.2 && visibleEnergyOfEventCMS>2.5',

                'mmin_1p7_1p85' : 'Mmin>=1.7 && Mmin<=1.85',
                'mmin_1p5_2p5' : 'Mmin>=1.5 && Mmin<=2.5',
                'mmin_1p5' : 'Mmin>=1.5',
                'ten_percent_mmin_1p7_1p85' : '(__event__%10==0) && Mmin>=1.7 && Mmin<=1.85',

                ## new cuts without invM stuff... can't be applied since not in the tuples yet!
                'new_sel'              :  new_sel,
                'new_sel_mmin_1p5_2' : f'{new_sel} && Mmin>=1.5 && Mmin<=2',
            }




def filter_module(infile, outfile, sequence, new_vars, cut='', drop_branches=[], keep_branches=[], overwrite_output=False, tree_name='tau3x1', n_max=-1, new_vector_vars=[],   *args, **kwargs):
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
    #new_vector_vars = []   # [ ('vname', nmax), (...) ] 

    strstr = \
    """\
    struct newvars_t {\n\
        %s\n\
        %s\n\
    \n}
    """ % ( '\n'.join(['Double_t {v};'.format(v=v) for v in new_vars]),
            '\n'.join(['Double_t {v}[{nMax}];'.format(v=v, nMax=nMax) for v, nMax in new_vector_vars])
         )
    print(strstr)

    files = glob.glob(infile)
    input_tree = ROOT.TChain(tree_name)
    for f in files:
        input_tree.Add(f)
    if not len(files):
        sys.exit("Exiting: No files were found at: %s"%infile)

    ##################################################
    ###  Applying Selection on the input tree     ####
    ##################################################

    n_total = input_tree.GetEntries()
    print('nevents', n_total)
    cut = cuts.get(cut, cut) 
 
    if cut:
        input_tree.Draw(">>eList", cut)
        eList = ROOT.eList
        if eList:
            n_selected = eList.GetN()
            print('selected events', n_selected)
        else:
            print("*****   WARNING... no events survived!")
            n_selelcted = 0
    else:
        n_selected = n_total

    #print(f"Input files added. nfiles={len(files)}. first input file is {files[:1]}")
    g = ROOT.TFile(outfile,"recreate")

    ##################################################
    ###       filtering branches if needed        ####
    ##################################################
    
    if drop_branches or keep_branches:
        from fnmatch import fnmatch
        initial_branch_list = [ x.GetTitle().split("/")[0] for x in input_tree.GetListOfBranches()]
        print(f"   starting with {len(initial_branch_list)} branches" )

        print("       first dropping some, ", drop_branches)
        branch_list = []
        branch_list_bl = []
        for branch in initial_branch_list:
            keep = True
            for pattern in drop_branches:
                if fnmatch(branch, pattern):
                    print("                    dropping branch", branch)
                    keep = False
                    if not pattern in ["*"]:
                        branch_list_bl.append(branch)
            if keep:
                branch_list.append(branch)

        print(f" now have {len(branch_list)}")
        print(f" now keeping some")
        for branch in initial_branch_list:
            keep = False
            for pattern in keep_branches:
                if fnmatch(branch, pattern):
                    #if branch in branch_list_bl:
                    #    print("           ", branch, "passes both balck and white lists... will drop it!")
                    #else:
                    #    keep = True
                    keep = True
            if keep and not branch in branch_list:
                branch_list.append(branch)
        print(f" Done with branch list! {len(branch_list)}")

        input_tree.SetBranchStatus("*", 0)
        for branch in branch_list:
            print(branch)
            input_tree.SetBranchStatus(branch, 1)
        input_tree.SetBranchStatus("__*__", 1)
        new_tree = input_tree.CloneTree(0)
        #input_tree.SetBranchStatus("*", 1)
    else:
        new_tree = input_tree.CloneTree(0)
    
    
    print ('adding new branches')
    ROOT.gROOT.ProcessLine(strstr)
    event = ROOT.newvars_t()
    for v in new_vars:
        new_tree.Branch(v, ROOT.AddressOf(event, v), f"{v}/D")
    for v, nMax in new_vector_vars:
        new_tree.Branch(v, ROOT.AddressOf(event, v), f"{v}[{nMax}]/D")
    
    
    if n_max > 0 and n_selected > n_max:
        n_selected = n_max
        print("\n WARNING: WILL ONLY RUN OVER %s EVENTS \n"%n_selected)
    
    print("---- starting the event loop...sit tight!\n")
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
    from SequenceFunctions import *    


    #DEFAULT_DROP_BRANCHES = ["tau_", "track*", "tau_*"]
    #DEFAULT_KEEP_BRANCHES = ["track*prong_E", "track*prong_pt*", "track*prong_px*", "track*prong_py*", "track*prong_pz*", "track*prong_EoverP", "track*prong_pionID", "track*prong_mc*", "tau_*prong_conv*", "tau_*prong_invM*", "tau_*prong_M*", "track*prong_*CMS", "track*prong*charge"]

    DEFAULT_DROP_BRANCHES = ["track*", "tau_*",]
    DEFAULT_KEEP_BRANCHES = ["track*prong_E", "track*prong_pt*", "track*prong_px*", "track*prong_py*", "track*prong_pz*", "track*prong_EoverP", "track*prong_pionID", "track*prong_mc*"]
    DEFAULT_KEEP_BRANCHES += ["tau_*prong_conv*", "tau_*prong_invM*", "tau_*prong_M*", "tau_*prong_p_CMS","tau_*prong_E_CMS", "invM*"]

    branch_list_for_fit = ['Mmin', 'nPhotons_1prong', 'track2_3prong_pt', 'nPi0s_1prong', 'tauPlusMCMode', 'track3_3prong_pt', 'invM_2e_min', 'nPhotons_3prong', 'tauMinusMCProng', 'visibleEnergyOfEventCMS', 'nPi0s_3prong', 'tauMinusMCMode', 'invM_2e_max', 'track1_3prong_pt', 'missingMomentumOfEventCMS', 'thrust', 'tauPlusMCProng', 'missingMomentumOfEventCMS_theta', 'lml*', 'ff*', ]
    branch_list_for_fit += ["track*_3prong_mcPX", "track*_3prong_mcPY", "track*_3prong_mcPZ", "track*_3prong_mcE"]
    DEFAULT_KEEP_BRANCHES = branch_list_for_fit
    DEFAULT_DROP_BRANCHES = ["*"]


    parser = argparse.ArgumentParser()
    #parser.add_argument("infile")
    #parser.add_argument("outfile")
    parser.add_argument("--input", default=[])
    parser.add_argument("--output")
    parser.add_argument("--tree_name", default="tau3x1")
    parser.add_argument("--cut", default='')
    parser.add_argument("--n", type=int, default=-1)
    parser.add_argument("--overwrite", default=False, action='store_true')
    parser.add_argument("--drop_branches", default=DEFAULT_DROP_BRANCHES, nargs="+", )
    parser.add_argument("--keep_branches", default=DEFAULT_KEEP_BRANCHES, nargs="+", )

    args = parser.parse_args()

    infile  = args.input
    outfile = args.output

    overwrite_output = args.overwrite
    drop_branches    = args.drop_branches
    keep_branches    = args.keep_branches
    cut   = args.cut
    n_max = args.n


    sequence = [tracksInvM, mMin]
    new_vars = [
                #'mmin2', 'mmin2_imag', 'mmin2_mag',
                #'mmax2', 'mmax2_imag', 'mmax2_mag',
                #'physical',   'IS_SIGNAL', 'IS_BKG', "IS_DATA",  'disc', 'coef_1', 'coef_2', 'coef_3',

                'm_min_mc', 'm_min',
       #         'mmin', 'mmax',
       #         'cosTheta_1prong', 'cosTheta_3prong',
       #         'cosTheta_1prong_up', 'cosTheta_3prong_up',
       #         'cosTheta_1prong_down', 'cosTheta_3prong_down',
    ]
    new_vars += [
                 'invM_2e_max',
                 'invM_2e_min',
                 'invM_2pi_max',
                 'invM_2pi_min',
    ]


    new_vector_vars = []

    filter_module(infile, outfile, sequence, new_vars, cut, drop_branches=drop_branches, keep_branches=keep_branches, overwrite_output=overwrite_output, tree_name=args.tree_name, n_max=n_max , new_vector_vars=new_vector_vars   )
