from PythonTools.SparseTools import SparseUp, np, arrayIsNone, getTHnSparse 
from PythonTools.SparseTools import default_kwargs
from PythonTools.standard import drawRatioPlot, th1Func, getHistMax
from PythonTools.Legend import fixLegend
from PythonTools.u_float import u_float
import ROOT

import time

default_intercept_cuts = {}
default_intercept_cuts.update(
                cut_track_pt=0.5,
                cut_n_track = 2,
                cut_pxd_intercept_dist = 0.1,
                cut_pxd_edge_dist = 0, 
                cut_pxd_in_intercept = True,
                cut_pxd_in_track = True,
                cut_intercept_inside = True,
                )

intercept_bins = np.array([
       ( 'z',      (250,   -5,    9) ),
       ( 'phi',    (200, -3.2,  3.2) ),
       ( 'u',      (200,   -0.7,    0.7) ),
       ( 'v',      (250,   -4,    4) ),
       ( 'layer',  (2,      1,    3) ),
       #( 'ladder', (8,      1,    9) ),
       ( 'ladder',  (12,      1,    13) ),
       ( 'sensor',  (2,      1,    3) ),
       ( 'inside',  (2,      1,    3) ),
       ( 'inside1', (2,      0,    2) ),
       ( 'inside2', (2,      0,    2) ),
       ( 'inside0p5', (2,      0,    2) ),

       ( 'dist',   (100,   0,    0.02) ),
       ( 'du',     (100,   -0.01,    0.01) ),
       ( 'dv',     (100,   -0.01,    0.01) ),
       ( 'sigu',   (100,   0,    0.005) ),
       ( 'sigv',   (100,   0,    0.01) ),
       #( 'nMatch', (5, 0, 5 ) ),
    ])
intercept_bins_dict = dict(intercept_bins)
    
    
class SparseUpPXD(SparseUp):
    ###
    ### functions for pxd eff
    ###
    intercept_bins = intercept_bins
    intercept_bins_dict = dict(intercept_bins)
    
    pxd_base = "PXDCluster"
    pxd_keys = [ "layer", "sensorID",  "intercept_idx",  "dist",  "inIntercept", "inTrack", "dist_to_edge" ]
    
    intercept_base = "Intercept"
    intercept_keys = [ "sensorID",  "z",  "phi",  "u",  "v",  "layer",  "ladder",  "sensor", 'pxd_idx']
    intercept_keys += ['du', 'dv', 'sigu', 'sigv', 'nMatch', 'dist', ]
    intercept_keys += ['inside', 'inside1', 'inside2', 'inside0p5']
    
    track_base = "RecoTracks"
    track_keys = [ "tan_lambda_estimate",  "phi0_estimate",  "pt_estimate", "n_pxd_hits", "first_pxd_layer", "last_pxd_layer" ]

    

    
    
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.passed_tot_dict = {}
    
    default_intercept_cuts = default_intercept_cuts
    #default_intercept_cuts.update(
    #            cut_track_pt=0.5,
    #            cut_n_track = 2,
    #            cut_pxd_intercept_dist = 0.1,
    #            cut_pxd_edge_dist = 0, 
    #            cut_pxd_in_intercept = True,
    #            cut_pxd_in_track = True,
    #            cut_intercept_inside = True,
    #)
    
    @default_kwargs(**default_intercept_cuts)
    def get_intercept_passed_tot(self, upr, **cuts):
        CUT_TRACK_PT            =   cuts['cut_track_pt']  
        CUT_N_TRACK             =   cuts['cut_n_track']  
        CUT_PXD_INTERCEPT_DIST  =   cuts['cut_pxd_intercept_dist']  
        CUT_PXD_EDGE_DIST       =   cuts['cut_pxd_edge_dist']  
        CUT_PXD_IN_INTERCEPT    =   cuts['cut_pxd_in_intercept']  
        CUT_PXD_IN_TRACK        =   cuts['cut_pxd_in_track']  
        CUT_INTERCEPT_INSIDE    =   cuts['cut_intercept_inside']  
        
        time_start = time.time()

        event_dict     = self.get_arrays(upr, ["nIntercept", "nRecoTracks", "nPXDCluster"])

        trks_dict        = self.get_arrays(upr, ['pt_estimate'], "RecoTracks" )
        trk_pt = trks_dict['pt_estimate']

        selected_evts  =      (event_dict['nRecoTracks'] == CUT_N_TRACK) \
                            & (trk_pt[trk_pt > CUT_TRACK_PT].counts == CUT_N_TRACK)

        ##
        ## Apply event cuts when getting arrays
        ##

        #pxd_dict       = self.get_arrays(upr,  self.pxd_keys,       self.pxd_base       , idx=selected_evts )

        intercept_dict = self.get_arrays(upr,  self.intercept_keys, self.intercept_base , idx=selected_evts )
        track_dict     = self.get_arrays(upr,  self.track_keys,     self.track_base     , idx=selected_evts )

        intercept_tot_idx =  (intercept_dict['inside']==1 if CUT_INTERCEPT_INSIDE else True)
        intercept_tot = {k:v[intercept_tot_idx].flatten() for k,v in intercept_dict.items()}

        #   ##
        #   ## Impose conditions on PXD Clusters ()
        #   ##

        #   pxd_selected_idx =    (pxd_dict['dist'] < CUT_PXD_INTERCEPT_DIST) \
        #                       & (pxd_dict['inIntercept'] == 1 if CUT_PXD_IN_INTERCEPT else True )\
        #                       & (pxd_dict['inTrack'] == 1 if CUT_PXD_IN_TRACK else True )\
        #                       & (pxd_dict['dist_to_edge'] >= CUT_PXD_EDGE_DIST )\


        #   intercept_passed_idx   = pxd_dict['intercept_idx'][pxd_selected_idx].astype(int)
        #   intercept_passed       = {k:v[intercept_passed_idx].flatten() for k,v in intercept_dict.items()}

        ##
        ## Impose conditions on PXD Clusters ()
        ##

        ## get intercepts which are matched to a pxd cluster
        intercept_matched_idx =  (intercept_tot_idx)\
                              & (intercept_dict['pxd_idx']>=0)
        intercept_matched_pxd_idx = intercept_dict['pxd_idx'][intercept_matched_idx].astype('int')


        ## selelct the corresponding pxd clusters and mpose additional requirements 
        pxd_dict_all     = self.get_arrays(upr,  self.pxd_keys,       self.pxd_base       , idx=selected_evts )
        pxd_dict_matched           = {k:v[intercept_matched_pxd_idx] for k,v in pxd_dict_all.items() }

        pxd_selected_idx =    (pxd_dict_matched['dist']          <  CUT_PXD_INTERCEPT_DIST) \
                            & (pxd_dict_matched['inIntercept']  == 1 if CUT_PXD_IN_INTERCEPT else True )\
                            & (pxd_dict_matched['inTrack']      == 1 if CUT_PXD_IN_TRACK else True )\
                            & (pxd_dict_matched['dist_to_edge'] >= CUT_PXD_EDGE_DIST )\

        # select the intercepts which correspond to the passed pxd clusters... (better way to do this?!)
        intercept_passed_idx   = pxd_dict_matched['intercept_idx'][pxd_selected_idx].astype(int)
        intercept_passed       = {k:v[intercept_passed_idx].flatten() for k,v in intercept_dict.items()}

        
        print(f"Got intercepts. It took { round(time.time() - time_start , 3) } s") 
        
        ret = {'intercept_tot'   :intercept_tot,  #flat
               'intercept_passed':intercept_passed, #flat
               'pxd':pxd_dict_matched, 
               'pxd_idx': pxd_selected_idx,
            
              }
        return ret
        
        
        
    def get_thn(self, di, bins):
        thn = getTHnSpare(di, bins)
        return thn
    
    def get_intercept_thn(self, upr, passed=None, tot=None, bins=None, **cuts):
        
        
        if arrayIsNone(passed) or arrayIsNone(tot):
            print("will use default passed, tot")
            if self.passed_tot_dict:
                print("  already made, will reuse")
                pass
            else:
                print("  getting passed_tot_dict")
                self.passed_tot_dict = self.get_intercept_passed_tot(upr, **cuts)
                self.cuts = cuts
        #else:
        #    self.passed_tot_dict = {}

        passed = self.passed_tot_dict['intercept_passed'] if arrayIsNone(passed) else passed
        tot    = self.passed_tot_dict['intercept_tot'] if arrayIsNone(tot) else tot
       
        bins = self.intercept_bins if arrayIsNone(bins) else bins
        
        start = time.time()
        print("Getting thns")
        self.thn_tot    = getTHnSparse( tot   , bins, name="InterceptsTotal")
        self.thn_passed = getTHnSparse( passed, bins, name="InterceptsPassed")
        print(f"   took, {round(time.time()-start, 3)} s")
        

        
def getRatio(denom, nom, name="Ratio", opt=None):
    r=denom.Clone()
    if not opt:
        r.Divide(nom)
    else:
        r.Divide(r, nom, 1,1,"B")
    r.SetName(name)
    r.SetTitle(name)
    return r


def fix_up_thn(thn):
    axes_labels = dict( [ (i,x.GetTitle()) for i,x in enumerate( thn.GetListOfAxes() ) ] )
    for iax, axname in axes_labels.items():
        setattr(thn, axname, thn.GetAxis(iax) )
        setattr(thn, f"index_{axname}", iax)
    setattr(thn, 'axes_labels', axes_labels )
    setattr(thn, 'axes_idx',  {v:k for k,v in axes_labels.items() })    
    
nCont=100

def get_proj_z_phi(thn, layer=1, ladder=None, sensor=None, name=None):
    if not hasattr(thn, 'layer'):
        fix_up_thn(thn)
    select_sensor(thn, layer=layer, ladder=ladder, sensor=sensor)
    #thn.layer.SetRange(  *(() if not layer else  (layer,layer)) ) 
    #thn.ladder.SetRange( *(() if not ladder else (ladder,ladder)) ) 
    #thn.sensor.SetRange( *(() if not sensor else (sensor,sensor)) ) 
    p = thn.Projection(thn.index_z, thn.index_phi)
    if name:
        p.SetTitle(name)
        p.SetName(name)
        
    p.SetContour(nCont)
    select_sensor(thn) # reset to original ranges
    return p

def get_tag_from_cuts(**cuts):
    tag = "_".join([ f"{k.replace('cut_','')}_{v}" for k,v in cuts.items()] )
    tag = tag if tag else "default"
    tag = f"cuts_{tag}"
    return tag
    


def select_sensor(thn, layer=None, ladder=None, sensor=None):
    thn.layer.SetRange(  *(() if not layer else  (layer,layer)) ) 
    thn.ladder.SetRange( *(() if not ladder else (ladder,ladder)) ) 
    thn.sensor.SetRange( *(() if not sensor else (sensor,sensor)) ) 


def get_proj_pass_tot(thn_pass,thn_tot, axis='z', eff_range=(0.5,1), x_range=None, doDraw=True, title=None, title_tot=None, title_pass=None, pass_title=None, leg_pos=(0.18,0.77,0.8,0.92), **sens):
    
    fix_up_thn(thn_tot)
    fix_up_thn(thn_pass)
    
    select_sensor(thn_tot, **sens)
    select_sensor(thn_pass, **sens)
    
    
    htot  = thn_tot.Projection(  getattr(thn_tot,  f'index_{axis}') )
    hpass = thn_pass.Projection( getattr(thn_pass, f'index_{axis}') )
    htot.Sumw2()
    hpass.Sumw2()
    if x_range:
        htot.GetXaxis().SetRangeUser(*x_range)
        hpass.GetXaxis().SetRangeUser(*x_range)

    line = hpass.Clone()
    line.Reset()

    hpass.SetLineColor(ROOT.kRed)
    hpass.SetLineWidth(2)
    htot.SetMaximum( getHistMax(htot)[1]*1.2 ) 
    htot.SetLineWidth(2)
 
    try:
        eff = ROOT.TEfficiency(hpass, htot)
        if not eff:
            raise Exception("TEfficiency Failed")
        c=ROOT.TCanvas()
        eff.Draw()
        c.Draw()
        heff = eff.GetPaintedGraph()
        if not heff:
            raise Exception("Something went wrong with getting TGraph from TEfficiency")
    except Exception as exp:
        print(e)
        eff = getRatio(hpass, htot)
        heff = eff

    heff.GetYaxis().SetRangeUser( *eff_range )

    #eff.SetLineColor(ROOT.kBlack)
    #eff.GetYaxis().SetTitle("Ratio")
    #eff.SetLineStyle(9)

    htot.drawOption = "hist"
    hpass.drawOption = "hist"

    eff.SetFillColor(40)
    eff.SetLineColor(40)
    eff.SetMarkerColor(46)
    eff.SetMarkerStyle(20)
    eff.SetMarkerSize(0.5)
    eff.drawOption = "5"
    eff2= eff.Clone()
    eff2.drawOption = "pX0"

    line.Reset()
    line = th1Func(line, lambda x, bc: 1 )
    line.SetName("line")
    line.SetTitle(" ")
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)
    line.SetLineStyle(3)
    line.SetStats(False)
    line.GetYaxis().SetTitle("Efficiency")
    line.GetYaxis().SetRangeUser( *eff_range )
    line.drawOption="hist X0 "

    #eff.Sumw2()=9
    #eff.drawOption="pE"
    #eff.SetMarkerSize(1)
    #eff.SetMarkerStyle(20)
    leg = ROOT.TLegend(*leg_pos)
    leg.AddEntry( hpass, title_pass if title_pass else hpass.GetTitle())
    leg.AddEntry( htot, title_tot if title_tot else htot.GetTitle())
    leg.SetMargin(0.1)
    fixLegend(leg)


    if title:
        htot.SetTitle(title)

    canvs = drawRatioPlot( [htot, hpass, leg], [line, eff, eff2], widths={'y_width':300,'x_width':800, 'y_ratio_width':300},
                        y_title_offset=1.5,
                        x_title_offset=3.0,
                        y_label_offset=0.02,
                        fix_overlapping_labels=False,
                        doDraw=doDraw
                     )
    #leg.Draw()

    select_sensor(thn_tot)
    select_sensor(thn_pass)
    return [ (htot,hpass,eff,line, eff2, heff),  canvs, leg ]
    
    
    
    #htot.Draw('hist')
    #hpass.Draw('hist same')
    #hpass.SetLineColor(ROOT.kRed)

    #h = th1Func(hpass, lambda x, bc : 100 if hpass.GetBinContent(x) > htot.GetBinContent(x) else 0 )
    #h.SetLineColor(ROOT.kGreen)
    #h.Draw('hist same')
    
    

