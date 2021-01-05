#from PythonTools.ROOTHelpers import *
#import standard
#from PythonTools.
#import PythonTools.standard as helpers 
import PythonTools.ROOTHelpers as ROOTHelpers
import math
import itertools
import uuid
import os
import ROOT



pxd_lay2_length = 61.440*1.E-3
pxd_lay2_length_small_pxl = 17.920*1.E-3
pxd_lay2_length_large_pxl = 43.520*1.E-3
pxd_lay1_length = 44.80*1.E-3
pxd_lay1_length_small_pxl = 14.080*1.E-3
pxd_lay1_length_large_pxl = 30.720*1.E-3
pxd_gap  = 0.85*1.E-3
pxd_u = 12.5*1E-3

pxd_lay2_v_total = 2*pxd_lay2_length + pxd_gap
pxd_lay1_v_total = 2*pxd_lay1_length + pxd_gap





binning_phi0 = (200,0,6.28)
binning_tan_lambda = (200,-2,4)
binning = binning_phi0+binning_tan_lambda

binning_phi0_coarse = (40,0,6.28)
binning_tan_lambda_coarse = (40,-2,4)
binning_coarse = binning_phi0_coarse+binning_tan_lambda_coarse

binning_phi0_overlap = (200,1,3)
binning_tan_lambda_overlap = (200,-2,4)
binning_overlap = binning_phi0_overlap+binning_tan_lambda_overlap



cut_first_pxd_layer = "(n_pxd_hits>0) & (first_pxd_layer==1)"
cut_only_first_pxd_layer = "(n_pxd_hits>0) & (last_pxd_layer==1)"
cut_second_pxd_layer = "(n_pxd_hits>0) & (last_pxd_layer==2)"
cut_only_second_pxd_layer = "(n_pxd_hits>0) & (first_pxd_layer==2)"
cut_any_pxd_layer = "(n_pxd_hits>0) & (first_pxd_layer>=1)"
cut_any_pxd_layer2 = "(n_pxd_hits>0) & (last_pxd_layer>=1)"
cut_both_pxd_layer = "(n_pxd_hits>0) & (last_pxd_layer==2) & (first_pxd_layer==1)"
cut_no_pxd_layer = "(n_pxd_hits==0)"

cut_any_svd_layer  = "last_svd_layer>=0"
cut_any_svd_layer2 = "first_svd_layer>=0"
cut_lotta_hits     = "n_svd_hits>=8 & n_cdc_hits>=40"


xvar="phi0_estimate"
yvar="tan_lambda_estimate"



scenarios = [
                {'cut': cut_first_pxd_layer        ,'name':'first_layer'       , 'title':'First Layer' },
                {'cut': cut_second_pxd_layer       ,'name':'second_layer'      , 'title':'Second Layer' },
                {'cut': cut_only_first_pxd_layer   ,'name':'only_first_layer'  , 'title':'Only First Layer' },
                {'cut': cut_only_second_pxd_layer  ,'name':'only_second_layer' , 'title':'Only Second Layer' },
                {'cut': cut_any_pxd_layer          ,'name':'any_layer'         , 'title':'Any Layer' },
                {'cut': cut_both_pxd_layer         ,'name':'both_layer'        , 'title':'Both Layers' },
                {'cut': cut_no_pxd_layer           ,'name':'no_layer'          , 'title':'No Layers' },
]


tracking_modes = ['default','no_pxd_layer_1','no_pxd_layer_2']


tracking_mode_opts = {
                        'default'  : {'name': "Default Tracking"},
                  'no_pxd_layer_1' : {'name': "No Layer 1 in Reco"},
                  'no_pxd_layer_2' : {'name': "No Layer 2 in Reco"},
        
    
}

runs_exp8 = [43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 59, 110, 111, 112, 114, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 132, 133, 155, 156, 157, 158, 161, 173, 174, 179, 181, 182, 183, 184, 185, 189, 190, 191, 192, 193, 194, 211, 212, 213, 214, 215, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 230, 231, 232, 234, 235, 236, 237, 238, 240, 241, 243, 244, 245, 246, 247, 276, 277, 278, 280, 282, 287, 288, 290, 291, 292, 294, 295, 296, 297, 299, 300, 301, 302, 303, 305, 306, 328, 329, 330, 331, 332, 333, 339, 340, 341, 342, 343, 344, 345, 347, 348, 349, 350, 351, 352, 355, 357, 358, 359, 360, 361, 362, 363, 364, 366, 367, 368, 369, 372, 373, 374, 375, 377, 378, 390, 393, 554, 555, 556, 557, 568, 569, 626, 631, 635, 641, 642, 647, 650, 654, 658, 659, 660, 665, 668, 671, 674, 676, 770, 771, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 790, 791, 792, 795, 796, 797, 798, 799, 825, 826, 827, 828, 829, 831, 832, 833, 834, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 1001, 1005, 1006, 1008, 1010, 1011, 1012, 1019, 1020, 1021, 1022, 1036, 1037, 1038, 1039, 1040, 1041, 1050, 1053, 1054, 1055, 1056, 1058, 1059, 1060, 1061, 1064, 1065, 1068, 1069, 1070, 1071, 1090, 1092, 1093, 1094, 1095, 1096, 1098, 1099, 1103, 1104, 1105, 1109, 1116, 1117, 1122, 1123, 1124, 1125, 1126, 1131, 1135, 1136, 1137, 1143, 1144, 1145, 1147, 1148, 1149, 1150, 1156, 1157, 1158, 1160, 1161, 1162, 1163, 1164, 1165, 1168, 1169, 1170, 1171, 1172, 1174, 1175, 1200, 1201, 1202, 1203, 1204, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1215, 1216, 1217, 1221, 1224, 1225, 1226, 1228, 1229, 1230, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1242, 1245, 1248, 1262, 1263, 1265, 1266, 1273, 1274, 1275, 1276, 1277, 1278, 1286, 1287, 1288, 1289, 1291, 1292, 1293, 1295, 1296, 1306, 1307, 1308, 1309, 1315, 1316, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1331, 1332, 1333, 1334, 1335, 1336, 1413, 1418, 1419, 1518, 1519, 1522, 1524, 1525, 1526, 1527, 1528, 1529, 1533, 1534, 1538, 1539, 1540, 1541, 1542, 1544, 1546, 1547, 1548, 1549, 1550, 1551, 1553, 1554, 1917, 1918, 1919, 1920, 1921, 1922, 1925, 1926, 1932, 1935, 1936, 1937, 1940, 1941, 1942, 1943, 1946, 1948, 1949, 1950, 1951, 1952, 1953, 1954, 1957, 1958, 1959, 1960, 1961, 1962, 1964, 1965, 1966, 1967, 1968, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1984, 1985, 1986, 1987, 1988, 1989, 1994, 1995, 1996, 2047, 2048, 2049, 2062, 2063, 2064, 2066, 2067, 2068, 2069, 2070, 2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080, 2081, 2082, 2084, 2085, 2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138, 2139, 2147, 2148, 2149, 2150, 2151, 2152, 2153, 2158, 2162, 2163, 2165, 2168, 2169, 2170, 2171, 2172, 2173, 2174, 2175, 2181, 2184, 2185, 2186, 2189, 2190, 2191, 2192, 2193, 2196, 2197, 2198, 2199, 2200, 2201, 2202, 2203, 2204, 2205, 2209, 2213, 2214, 2215, 2216, 2217, 2221, 2222, 2223, 2228, 2229, 2233, 2234, 
            2235, 2236, 2237, 2239, 2240, 2241, 2242, 2243, 2244, 2246, 2249, 2250, 2251, 2252, 2253, 2254, 2259, 2261, 2265, 2266, 2267, 2269, 2270, 2273, 2279, 2280, 2281, 2412, 2413, 2414, 2415, 2416, 2417, 2418, 2419, 2420, 2421, 2422, 2423, 2426, 2427, 2430, 2431, 2432, 2523, 2524, 2525, 2526, 2538, 2539, 2540, 2542, 2543, 2545, 2546, 2547, 2549, 2550, 2551, 2552, 2556, 2559, 2560, 2561, 2562, 2563, 2564, 2565, 2567, 2568, 2569, 2575, 2578, 2579, 2580, 2581, 2584, 2585, 2586, 2588, 2592, 2593, 2594, 2595, 2599, 2600, 2605, 2606, 2607, 2608, 2609, 2610, 2611, 2612, 2614, 2615, 2616, 2617, 2618, 2620, 2621, 2622, 2623, 2624, 2625, 2627, 2628, 2630, 2631, 2632, 2633, 2634, 2635, 2636, 2637, 2638, 2639, 2640, 2641, 2642, 2643, 2644, 2645, 2646, 2647, 2648, 2652, 2653, 2654, 2655, 2656, 2657, 2658, 2659, 2660, 2661, 2662, 2663, 2671, 2672, 2734, 2736, 2737, 2738, 2739, 2741, 2742, 2743, 2744, 2747, 2749, 2750, 2758, 2759, 2760, 2761, 2762, 2783, 2784, 2785, 2786, 2788, 2789, 2791, 2792, 2793, 2795, 2797, 2798, 2799, 2802, 2803, 2804, 2808, 2888, 2897, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2906, 2909, 2934, 2935, 2939, 2940, 2941, 2942, 2943, 2944, 2945, 2949, 2950, 2951, 2952, 2953, 2954, 2970, 2971, 2972, 2973, 2974, 3087, 3088, 3089, 3092, 3115, 3118, 3119, 3120, 3121, 3122, 3123]

run_ranges_10 = [
 (44, 367),
 (368, 1062),
 (1065, 1243),
 (1246, 1540),
 (1541, 2068),
 (2069, 2250),
 (2251, 2576),
 (2579, 2633),
 (2634, 2786),
 (2787, 3124),
 ]

run_ranges_2=[
                (44,2060),
                (2069,3124),
]

run_ranges_all = run_ranges_10 + run_ranges_2

#####  getting chains for PXD samples



# base_dir = "/group/belle/users/nrad/pxd_eff/"
base_dir = "/pnfs/desy.de/belle/local/user/nrad/PXD"
file_name_base = base_dir + "/skimmed_v2/git_hash=None/global_tag=data_reprocessing_proc9/experiment_number=8/tracking_mode={tracking_mode}/loose/trackLevelInformation.root"


cut_good_track = "n_cdc_hits >=10 && n_svd_hits>=6"


binning_u = (100,-1.2,1.2)
binning_v = (100,-2.5,2.5)
binning_uv = binning_u + binning_v
binning_uid = (250, 0, 250)
binning_vid = (756, 0, 756)

binning_res_u = (100,-0.05,0.05)
binning_res_v = (100,-0.05,0.05)
binning_res_uid = (40,-20,20)
binning_res_vid = (40,-20,20)

#binning_res_v = (100,-0.01*1E4,0.01*1E4)
#binning_res_u = (100,-0.01*1E4,0.01*1E4)
binning_sig_u = (100,0,0.05)
binning_sig_v = (100,0,0.05)

#binning_sig_u = (100,0,50)
#binning_sig_v = (100,0,50)

binning_res = (100,0,20)

sensorIDvar="layer*10000+ladder*100+sensor"

#res_u = "U Residual [#mum]"
#res_v = "V Residual [#mum]"
#res   =   "Residual [#mum]"
label_res_u = "U Residual [cm]"
label_res_v = "V Residual [cm]"
label_res_uid = "U Residual [pixel]"
label_res_vid = "V Residual [pixel]"
label_res   =   "Residual [cm]"

#sig_u = "Sigma U [#mum]"
#sig_v = "Sigma V [#mum]"
label_sig_u = "Sigma U [cm]"
label_sig_v = "Sigma V [cm]"

label_u     = "U [cm]"
label_v     = "V [cm]"

label_uid     = "U [pixel]"
label_vid     = "V [pixel]"


#cuts="abs(z0_estimate)<1 && abs(d0_estimate)<0.3 && n_cdc_hits >=10 && n_svd_hits>=6"
goodtrk="abs(z0_estimate)<1 && abs(d0_estimate)<0.3 && n_cdc_hits >=10 && n_svd_hits>=6"
#uts = "(1)"
#uts = "n_pxd_hits==0"

cutname,cut = ( "no_pxd_hits", "n_pxd_hits==0" )
cutname,cut = ( "goodTrk", goodtrk )
#cutname,cut = ("nocut", "(1)")

pxdResVars = {
                'pxd_resv_vs_v'     : {'prefix':"", 'var':["v","distv"] ,        'cut':cut,  'binning':binning_res_v+binning_v, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_res_v, 'ytitle':label_v  , 'nPXDs':1},
                'pxd_resu_vs_u'     : {'prefix':"", 'var':["u","distu"] ,        'cut':cut,  'binning':binning_res_u+binning_u, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_res_u, 'ytitle':label_u , 'nPXDs':1},
                'pxd_sigv_vs_v'     : {'prefix':"", 'var':["v","sigv"] ,         'cut':cut,  'binning':binning_sig_v+binning_v, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_sig_v, 'ytitle':label_v  , 'nPXDs':1},
                'pxd_sigu_vs_u'     : {'prefix':"", 'var':["u","sigu"] ,         'cut':cut,  'binning':binning_sig_u+binning_u, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_sig_u, 'ytitle':label_u , 'nPXDs':1},
                'pxd_resu_vs_resv'  : {'prefix':"", 'var':["distv","distu"] ,    'cut':cut,  'binning':binning_res_v+binning_res_u, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_res_u, 'ytitle':label_res_v , 'nPXDs':1},
                'pxd_resu_vs_sigu'  : {'prefix':"", 'var':["distu","sigu"]  ,    'cut':cut,  'binning':binning_sig_u+binning_res_u, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_sig_u, 'ytitle':label_res_u , 'nPXDs':1},
                'pxd_resv_vs_sigv'  : {'prefix':"", 'var':["distv","sigv"]  ,    'cut':cut,  'binning':binning_sig_v+binning_res_v, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_sig_v, 'ytitle':label_res_v , 'nPXDs':1},
                'pxd_res_vs_v'      : {'prefix':"", 'var':["distv","v"]  ,       'cut':cut,  'binning':binning_v+binning_res_v,     'sensorID': sensorIDvar , 'title':'', 'xtitle':label_v, 'ytitle':label_res , 'nPXDs':1},
                'pxd_res_vs_u'      : {'prefix':"", 'var':["distu","u"]  ,       'cut':cut,  'binning':binning_u+binning_res_u,     'sensorID': sensorIDvar , 'title':'', 'xtitle':label_u, 'ytitle':label_res , 'nPXDs':1},

                "pxd_resv_vs_v_id"  : {'prefix':"", 'var':["vid_clus","vid_clus-vid_fit"] ,        'cut':cut,  'binning':binning_res_v+binning_v, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_res_vid, 'ytitle':label_v  , 'nPXDs':1},
                "pxd_resu_vs_u_id"  : {'prefix':"", 'var':["uid_clus","uid_clus-uid_fit"] ,        'cut':cut,  'binning':binning_res_u+binning_u, 'sensorID': sensorIDvar , 'title':'', 'xtitle':label_res_uid, 'ytitle':label_u  , 'nPXDs':1},

             }


#chains = {}
icolor = 0
my_colors = [12,46,30,20, 40,41,42,46,49,30,33,38,20]
my_styles  = [2,5,9,6]
#####




def getChains():
    chains = {}
    icolor = 0
    for tm in tracking_modes:
        file_name = file_name_base.format(tracking_mode=tm)
        chains[tm] = ROOT.TChain("my_tree")
        chains[tm].Add(file_name)
        chains[tm].SetLineColor( my_colors[icolor] )
        icolor +=1
    return chains




####



###
def decorateHist(h,**d):
    xaxis = h.GetXaxis()
    yaxis = h.GetYaxis()
    xaxis.SetTitle(d.get("xtitle","#phi_{0}[rad]"))
    yaxis.SetTitle(d.get("ytitle","tan(#lambda)"))
    xaxis.SetTitleSize(0.07)
    yaxis.SetTitleSize(0.07)
    xaxis.SetLabelSize(0.05)
    yaxis.SetLabelSize(0.05)
    xaxis.SetTitleOffset(d.get('xTitleOffset',0.8))
    yaxis.SetTitleOffset(d.get('yTitleOffset',0.5))

        

def getPXDth2d(tree, name, var = "%s:%s"%(yvar,xvar), binning=binning_x_phi, cut="(1)", weight="(1)", title=None, decorateDict={}):  
    h = ROOTHelpers.getPlotFromChain( tree, var, binning, cut, weight=weight, name=name)
    decorateHist(h,**decorateDict)
    if title:
        h.SetTitle(title)
    return h


def getE1pMinusE1(e1p,e1, label=""):
    e1pe1 = th2Func( e1p, lambda x,y,z: e1p.GetBinContent(x,y) - e1.GetBinContent(x,y) 
                    if  ( e1p.GetBinContent(x,y) and e1.GetBinContent(x,y) ) else 0)
    e1pe1.GetZaxis().SetRangeUser(-0.1,0.1)
    e1pe1.SetTitle("#epsilon_{1}'-#epsilon_{1} %s"%label)
    e1pe1.GetZaxis().SetTitle("difference")
    return e1pe1

def makeExplicateFunc( hOnly1st, hOnly2nd, hBoth):
    def explicateEffs(x,y,bc):
        X = hOnly1st.GetBinContent(x,y) # first_layer
        Y = hOnly2nd.GetBinContent(x,y) # second_layer        
        Z = hBoth.GetBinContent(x,y) # both

        disc = (-1 + X - Y)**2 - 4*Y
        sqrt_disc = math.sqrt(disc) if disc>=0 else -999

        #e1  = 1/2 * (1 + X - sqrt_disc - Y) 
        #e2  = 1/2 * (1 - X - sqrt_disc + Y) 
        #e1p = 1/2 * (Z - X*Z + sqrt_disc* Z + Y*Z)/(Y) if Y else 0
        #ret = [[e1,e2,e1p],] 

        e1__  = 1/2 * (1 + X + sqrt_disc - Y) if sqrt_disc != -999 and Y and X else 0
        e2__  = 1/2 * (1 - X + sqrt_disc + Y) if sqrt_disc != -999 and Y and X else 0
        e1p__ = 1/2*Z*(1 - X - sqrt_disc + Y)/(Y) if sqrt_disc != -999 and Y and X else 0
        ret = [e1__,e2__,e1p__]
        return ret
    
    ret = [ th2Func( hOnly1st, lambda x,y,z: explicateEffs(x,y,z)[i] )  for i in range(3)]
    #if title:
    #    ret.SetTitle(title)
    
    return ret



def getE1p( hE2, hBoth, label):
    def func(x,y,bc):
        bc_both = hBoth.GetBinContent(x,y)
        bc_e2   = hE2.GetBinContent(x,y)
        if bc_both and bc_e2:
            ratio = bc_both/bc_e2
            if ratio >1:
                return 0
            else:
                return ratio
        else:
            return 0
    
    e1p = th2Func( hE2, func )
    hname = "e1prime" + label if label else ""
    e1p.GetZaxis().SetRangeUser(0,1)
    e1p.SetName(hname)
    e1p_ = hBoth.Clone(hname)
    e1p_.Divide(hE2)
    e1p_.GetZaxis().SetRangeUser(0,1)
    return [e1p,e1p_]



# def explicateEffs(x,y,bc):
#     Y = lay2Only.GetBinContent(x,y) # second_layer
#     X = lay1Only.GetBinContent(x,y) # first_layer
#     Z = lay1a2.GetBinContent(x,y) # both
#     
#     disc = (-1 + X - Y)**2 - 4*Y
#     sqrt_disc = math.sqrt(disc) if disc>=0 else -999
#     
#     #e1  = 1/2 * (1 + X - sqrt_disc - Y) 
#     #e2  = 1/2 * (1 - X - sqrt_disc + Y) 
#     #e1p = 1/2 * (Z - X*Z + sqrt_disc* Z + Y*Z)/(Y) if Y else 0
#     #ret = [[e1,e2,e1p],] 
#     
#     e1__  = 1/2 * (1 + X + sqrt_disc - Y) if sqrt_disc != -999 and Y and X else 0
#     e2__  = 1/2 * (1 - X + sqrt_disc + Y) if sqrt_disc != -999 and Y and X else 0
#     e1p__ = 1/2 * (Z - X*Z - sqrt_disc*Z + Y*Z)/(Y) if sqrt_disc != -999 and Y and X else 0
#     ret = [e1__,e2__,e1p__]
#     return ret







def decorateHist2(h, xrange=(1,3)):

    ROOT.gStyle.SetPadLeftMargin(0.2)
    ROOT.gStyle.SetPadRightMargin(0.2)
    ROOT.gStyle.SetTitleFontSize(0.1)
    ROOT.gStyle.SetTitleAlign(23)
    ROOT.gStyle.SetTitleFont(63)

    h.Draw("COLZ")
    h.Modify()
    
    xaxis = h.GetXaxis()
    yaxis = h.GetYaxis()
    zaxis = h.GetZaxis()
    #h.SetTitleFont(63,"xyz")
    #h.SetLabelFont(63,"xyz")
    xaxis.SetTitle("#phi_{0}[rad]")
    yaxis.SetTitle("tan(#lambda)")

    for axis in [xaxis,yaxis,zaxis]:
        #axis.SetTitleFont(63)
        #axis.SetLabelFont(63)
        axis.SetTitleFont(42)
        axis.SetLabelSize(0.08)
        axis.SetTitleSize(0.09)
        axis.SetTitleOffset(0.8)
    h.SetTitleFont(63,"t")
    #xaxis.SetTitleOffset(1.2)
    xaxis.SetTitleOffset(0.6)
    xaxis.SetLabelOffset(-0.025)
    zaxis.SetLabelSize(0.05)
    zaxis.SetTitleSize(0.06)
    zaxis.SetTitleOffset(1.2)
    #yaxis.SetTitleOffset(0.0)
    
    if xrange:
        h.GetXaxis().SetNdivisions(4,8,0)
        h.GetXaxis().SetRangeUser(1,3)
        
    
    palette = h.GetListOfFunctions().FindObject("palette")
    if not palette:
        h.Draw("Z")
        h.Modify()
        ROOT.gPad.Update()
        ROOT.gPad.Modify()
        palette = h.GetListOfFunctions().FindObject("palette")
        
        if not palette:
            
            print("failed to get color palette %s"%palette)
            return
    palette.SetX1NDC(0.81);
    palette.SetX2NDC(0.86);
    palette.SetY1NDC(0.145);
    palette.SetY2NDC(0.9);
    

    

    #h.GetXaxis().SetRangeUser(1.3,2.6)
##
##
##









def fixHistLabels(h, font=42, size=12):
    x = h.GetXaxis()
    y = h.GetYaxis()
    z = h.GetZaxis()
    for axis in [x,y,z]:
        axis.SetLabelFont(font)
        axis.SetLabelSize(size)
        axis.SetTitleFont(font)
        axis.SetTitleSize(size)
    h.SetTitleFont(font)
    h.SetTitleSize(font)

def getOnly(hThis, hThat, label="", title=""):
    """
        eff_this*(1-eff_that)
    """
    def func(x,y,bc):
        bc_this = hThis.GetBinContent(x,y)
        bc_that = hThat.GetBinContent(x,y)
        return bc_this*(1-bc_that)
    
    onlyThis = th2Func( hThis, func )
    if title: onlyThis.SetTitle(title)
    return onlyThis
    

__canvs__ = []

def_widths = (700,500)
def draw(plot,  path="", name="", width=(250,500) , dOpt="COLZ"):
    cname=uuid.uuid4().hex
    c = ROOT.TCanvas(cname,cname,*width)
    plot.Draw(dOpt)
    c.Draw()
    __canvs__.append(c)
    if  path and name:
        saveCanvas(c,path,name,["png","pdf"],[])
    return 



def getRunRangeTagAndCut(run_range=[], tag_template="_runs_%s_to_%s", run_var="run_number"):
    #tag_run_range = "_runs_%s_to_%s"%(run_range[0],run_range[1]) if run_range else ""
    tag_run_range = tag_template%(run_range[0],run_range[1]) if run_range else ""
    cut_run_range = f"{run_var}>=%s && {run_var}<=%s"%(run_range[0],run_range[1]) if run_range else "(1)"
    return tag_run_range, cut_run_range





binning_uID = (250,0,250)
binning_vID = (768,0,768)
binning_uvID  = binning_uID + binning_vID

binning_uDelta = (50,-50,50)
binning_vDelta = (50,-50,50)
binning_uvDelta = binning_uDelta + binning_vDelta
binning_delta  = (40,0,40)


res_in_u = "Residual in U [pixel]"
res_in_v = "Residual in V [pixel]"
resid = "Residual [pixel]"
pos_u = "U [pixel]"
pos_v = "V [pixel]"


pxdROIvars = {
                'pxd_res_vs_u'     : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_uStart","{prefix}_roi_uDelta"] ,     'binning':binning_uDelta + binning_uID, 'sensorID': "{prefix}_sensorID" , 'title':'', 'xtitle':res_in_u, 'ytitle':pos_u },
                'pxd_res_vs_v'     : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_vStart","{prefix}_roi_vDelta"] ,     'binning':binning_vDelta + binning_vID, 'sensorID': "{prefix}_sensorID" , 'title':'', 'xtitle':res_in_v, 'ytitle':pos_v },
                'pxd_res_tot_vs_u' : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_uStart","{prefix}_roi_delta"] ,      'binning':binning_delta + binning_uID, 'sensorID': "{prefix}_sensorID" , 'title':'', 'xtitle':resid, 'ytitle':pos_u },
                'pxd_res_tot_vs_v' : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_vStart","{prefix}_roi_delta"] ,      'binning':binning_delta + binning_vID, 'sensorID': "{prefix}_sensorID" , 'title':'', 'xtitle':resid, 'ytitle':pos_v },
                'pxd_res'          : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_roi_vDelta","{prefix}_roi_uDelta"] , 'binning':binning_uvDelta  , 'sensorID': "{prefix}_sensorID"            , 'title':'', 'xtitle':res_in_u, 'ytitle':res_in_v },
                'pxd_hits_in_rois' : {'prefix':"pxd_in_roi_%(iPXD)s", 'var':["{prefix}_vStart","{prefix}_uStart"] ,         'binning':binning_uvID     , 'sensorID': "{prefix}_sensorID"            , 'title':'', 'xtitle':pos_u, 'ytitle':pos_v },
                'pxd_hits_trks'    : {'prefix':"%(iPXD)s"           , 'var':["vStart_{prefix}","uStart_{prefix}"] ,         'binning':binning_uvID     , 'sensorID':'sensorID_{prefix}'           , 'title':'', 'xtitle':'', 'ytitle':'' },
             }




def makePXDModulePlots(chain, prefix="pxd_in_roi_%(iPXD)s", var=["{prefix}_roi_vDelta","{prefix}_roi_uDelta"] , cut="run_number <= 2500", binning=(200,-100,100,200,-100,100), sensorID= "{prefix}_sensorID", title="", xtitle="", ytitle="", nPXDs=1, iLayNum=1):

    from PythonTools import setstyle
    setstyle.setMyStyle()

    #ROOTHelpers.getAndSetEList(chain, cut ) 
    chain.SetEventList(0)
    eListBase = ROOTHelpers.retrieveEList(chain,cut)

    c=ROOT.TCanvas("c1","c1", 1400,800)
    c.SetLeftMargin(0.3)
    c.SetRightMargin(0.3)
    #c.Divide(8,2,0.0001, 0.0001)
    
    c, pads = ROOTHelpers.CanvasPartition(c, 8,2 )
    junks = [c,pads]
    
 
    iPlot = 0

    #var_ = var.format(prefix=prefix) if isinstance(var,str) else ":".join( [x.format(prefix=prefix) for x in var] )
    hists = []
    hist_dicts = []

    for iLadNum in range(1,9,1):
        for iSenNum in range(1,3,1):
            sensorID_ = "%i0%i0%i"%(iLayNum,iLadNum, iSenNum)
            sensorIDCut = "(%s==%s)"%(sensorID.format(prefix=prefix) , sensorID_)

            ## retrieve and set eList to chain.... gotta be careful to remove any previously applied eList
            chain.SetEventList(0)
            eListSensorID = ROOTHelpers.retrieveEList(chain, sensorIDCut )
            eListTot = eListBase.Clone().Intersect(eListSensorID)
            chain.SetEventList(eListTot)
            ## 

            iPad = (iSenNum%2-1)+ 2*iLadNum-1 
            #print(iPad)
            pad = pads[iPad]
            h_ = []
            for iPXD in range(nPXDs):

                var_ = var.format(prefix=prefix) if isinstance(var,str) else ":".join( [x.format(prefix=prefix) for x in var] )   
                #print(var_)
                var_ = var_%{'iPXD':iPXD}
                #print(var_)
                cutString = "%s &&(%s==%s)"%(cut, sensorID.format(prefix=prefix) , sensorID_)
                cutString = cutString%{'iPXD':iPXD}
                h_.append( ROOTHelpers.getPlotFromChain( chain, var_, binning, cutString, "(1)" ) )
                #h_.append( hist )
                h_[-1].SetDirectory(0)
            h = h_[0]
            for h__ in h_[1:]:
                h.Add(h__)
            hists.append(h)
            hist_dicts.append( {'hist':h, 'pad':pad, 'iPad':iPad, 'sensorID':sensorID_ ,'sensorIDList':[iLayNum,iLadNum,iSenNum]  } ),

    junks.append(hist_dicts)
    maxVal = max( [ hdict['hist'].GetMaximum() for hdict in hist_dicts] )
    minVal = min( [ hdict['hist'].GetMinimum() for hdict in hist_dicts] )
    #print(maxVal)

    for hdict in hist_dicts:
        h        = hdict['hist']
        pad      = hdict['pad']
        sensorID = hdict['sensorID']
        iLayNum,iLadNum,iSenNum = hdict['sensorIDList']
        pad.cd()
        pad.SetLogz(1)
        pad.Update()


        hFrame = h.Clone()
        hdict['frame'] = hFrame
        hists.append(hFrame)


        h.SetMaximum(maxVal)
        hFrame.SetMaximum(maxVal)
        h.SetMinimum(minVal*0.9)
        hFrame.SetMinimum(minVal*0.9)
    
        h.GetZaxis().SetRangeUser(minVal*0.9,maxVal)
        hFrame.GetZaxis().SetRangeUser(minVal*0.9,maxVal)


        #hFrame.Reset();
        dOpt = ("COLZ" if iLadNum==8 else "COL")
        hFrame.Draw(dOpt);
        xFactor = pads[0].GetAbsWNDC()/pad.GetAbsWNDC();
        yFactor = pads[0].GetAbsHNDC()/pad.GetAbsHNDC();

        hFrame.GetXaxis().ChangeLabel(1,-1,0.)  # remove the first label
        hFrame.GetXaxis().ChangeLabel(-1,-1,0.) # remove the last label

        
        
        #Format for y axis
        hFrame.GetYaxis().SetLabelFont(43);
        hFrame.GetYaxis().SetLabelSize(16);
        hFrame.GetYaxis().SetLabelOffset(0.02);
        hFrame.GetYaxis().SetTitleFont(43);
        hFrame.GetYaxis().SetTitleSize(16);
        hFrame.GetYaxis().SetTitleOffset(5);
        
        hFrame.GetYaxis().CenterTitle();
        hFrame.GetYaxis().SetNdivisions(505);
        
        #TICKS Y Axis
        hFrame.GetYaxis().SetTickLength(xFactor*0.04/yFactor);
        
        #Format for x axis
        hFrame.GetXaxis().SetLabelFont(43);
        hFrame.GetXaxis().SetLabelSize(16);
        hFrame.GetXaxis().SetLabelOffset(0.02);
        hFrame.GetXaxis().SetTitleFont(43);
        hFrame.GetXaxis().SetTitleSize(16);
        hFrame.GetXaxis().SetTitleOffset(5);
        hFrame.GetXaxis().CenterTitle();
        hFrame.GetXaxis().SetNdivisions(505);
        
        #TICKS X Axis
        hFrame.GetXaxis().SetTickLength(yFactor*0.06/xFactor);
        # #h.GetXaxis().SetLabelFont(63)
        # h.GetXaxis().SetTitleFont(63)
        # h.GetYaxis().SetLabelFont(63)
        # h.GetYaxis().SetTitleFont(63)

        # h.GetXaxis().SetLabelSize(0.1)
        # h.GetXaxis().LabelsOption("v")
        # #h.GetXaxis().SetTitleSize(12)
        # h.GetYaxis().SetLabelSize(12)
        # h.GetYaxis().SetTitleSize(12)

        # h.GetXaxis().SetLabelOffset(0.01)
        # h.GetXaxis().SetTitleOffset(3)
        # h.GetYaxis().SetLabelOffset(0)
        # h.GetYaxis().SetTitleOffset(3)

        #dOpt = "SAME" + ("COLZ" if iLadNum==8 else "COL")
        pad.Update()

        #h.Draw(dOpt )
        if iSenNum==1:
            mod_label = sensorID.replace("0",".")
            mod_label = ".".join( mod_label.split(".")[:2]+["X"] )
            mod_label_xloc =  (pad.GetUxmax()+pad.GetUxmin())/2 
            #mod_label_yloc =   pad.GetUymax()*(0.9 if iSenNum==1 else 0.1)
            mod_label_yloc =   pad.GetUymax() - ( pad.GetUymax() - pad.GetUymin() )*(-0.05 if iSenNum==1 else 0.93 )
            #print( mod_label, mod_label_xloc, mod_label_yloc, pad.GetUymin(), pad.GetUymax() )
            hdict['label']=ROOTHelpers.drawLatex(mod_label , mod_label_xloc, mod_label_yloc , size=0.1*xFactor  , align=22, font=42,setNDC=False) 


    c.cd()
    c.Update()
    xtitle_xloc   =  0.8
    xtitle_yloc   =  0.01

    ytitle_xloc   =  0.03
    ytitle_yloc   =  0.5

    p = ROOT.TPad("axis","axis",0,0,1,1)
    p.SetFillStyle(4000)
    p.Draw()
    p.cd()
    lx      =  ROOTHelpers.drawLatex(xtitle, xtitle_xloc, xtitle_yloc, size=0.04 , font=42 ) 
    ly      =  ROOTHelpers.drawLatex(ytitle, ytitle_xloc, ytitle_yloc, size=0.04 , font=42 , align=22, angle=90) 
    ly.Draw()


    junks.extend([lx,ly,p])

    p.Update()
    c.Update()
    c.cd()
    #return hists, c,pads, hist_dicts
    chain.SetEventList(0)
    return  junks

def makePXDModulePlotsRDF(chain, prefix="pxd_in_roi_%(iPXD)s", var=["{prefix}_roi_vDelta","{prefix}_roi_uDelta"] , cut="run_number <= 2500", binning=(200,-100,100,200,-100,100), sensorID= "{prefix}_sensorID", title="", xtitle="", ytitle="", nPXDs=1):

    from PythonTools import setstyle
    setstyle.setMyStyle()
    rdf = ROOT.RDataFrame(chain)
    rdf_ = rdf.Filter(cut)

    is1DPlot = isinstance(var,str) 
    vs_ = [var] if is1DPlot else var
    ops = ["*","+","/","-"]
    newvs_ = []
    for v in vs_:
        vname  = v
        if any([op in v for op in ops]):
            for op in ops:
                vname = vname.replace(op,"__")
            rdf_=rdf_.Define(vname,v)
        newvs_.append(vname)
    #print(newvs_)
    var = newvs_
    

    c=ROOT.TCanvas("c1","c1", 1400,800)
    c.SetLeftMargin(0.3)
    c.SetRightMargin(0.3)
    c, pads = ROOTHelpers.CanvasPartition(c, 8,2 )
    junks = [c,pads]
 
    iPlot = 0
    iLayNum=1

    hists = []
    hist_dicts = []

    for iLadNum in range(1,9,1):
        for iSenNum in range(1,3,1):
            sensorID_ = "%i0%i0%i"%(iLayNum,iLadNum, iSenNum)
            sensorIDCut = "(%s==%s)"%(sensorID.format(prefix=prefix) , sensorID_)
            iPad = (iSenNum%2-1)+ 2*iLadNum-1 
            pad = pads[iPad]
            h_ = []
            for iPXD in range(nPXDs):
                #print(sensorID_, "iPXD")
                var_ = var.format(prefix=prefix) if isinstance(var,str) else ":".join( [x.format(prefix=prefix) for x in var] )   
                var_ = var_%{'iPXD':iPXD}
                cutString = "%s &&(%s==%s)"%(cut, sensorID.format(prefix=prefix) , sensorID_)
                cutString = cutString%{'iPXD':iPXD}
                #h_.append( ROOTHelpers.getPlotFromChain( chain, var_, binning, cutString, "(1)" ) )
                histo2d_args = ( (sensorID_,sensorID_) +binning, var_.split(":")[0], var_.split(":")[1]  )
                #print(histo2d_args)
                h_.append( rdf_.Filter( sensorIDCut ).Histo2D( *histo2d_args )  )
                #h_.append( hist )
                #print("set dir")
                #h_[-1].SetDirectory(0)
                #print("set dir done")
            h = h_[0]
            #print("got hists")
            for h__ in h_[1:]:
                h.Add(h__)
            hists.append(h)
            hist_dicts.append( {'hist':h, 'pad':pad, 'iPad':iPad, 'sensorID':sensorID_ ,'sensorIDList':[iLayNum,iLadNum,iSenNum]  } ),

    junks.append(hist_dicts)
    maxVal = max( [ hdict['hist'].GetMaximum() for hdict in hist_dicts] )
    minVal = min( [ hdict['hist'].GetMinimum() for hdict in hist_dicts] )
    #print(maxVal)

    for hdict in hist_dicts:
        h        = hdict['hist']
        pad      = hdict['pad']
        sensorID = hdict['sensorID']
        iLayNum,iLadNum,iSenNum = hdict['sensorIDList']
        pad.cd()
        pad.SetLogz(1)
        pad.Update()


        hFrame = h.Clone()
        hdict['frame'] = hFrame
        hists.append(hFrame)


        h.SetMaximum(maxVal)
        hFrame.SetMaximum(maxVal)
        h.SetMinimum(minVal*0.9)
        hFrame.SetMinimum(minVal*0.9)
    
        h.GetZaxis().SetRangeUser(minVal*0.9,maxVal)
        hFrame.GetZaxis().SetRangeUser(minVal*0.9,maxVal)


        #hFrame.Reset();
        dOpt = ("COLZ" if iLadNum==8 else "COL")
        hFrame.Draw(dOpt);
        xFactor = pads[0].GetAbsWNDC()/pad.GetAbsWNDC();
        yFactor = pads[0].GetAbsHNDC()/pad.GetAbsHNDC();

        hFrame.GetXaxis().ChangeLabel(1,-1,0.)  # remove the first label
        hFrame.GetXaxis().ChangeLabel(-1,-1,0.) # remove the last label

        
        
        #Format for y axis
        hFrame.GetYaxis().SetLabelFont(43);
        hFrame.GetYaxis().SetLabelSize(16);
        hFrame.GetYaxis().SetLabelOffset(0.02);
        hFrame.GetYaxis().SetTitleFont(43);
        hFrame.GetYaxis().SetTitleSize(16);
        hFrame.GetYaxis().SetTitleOffset(5);
        
        hFrame.GetYaxis().CenterTitle();
        hFrame.GetYaxis().SetNdivisions(505);
        
        #TICKS Y Axis
        hFrame.GetYaxis().SetTickLength(xFactor*0.04/yFactor);
        
        #Format for x axis
        hFrame.GetXaxis().SetLabelFont(43);
        hFrame.GetXaxis().SetLabelSize(16);
        hFrame.GetXaxis().SetLabelOffset(0.02);
        hFrame.GetXaxis().SetTitleFont(43);
        hFrame.GetXaxis().SetTitleSize(16);
        hFrame.GetXaxis().SetTitleOffset(5);
        hFrame.GetXaxis().CenterTitle();
        hFrame.GetXaxis().SetNdivisions(505);
        
        #TICKS X Axis
        hFrame.GetXaxis().SetTickLength(yFactor*0.06/xFactor);
        # #h.GetXaxis().SetLabelFont(63)
        # h.GetXaxis().SetTitleFont(63)
        # h.GetYaxis().SetLabelFont(63)
        # h.GetYaxis().SetTitleFont(63)

        # h.GetXaxis().SetLabelSize(0.1)
        # h.GetXaxis().LabelsOption("v")
        # #h.GetXaxis().SetTitleSize(12)
        # h.GetYaxis().SetLabelSize(12)
        # h.GetYaxis().SetTitleSize(12)

        # h.GetXaxis().SetLabelOffset(0.01)
        # h.GetXaxis().SetTitleOffset(3)
        # h.GetYaxis().SetLabelOffset(0)
        # h.GetYaxis().SetTitleOffset(3)

        #dOpt = "SAME" + ("COLZ" if iLadNum==8 else "COL")
        pad.Update()

        #h.Draw(dOpt )
        if iSenNum==1:
            mod_label = sensorID.replace("0",".")
            mod_label = ".".join( mod_label.split(".")[:2]+["X"] )
            mod_label_xloc =  (pad.GetUxmax()+pad.GetUxmin())/2 
            #mod_label_yloc =   pad.GetUymax()*(0.9 if iSenNum==1 else 0.1)
            mod_label_yloc =   pad.GetUymax() - ( pad.GetUymax() - pad.GetUymin() )*(-0.05 if iSenNum==1 else 0.93 )
            #print( mod_label, mod_label_xloc, mod_label_yloc, pad.GetUymin(), pad.GetUymax() )
            hdict['label']=ROOTHelpers.drawLatex(mod_label , mod_label_xloc, mod_label_yloc , size=0.1*xFactor  , align=22, font=42,setNDC=False) 


    c.cd()
    c.Update()
    xtitle_xloc   =  0.8
    xtitle_yloc   =  0.01

    ytitle_xloc   =  0.03
    ytitle_yloc   =  0.5

    p = ROOT.TPad("axis","axis",0,0,1,1)
    p.SetFillStyle(4000)
    p.Draw()
    p.cd()
    lx      =  ROOTHelpers.drawLatex(xtitle, xtitle_xloc, xtitle_yloc, size=0.04 , font=42 ) 
    ly      =  ROOTHelpers.drawLatex(ytitle, ytitle_xloc, ytitle_yloc, size=0.04 , font=42 , align=22, angle=90) 
    ly.Draw()


    junks.extend([lx,ly,p])

    p.Update()
    c.Update()
    c.cd()
    #return hists, c,pads, hist_dicts
    chain.SetEventList(0)
    return  junks



if False:
    save_dir = f"/afs/desy.de/user/n/nrad/www/PXD/{proc_tag}/"
    for run_range in run_ranges_all:
        run_tag, run_cut = getRunRangeTagAndCut( run_range )
        for pxdvarname, pxdvardict in pxdROIvars.items():
            ret = makePXDModulePlots( **pxdvardict , cut=run_cut, nPXDs=5)
            save_dir_ = os.path.join( save_dir, run_tag ) 
            saveCanvas( ROOT.gPad, save_dir_ , pxdvarname)




####

###
### RDF And Module Plots
###



#def makePXDModulePlotsRDF(chain, prefix="pxd_in_roi_%(iPXD)s", var=["{prefix}_roi_vDelta","{prefix}_roi_uDelta"] , cut="run_number <= 2500", binning=(200,-100,100,200,-100,100), sensorID= "{prefix}_sensorID", title="", xtitle="", ytitle="", nPXDs=1):
#
#xtitle = "xtitle"
#ytitle = "ytitle"


#chain = chains['trks']
isTH2D = True


def getModuleRDF(rdf, cut, sensorIDvar="sensorID", iLayNum=2):
    """
    Returns a dict with an RDataFrame corresponding to each module
    """
    rdf_dict = {}

    #rdf = rdf if isinstance(rdf,ROOT.RDataFrame) else ROOT.RDataFrame(chain) 
    rdf_filtered = rdf.Filter(cut)
    for iLadNum in range(1,9,1):
        for iSenNum in range(1,3,1):
            sensorID_ = "%i0%i0%i"%(iLayNum,iLadNum, iSenNum)
            #sensorIDCut = "(%s==%s)"%(sensorIDvar.format(prefix=prefix) , sensorID_)
            sensorIDCut = "(%s==%s)"%(sensorIDvar, sensorID_)
            #hname = uniqueName("h%s"%sensorId_)
            #cutString = f"{cut} && {sensorIDCut}"
            rdf_module = rdf_filtered.Filter(sensorIDCut)
            rdf_dict[sensorID_] = rdf_module
    return rdf_dict


#def getTH2DFromRDF(rdf, var, binning, name=None, title=None):#, xtitle=None, ytitle=None): 
#    var_ = ":".join(var) if isinstance(var,(tuple,list)) else var
#    if ":" in var:
#        var_ = list(reversed(var_.split(":"))) #RDataFrame wants x:y
#    
#    name = name if name else uniqueName("Hist")
#    title = title if title else name
#    histo2d_args = ( (name,title) +binning ,) + tuple( var_.split(":") )
#    print(histo2d_args)
#    h = rdf.Histo2D(*histo2d_args)
#    #if xtitle: h.GetXaxis().SetTitle(xtitle)
#    #if ytitle: h.GetYaxis().SetTitle(ytitle)
#    return h
    
getTH2DFromRDF = ROOTHelpers.getTH2DFromRDF    



def drawPXDModulePlots( hist_dicts, xtitle='', ytitle='', logz=1 ):
    #from PythonTools import setstyle
    #setstyle.setMyStyle()
    
    junks = []
    c=ROOT.TCanvas("c1","c1", 1400,800)
    c.SetLeftMargin(0.3)
    c.SetRightMargin(0.3)
    c, pads = ROOTHelpers.CanvasPartition(c, 8,2 )

    junks.extend([c,pads])

    maxVal = max( [ h.GetMaximum() for h in hist_dicts.values()] )
    minVal = min( [ h.GetMinimum() for h in hist_dicts.values()] )

    for sensorID, h in hist_dicts.items():
        iLayNum,iLadNum,iSenNum = map(int, sensorID.split("0"))
        iPad = (iSenNum%2-1)+ 2*iLadNum-1

        pad  = pads[iPad]
        pad.cd()
        pad.SetLogz(logz)
        pad.Update()


        hFrame = h.Clone()
        hFrame.SetTitle("")
        junks.append(hFrame)


        for h_ in [h,hFrame]:
            h_.SetMaximum(maxVal)
            h_.SetMinimum(minVal*0.9)
            h_.GetZaxis().SetRangeUser(minVal*0.9,maxVal)


        dOpt = ("COLZ" if iLadNum==8 else "COL")
        hFrame.Draw(dOpt);
        xFactor = pads[0].GetAbsWNDC()/pad.GetAbsWNDC();
        yFactor = pads[0].GetAbsHNDC()/pad.GetAbsHNDC();

        hFrame.GetXaxis().ChangeLabel(1,-1,0.)  # remove the first label
        hFrame.GetXaxis().ChangeLabel(-1,-1,0.) # remove the last label

        #Format for y axis
        hFrame.GetYaxis().SetLabelFont(43);
        hFrame.GetYaxis().SetLabelSize(16);
        hFrame.GetYaxis().SetLabelOffset(0.02);
        hFrame.GetYaxis().SetTitleFont(43);
        hFrame.GetYaxis().SetTitleSize(16);
        hFrame.GetYaxis().SetTitleOffset(5);

        hFrame.GetYaxis().CenterTitle();
        hFrame.GetYaxis().SetNdivisions(505);

        #TICKS Y Axis
        hFrame.GetYaxis().SetTickLength(xFactor*0.04/yFactor);

        #Format for x axis
        hFrame.GetXaxis().SetLabelFont(43);
        hFrame.GetXaxis().SetLabelSize(16);
        hFrame.GetXaxis().SetLabelOffset(0.02);
        hFrame.GetXaxis().SetTitleFont(43);
        hFrame.GetXaxis().SetTitleSize(16);
        hFrame.GetXaxis().SetTitleOffset(5);
        hFrame.GetXaxis().CenterTitle();
        hFrame.GetXaxis().SetNdivisions(505);

        #TICKS X Axis
        hFrame.GetXaxis().SetTickLength(yFactor*0.06/xFactor);

        #dOpt = "SAME" + ("COLZ" if iLadNum==8 else "COL")
        pad.Update()

        #h.Draw(dOpt )
        if iSenNum==1:
            mod_label = sensorID.replace("0",".")
            mod_label = ".".join( mod_label.split(".")[:2]+["X"] )
            mod_label_xloc =  (pad.GetUxmax()+pad.GetUxmin())/2
            #mod_label_yloc =   pad.GetUymax()*(0.9 if iSenNum==1 else 0.1)
            mod_label_yloc =   pad.GetUymax() - ( pad.GetUymax() - pad.GetUymin() )*(-0.05 if iSenNum==1 else 0.93 )
            #print( mod_label, mod_label_xloc, mod_label_yloc, pad.GetUymin(), pad.GetUymax() )
            l = ROOTHelpers.drawLatex(mod_label , mod_label_xloc, mod_label_yloc , size=0.1*xFactor  , align=22, font=42,setNDC=False)
            junks.append(l)

    c.cd()
    c.Update()
    
    xtitle_xloc   =  0.8
    xtitle_yloc   =  0.01

    ytitle_xloc   =  0.03
    ytitle_yloc   =  0.5

    xtitle = xtitle if xtitle else h.GetXaxis().GetTitle()
    ytitle = ytitle if ytitle else h.GetYaxis().GetTitle()

    
    p = ROOT.TPad("axis","axis",0,0,1,1)
    p.SetFillStyle(4000)
    p.Draw()
    p.cd()
    lx      =  ROOTHelpers.drawLatex(xtitle, xtitle_xloc, xtitle_yloc, size=0.04 , font=42 )
    ly      =  ROOTHelpers.drawLatex(ytitle, ytitle_xloc, ytitle_yloc, size=0.04 , font=42 , align=22, angle=90)
    ly.Draw()


    junks.extend([lx,ly,p])
     
    p.Update()
    c.Update()
    c.cd()
    c.Draw()
    return junks



