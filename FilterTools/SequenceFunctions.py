"""
    Common functions which can be used as a sequence function in post processing
    functions have to take arugments f(chain, event, *args, **kwargs)
    chain is the input ROOT file
    event object is used to store the new variables 

"""
from Helpers import *





def isSignal(chain, event, *args, **kwargs):
    event.IS_SIGNAL = 1.0
    event.IS_BKG    = 0.0
    event.IS_DATA   = 0.0

def isBackground(chain, event, *args, **kwargs):
    event.IS_SIGNAL = 0.0
    event.IS_BKG    = 1.0
    event.IS_DATA   = 0.0

def isData(chain, event, *args, **kwargs):
    event.IS_SIGNAL = 0.0
    event.IS_BKG    = 0.0
    event.IS_DATA   = 1.0

def mTauMinMax(chain, event, *args, **kwargs):
    """
        Calculates mmin and mmax according to the CLEO Method
        also stores the coefficients of the quadratic equation
        and the imaginary solutions

        kwarg['mode'] can be '1x1', '3x1' or a list of track_names like
        track_names = ['track_tag_%s', 'pion1_signal_%s', 'pion2_signal_%s', 'lepton_signal_%s' ] but only fo 3x1  

    """
    mode = kwargs.get('mode', '3x1')
    p4s = getTrackP4s(chain, mode=mode)
    P1 = p4s[0]
    P2 = sumP4s(p4s[1:])
    m_min, m_max = calcMTauSquaredCLEO(P1, P2)

    a, b, c = calcMTauSquaredCLEO(P1, P2, ret='coef')
    disc    = b*b-4*a*c

    #print(m_min, m_max)
    event.mmin2 = np.real(m_min)
    event.mmax2 = np.real(m_max)
    event.mmin2_imag = np.imag(m_min)
    event.mmax2_imag = np.imag(m_max)
    event.mmin2_mag  = np.absolute(m_min)
    event.mmax2_mag  = np.absolute(m_max)
    event.mmin       = np.sqrt(np.absolute(m_min))
    event.mmax       = np.sqrt(np.absolute(m_max))



    event.disc = disc
    event.coef_1 = a
    event.coef_2 = b
    event.coef_3 = c

    
    event.cosTheta_3prong = getCosTheta(P2, 1.777)
    event.cosTheta_1prong = getCosTheta(P1, 1.777)
    event.cosTheta_3prong_down = getCosTheta(P2, 1.776)
    event.cosTheta_1prong_down = getCosTheta(P1, 1.776)
    event.cosTheta_3prong_up = getCosTheta(P2, 1.778)
    event.cosTheta_1prong_up = getCosTheta(P1, 1.778)



    sqrt_disc  = sqrt(abs(disc)) 
    event.physical = disc >= 0

    event.mmin2_man  =  (-b - sqrt_disc) / (2*a)  if a else -999
    event.mmax2_man  =  (-b + sqrt_disc) / (2*a)  if a else -999


import functools
from scipy import optimize
def ms2(chain, event, *args, **kwargs):
    """

    """
    p4s = getTrackP4s(chain, p4_vars=p4_vars_CMS)
    p4_1pr = p4s[0]
    p4_3pr = p4s[1:]
    p4_a1  = sumP4s(p4_3pr)
    p4_miss = getPmiss(chain)

    f = functools.partial( ms2_func, p4_1pr, p4_a1, p4_miss ) 
    res = optimize.minimize(f, 
                 np.array([ [1.0],[1.0] ]),
                 method="SLSQP",
                 bounds = ((-4,4),(-4,4)),
                 )
    
    px_inv, py_inv = res.x
    ms2 = res.fun
    ms = math.sqrt(ms2) if ms2>=0 else 9999

    try:
        pz_inv  = math.sqrt( (Eb-p4_a1.E())**2 - px_inv*px_inv - py_inv*py_inv )
    except ValueError:
        pz_inv  = 9999

    event.inv_px = px_inv
    event.inv_py = py_inv
    event.inv_pz = pz_inv
    event.ms2 = ms2
    event.ms  = ms
    m2base_1prong, m2base_3prong = ms2_func(p4_1pr, p4_a1, p4_miss, (px_inv, py_inv), return_max=False)
    event.ms2_1prong = m2base_1prong
    event.ms2_3prong = m2base_3prong
    event.success = res.success
    #print(p4_1pr.M(), p4_a1.M(), ms2, ms, event.ms2, event.success, res)
    #"inv_px", "inv_py", "inv_pz", "ms2", "ms", "ms2_1prong", "ms2_3prong"
    #    #        pz_inv  = math.sqrt( (Eb-p4_3prong.E())**2 - px_inv*px_inv - py_inv*py_inv )



def tracksInvM(chain, event, *args, **kwargs):
    event.invM_2e_max = max(max(chain.tau_3prong_invM12_e,chain.tau_3prong_invM13_e),chain.tau_3prong_invM23_e)
    event.invM_2e_min = min(min(chain.tau_3prong_invM12_e,chain.tau_3prong_invM13_e),chain.tau_3prong_invM23_e)
    event.invM_2pi_max = max(max(chain.tau_3prong_invM12,chain.tau_3prong_invM13),chain.tau_3prong_invM23)
    event.invM_2pi_min = min(min(chain.tau_3prong_invM12,chain.tau_3prong_invM13),chain.tau_3prong_invM23)


def mMin(chain, event, *args, **kwargs):
    event.m_min_mc = getMmin(chain, p4_vars=p4_vars_MC, boost=cmsBoostVector)[0]  
    #event.m_min    = getMmin(chain, p4_vars=p4_vars_CMS, boost=cmsBoostVector)[0]
    event.m_min    = getMmin(chain, p4_vars=p4_vars_CMS)[0]


#rand = ROOT.TRandom(123456789)
rand = ROOT.TRandom()

def mMinSmearedReco(chain, event, *args, **kwargs):
    p4s = get3prongP4s(chain, p4_vars=p4_vars_CMS)
    resols     = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0] 
    mmin_names = ["m_min_"+str(sigma).replace(".","p") for sigma in resols]

    mmins = []
    for i, resol in enumerate(resols):
        smeared_p4s = [smearP4ByResol(x.Clone(), resol, rand=rand) for x in  p4s]
        assert len(smeared_p4s)==3
        smeared_a1 = sumP4s(smeared_p4s)
        smeared_mmin = getMminFromP4(smeared_a1)
        mmin_name = mmin_names[i]
        setattr(event, mmin_name, smeared_mmin)
    
def mMinSmeared(chain, event, *args, **kwargs):
    p4s_3x1 = getTrackP4s(chain, mode="3x1", p4_vars=p4_vars_MC, boostToCMS=True)
    p4s = p4s_3x1[1:] #
    resols     = [0.001, 0.005, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2]
    mmin_names = ["m_min_mc_"+str(sigma).replace(".","p") for sigma in resols]
     
    mmins = []
    for i, resol in enumerate(resols):
        smeared_p4s = [smearP4ByResol(x.Clone(), resol, rand=rand) for x in  p4s]
        assert len(smeared_p4s)==3
        smeared_a1 = sumP4s(smeared_p4s)
        smeared_mmin = getMminFromP4(smeared_a1)
        mmin_name = mmin_names[i]
        setattr(event, mmin_name, smeared_mmin)




def mMinSFs(chain, event, *args, **kwargs):
   
    sf_name = "bucket16-proc12-pi_2" 
    p4s = get3prongP4s(chain, p4_vars=p4_vars_CMS)

    sf_list       = [ findCosThetaSF(p4, sf_dict=sf_dicts[sf_name], key='central'  ) for p4 in p4s]
    sf_list_up    = [ findCosThetaSF(p4, sf_dict=sf_dicts[sf_name], key='up'       ) for p4 in p4s]
    sf_list_down  = [ findCosThetaSF(p4, sf_dict=sf_dicts[sf_name], key='down'     ) for p4 in p4s]
   
    #print(sf_list)
    #print(sf_list_up)
    #print(sf_list_down)
 
    event.m_min_sf        = getMminCorrected(p4s, sf_list)
    event.m_min_sf_up     = getMminCorrected(p4s, sf_list_up)
    event.m_min_sf_down   = getMminCorrected(p4s, sf_list_down)

def evt(chain, event, *args, **kwargs):
    event.evt = chain.__event__

def mMinSFsCharged(chain, event, *args, **kwargs):
    
    p4s = get3prongP4s(chain, p4_vars=p4_vars_CMS)
    charges = get3prongVar(chain, 'charge')
    sf_charge_dict = {-1:sf_dicts['bucket16-proc12-neg'], 1:sf_dicts['bucket16-proc12-pos']}
    sf_dict_list = [ sf_charge_dict[charge] for charge in charges]

    sf_list       = [ findCosThetaSF(p4, sf_dict=sf_dict_list[i], key='central'  ) for i,p4 in enumerate(p4s)]
    sf_list_up    = [ findCosThetaSF(p4, sf_dict=sf_dict_list[i], key='up'       ) for i,p4 in enumerate(p4s)]
    sf_list_down  = [ findCosThetaSF(p4, sf_dict=sf_dict_list[i], key='down'     ) for i,p4 in enumerate(p4s)]
   
    #print(sf_list)
    #print(sf_list_up)
    #print(sf_list_down)
 
    event.m_min_sf        = getMminCorrected(p4s, sf_list)
    event.m_min_sf_up     = getMminCorrected(p4s, sf_list_up)
    event.m_min_sf_down   = getMminCorrected(p4s, sf_list_down)


def ArgusGKK(chain, event, *args, **kwargs):
    p4s = get3prongP4s(chain, p4_vars=p4_vars_CMS)

  
    n_theta  = 1000
    in_theta = np.linspace(0,2*np.pi,n_theta)

    p4 = sumP4s(p4s)
    for i_theta, theta in enumerate(in_theta):
        event.m_min_gkk[i_theta] = getMminGKKFromP4(p4, theta)


