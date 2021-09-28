import os, sys, time
import math
from math import *
import itertools
import functools 
import ROOT
import numpy as np



p4_vars_MC = ["mcPX", "mcPY", "mcPZ", "mcE"]
p4_vars_CMS   = ['px_CMS','py_CMS','pz_CMS','E_CMS']
p4_vars   = ['px','py','pz','E']

electron_mass = 0.00051099895000
pion_mass     = 0.1395701766014103


Eb = 10.58/2.
Etau = Eb
mnu = 0 


def getlistofbranches(filename):
    f = open(filename,'r')
    outlist = []
    for line in f:
        if ":" in line:
            ll = line.split(':')
            print(ll)
            if len(ll) == 3: 
                outlist.append(ll[1].strip())
        else:
            outlist.append(line.strip())
    f.close()
    return outlist


def Mmin( m,e,p,be):
    return math.sqrt(  m*m + 2*(be/2-e)*(e-p) ) 


def getMmin(chain, p4_vars=p4_vars,  boost=False):
    a1 = sumP4s(get3prongP4s(chain, p4_vars=p4_vars, ))
    if boost:
        a1.Boost(boost)
    mmin = getMminFromP4(a1)
    return mmin, a1

def getMminFromP4(p4):
    mmin = -999
    if (not p4.E()>0) or (not p4.M()>0):
        mmin = -999
    else:
        try:
            mmin = math.sqrt( p4.M()**2 + 2*(Eb-p4.E())*(p4.E()-p4.P()) )
        except ValueError:
            mmin = -999
    return mmin



def getMminGKKFromP4(p4, theta):
    mmin = -999
    if (not p4.E()>0) or (not p4.M()>0):
        mmin = -999
    else:
        try:
            mmin = math.sqrt( p4.M()**2 + 2*(Eb-p4.E())*(p4.E()-p4.P()*cos(theta)) )
        except ValueError:
            mmin = -999
    return mmin

##
## cos theta SF
## 

sfs = {  
        ( -99, -0.5)   : 1.00054,  
        (-0.5, 0.6  )  : 1.00097,  
        (0.6, 0.8   )  : 1.00022, 
        (0.8, 99  )    : 0.99885,
      }

sfs_down = {  
        ( -99, -0.5)   : 1.00042,  
        (-0.5, 0.6  )  : 1.00094,  
        (0.6, 0.8   )  : 1.00017, 
        (0.8, 99  )    : 0.99878,
      }

sfs_up = {  
        ( -99, -0.5)   : 1.00066,  
        (-0.5, 0.6  )  : 1.00099,  
        (0.6, 0.8   )  : 1.00027, 
        (0.8, 99  )    : 0.99892,
      }


sf_dicts = { 
              #  'bucket16-proc12-neg':{
              #   ( -1.1, -0.5) : {'central':1.00059, 'down':1.00043, 'up':1.00065},
              #   (-0.5, 0.6)   : {'central':1.00015, 'down':1.00014, 'up':1.00017},
              #   (0.6, 0.8)    : {'central':0.99944, 'down':0.99932, 'up':0.99946},
              #   (0.8, 1)      : {'central':0.99832, 'down':0.99784, 'up':0.99835},
              #                     },
            'proc12':{
                   ( -1.1, -0.5)  : {'central':1.00046 , 'down': 1.00036,  'up':1.00055} ,
                   (-0.5, 0.6)    : {'central':1.00006 , 'down': 1.00003,  'up':1.00008} ,
                   (0.6, 0.8)     : {'central':0.99931 , 'down': 0.99926,  'up':0.99936} ,
                   (0.8, 1)       : {'central':0.99784 , 'down': 0.99776,  'up':0.99791} ,
                },
                'bucket16-proc12-pi':{
                 ( -1.1, -0.5) : {'central':1.00043,  'down':1.00039, 'up':1.00059},  
                 (-0.5, 0.6)   : {'central':1.00015,  'down':1.00014, 'up':1.00016},  
                 (0.6, 0.8)    : {'central':0.99932,  'down':0.99929, 'up':0.99944},  
                 (0.8, 1)      : {'central':0.99784,  'down':0.99780, 'up':0.99832},  
                                   },
                'bucket16-proc12-pi_2':{
                 ( -1.1, -0.5) : {'central':1.00043,  'down':1.00039, 'up':1.00059},  
                 (-0.5, 0.6)   : {'central':1.00015,  'down':1.00012, 'up':1.00027},  
                 (0.6, 0.8)    : {'central':0.99932,  'down':0.99929, 'up':0.99944},  
                 (0.8, 1)      : {'central':0.99784,  'down':0.99780, 'up':0.99832},  
                                   }
}

  
  
  




def findCosThetaSF(p4, sf_dict=sf_dicts['proc12'], key='central'):
    cosTheta = cos(p4.Theta())
    for (low,high), sfs in sf_dict.items():
        if low < cosTheta <= high:
            return sfs[key]
    return False

def applyMomentumSFonP4(p4, sf=1):
    scale_fact = [sf]*3 + [1]
    p4_corrected = ROOT.TLorentzVector()
    m = p4.M() 
    vect = p4.Vect()
    vect_corrected = vect*sf
    e_corrected    = math.sqrt(m*m + vect_corrected.Mag2() ) 
    p4_corrected.SetPxPyPzE(p4.Px()*sf, p4.Py()*sf, p4.Pz()*sf, e_corrected) 
    return p4_corrected


def getMminCorrected(p4s, sf_list=[1,1,1]):
    assert len(sf_list)==len(p4s)
    corrected_p4s = [ applyMomentumSFonP4(p4, sf_list[i]) for i, p4 in enumerate(p4s) ]
    corrected_mmin = getMminFromP4(sumP4s(corrected_p4s))
    return corrected_mmin




#######################################################################
#######################################################################

def getP1P2(chain, base_name="track_1prong{iTrk}_%s", p4_vars=p4_vars_CMS):
    P1, P2 = get1prongP4s(chain, p4_vars=p4_vars, base_name=base_name )
    return P1,P2

def getCosTheta(P, mtau): #four momentum of the pion and the given mtau
    """
        where theta is the angle between the track and the tau
    """
    ptau = np.sqrt(Eb**2 - mtau**2)
    #return (2*Eb*P.E() - mtau**2 - pion_mass**2)/ (2*P.P()*ptau)
    cosTheta = (2*Eb*P.E() - mtau**2 - P.M2() )/ np.abs(2*P.P()*ptau)
    #print(cosTheta)
    return cosTheta


def calcMTauSquaredCLEO(P1, P2, ret='roots'): #from four momentum
    """
            Use the analytical equations to find the mtau min
            ret :
                    'roots' : returns the solutions
                    'disc'  : returns the discriminant b^2-4ac
                    'coef'  : returns coefficients a,b,c
    
    """

    p1 = P1.Vect()
    p2 = P2.Vect()

    p1_dot_p2 = p1.Dot(p2)
    p1_cross_p2 = p1.Cross(p2)
    p1_cross_p2_2 = p1_cross_p2.Mag2()


    q1 = np.sqrt( 2*Etau*P1.E() - P1.M2() ) 
    q2 = np.sqrt( 2*Etau*P2.E() - P2.M2() ) 
    if np.isnan(q1) or np.isnan(q2):
        print('invalid q1 or q2', q1, q2)
        #return [float('nan'), float('nan')]
        return []

    p1_2 = p1.Mag2()
    p2_2 = p2.Mag2()


    a = p1_2 + p2_2 + 2*p1_dot_p2 
    b = -( 2*( p1_2*q2**2 + p2_2*q1**2 + p1_dot_p2*(q1**2+q2**2) - 2*p1_cross_p2_2 )  )
    c = p2_2*q1**4+p1_2*q2**4 + 2*p1_dot_p2*q1**2*q2**2-4*Eb**2*p1_cross_p2_2 + 4*mnu**2*p1_cross_p2_2
    if ret.lower() == 'coef':
        return (a, b, c)    
    elif ret.lower() == 'disc':
        return b*b-4*a*c
    elif ret.lower() == 'roots':

        try:
            roots = np.roots( [a,b,c] )
        except:
            #print(q1,q2)
            print('----------- WARNING no solution', [a,b,c], P1.M(), P2.M() , q1, q2, p1_cross_p2_2, p1_dot_p2)
            return [-999,-999]    
    
        #ret = np.array( [ (np.sqrt(x) if np.isreal(x) and np.real(x)>=0 else  - np.sqrt( np.real(x*np.conj(x)) ) ) for x in  roots ] , dtype='float64') ## if complex, still keep it as a crazy number
        #ret = np.array( [ (np.sqrt(x) if np.isreal(x) and np.real(x)>=0 else  - np.sqrt( np.real(x*np.conj(x)) ) ) for x in  roots ] , dtype='float64') ## if complex, still keep it as a crazy number
        ret  = np.array( roots )
        ret.sort()
        return ret
    else:
        raise ValueError("Ret option not recognized! %s"%ret)
    


p4_beam = ROOT.TLorentzVector()
p4_beam.SetPxPyPzE(0.4566179, 0, 2.9994152, 11.006)
p_beam = p4_beam.Vect()
e_beam = p4_beam.E()
cmsBoostVector = -ROOT.TVector3(p_beam.X()/e_beam, p_beam.Y()/e_beam, p_beam.Z()/e_beam);

def boostToFrame(p4, boost=cmsBoostVector, inplace=True):
    p4_ = p4.Clone() if not inplace else p4
    p4_.Boost(boost.X(), boost.Y(), boost.Z())
    return p4_



def getP4(chain,  base_name="track_1prong0_%s", p4_vars=p4_vars_CMS, mass=None):
    """ 
          p4_vars correspond to "px, py, pz, E"
          if mass given, the Energy will be recalculated accordingly 
    """

    var_names = [base_name%x for x in p4_vars]
    p4 = ROOT.TLorentzVector()
    p4.SetPxPyPzE(*[ getattr(chain, v) for v in var_names] )
    m = mass
    if m:
        new_E = math.sqrt( m*m + p4.P()*p4.P() ) #this is the corrected momentum already
        p4.SetE(new_E)
    return p4


def getVar(chain, var_names, base_name="track_1prong_%s", di=False):
    var_names = var_names if isinstance(var_names, (list,tuple)) else [var_names]
    ret = [getattr(chain, base_name%vname) for vname in var_names]
    if di:
        return dict(zip(var_names, ret))
    else:
        return ret[0] if len(var_names)==1 else ret    


def get1prongP4s(chain, base_name="track_1prong{iTrk}_%s",  p4_vars=p4_vars_CMS, mass=pion_mass):
    p4s = [ getP4(chain, base_name=basee_name.format(iTrk=iTrk), p4_vars=p4_vars, mass=mass) for iTrk in [0,1] ]
    return p4s



def get3prongP4s(chain, p4_vars=p4_vars, mass=None):
    p4s = [ getP4(chain, base_name="track{iTrk}_3prong_%s".format(iTrk=iTrk), p4_vars=p4_vars, mass=mass) for iTrk in range(1,4)]
    return p4s

def get3prongVar(chain, var):
    ret = [getVar(chain, var, base_name="track{iTrk}_3prong_%s".format(iTrk=iTrk)) for iTrk in range(1,4)]
    return ret

def get3prongInvMs(chain, p4_vars=p4_vars, mass=pion_mass):
    p4s = get3prongP4s(chain, p4_vars=p4_vars, mass=mass )
    p12 = p4s[0] + p4s[1]
    p23 = p4s[1] + p4s[2]
    p13 = p4s[0] + p4s[2]
    return (p12.M(), p23.M(), p13.M() )



def getTrackP4s(chain, mode='3x1', p4_vars=p4_vars_CMS, mass=None, boostToCMS=False):
    if mode == '3x1':
        base_names = ['track_1prong_%s'] + ['track{iTrk}_3prong_%s'.format(iTrk=iTrk)  for iTrk in [1,2,3] ]
    elif mode == '1x1':
        base_names = ['track_1prong{iTrk}_%s'.format(iTrk=iTrk)  for iTrk in [0,1] ]
    elif isinstance(mode, (list, tuple) ):
        base_names = mode
    else:
        raise Exception("Tau Decay Mode (%s) Not recognized!"%mode)

    p4s = []
    for base_name in base_names:
        p4s.append( getP4(chain, base_name, p4_vars, mass) )

    if boostToCMS:
        for p4 in p4s:
            boostToFrame(p4, boost=cmsBoostVector)

    return p4s



def sumP4s(p4s):
    tot = ROOT.TLorentzVector()
    for p4 in p4s:
        tot += p4
    return tot


