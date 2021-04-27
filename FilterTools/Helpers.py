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


def getP1P2(chain, base_name="track_1prong{iTrk}_%s", p4_vars=p4_vars_CMS):
    P1, P2 = get1prongP4s(chain, p4_vars=p4_vars, base_name=base_name )
    return P1,P2

def getCosTheta(P, mtau): #four momentum of the pion and the given mtau
    ptau = np.sqrt(Eb**2 - mtau**2)
    #return (2*Eb*P.E() - mtau**2 - pion_mass**2)/ (2*P.P()*ptau)
    cosTheta = (2*Eb*P.E() - mtau**2 - P.M2() )/ (2*P.P()*ptau)
    #print(cosTheta)
    return cosTheta


def calcMTau(P1, P2, ret='roots'): #from four momentum
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
            print('no solution', [a,b,c], P1.M(), P2.M() , q1, q2, p1_cross_p2_2, p1_dot_p2)
            return []    
    
        ret = np.array( [ (np.sqrt(x) if np.isreal(x) and np.real(x)>=0 else  - np.sqrt( np.real(x*np.conj(x)) ) ) for x in  roots ] , dtype='float64') ## if complex, still keep it as a crazy number
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

def get1prongP4s(chain, base_name="track_1prong{iTrk}_%s",  p4_vars=p4_vars_CMS, mass=pion_mass):
    p4s = [ getP4(chain, base_name=basee_name.format(iTrk=iTrk), p4_vars=p4_vars, mass=mass) for iTrk in [0,1] ]
    return p4s



def get3prongP4s(chain, p4_vars=p4_vars, mass=pion_mass):
    p4s = [ getP4(chain, base_name="track{iTrk}_3prong_%s".format(iTrk=iTrk), p4_vars=p4_vars, mass=mass) for iTrk in range(1,4)]
    return p4s

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


