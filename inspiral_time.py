import numpy as np

from scipy import special
from scipy import integrate
from scipy.integrate import ode
from scipy.optimize import brentq
pi = np.pi

import pdb
import matplotlib.pyplot as plt
import scipy.constants as ct

from scipy.integrate import odeint, solve_ivp
from scipy.special import jv #bessel function of the first kind
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import interp1d


Msun = 1.99e30
G = 6.67e-11
c = 3.0e8

def e_factor(e):
    return (1+e)**1.1954 / (1-e**2.)**1.5

def chirp_mass(m1,m2,z=0):
    return (1+z)*(m1*m2)**0.6 / (m1+m2)**0.2

def eta(m1,m2):
    return (m1*m2)/(m1+m2)**2

def deda_peters(a,e):
    num = 12*a*(1+(73./24)*e**2 + (37./96)*e**4)
    denom = 19*e*(1-e**2)*(1+(121./304)*e**2)
    return denom/num

def inspiral_time_peters(a0,e0,m1,m2,af=0):
    """
    Computes the inspiral time, in Gyr, for a binary
    a0 in Au, and masses in solar masses
    
    if different af is given, computes the time from a0,e0
    to that af
    
    for af=0, just returns inspiral time
    for af!=0, returns (t_insp,af,ef)
    """
    coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        print(e0,a0)
        if not af == 0:
            print("ERROR: doesn't work for circular binaries")
            return 0
        return a0**4 / (4*beta)
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)
    
    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0)
        r.integrate(af)
        if not r.successful():
            print("ERROR, Integrator failed!")
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)
    
    if af==0:
        return integral * (12./19.) * c0**4. / beta
    else:
        return (integral * (12./19.) * c0**4. / beta), af, eFinal
