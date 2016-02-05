"""Calculate physical quantities which are sometimes handy.

All return values in cgs!
"""
import scipy.constants as C ## Careful: in SI!
import numpy as np

from scipy import special

## Some constants ###    
Lya_F12 = 0.4162
Delta_nu_L = 9.936e7
nu0 = 2.466e15
lambda0 = 1215.67
## cgs
echarge = 4.803205e-10
me = C.electron_mass * 1000
mH  = (C.m_p + C.m_e) * 1000
c = C.c * 100.
kb = 1.38065e-16
pc = 3.08567758e18

def thermal_velocity(T):
    """Returns v_th or hydrogen in cm/s given a temperature `T` in K
    """
    cc = C.k / (C.m_p + C.m_e)  # k_b / m_hydrogen in m^2/s^2 / K
    return np.sqrt(2. * cc * T) * 100.


def hydrogen_mass():
    """Mass of hydrogen atom in g
    """
    return (C.m_p + C.m_e) * 1e3



def voigt(x, sigma, gamma):
    """
    Voigt profile

    V(x, sigma, gamma) = \int_-inf^inf G(x', sigma)L(x-x', gamma) dx'
    where G is the gauss distribution with sd sigma and L the lorentzian with
    shape parameter gamma

    Computed using
    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))
    """
    try:
        lx = len(x)
    except TypeError:
        lx = 1
    
    z = (x + 1j * np.ones(lx) * gamma)/(sigma * np.sqrt(2))
    return special.wofz(z).real/(sigma * np.sqrt(2 * np.pi))

def Delta_nu_D(T):
    """Returns the frequency shift due to hydrogen gas at temperature T
    """
    return nu0/c * np.sqrt(2 * kb * T/mH)


def a(T):
    """Returns ratio between natural line width and doppler width
    """
    return Delta_nu_L/(2.0 * Delta_nu_D(T))

    
def H(T, x):
    """pi * Delta_nu_D * voigt

    Note: is not voigt_function_approx in tlac
    voigt_function_approx = H * Delta_nu_D
    """
    mya = a(T)
    return np.sqrt(np.pi) * Delta_nu_D(T) * voigt(x, 1/np.sqrt(2), mya)


def sigma_HI(T, x):
    """Hydrogen cross section
    """
    c1 = Lya_F12 * np.sqrt(np.pi) * echarge**2 / (me * c)
    return c1 * H(T, x) / Delta_nu_D(T)**2



def redshift_velocity(T, x):
    """ Returns the to the frequency corresponding velocity in km/s
    Keyword Arguments:
    T -- Temperature in K
    x -- frequency (dimensionless)
    """
    D = Delta_nu_D(T)
    return - x * D / ( x * D + nu0) * C.c / 1e3



def xcw(a):
    """Seperation between core and wing according to Laursen+10
    """
    return 1.59 - 0.6 * np.log(a) - 0.03 * np.log(a)**2


def x_to_v(x, T):
    """Converts frequency units x to v in cm/s
    """
    v_th = thermal_velocity(T)
    return - c * v_th * x / (v_th * x + c)

def v_to_lambda(v, z):
    """Converts frequency units from v (in cm/s) to lambda (in Angstrom)
    """
    return (v / c + 1.) * lambda0 * (z + 1)
    
    
