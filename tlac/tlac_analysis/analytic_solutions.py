"""Analytic solutions to various Lya radiative transfer settings.
"""
import numpy as np
from physics import *

def neufeld_solution(tau0, xlst, T, xi = 0):
    """Neufeld slab solution.
    """
    atau = a(T) * tau0
    # Eq. (2.24) in Neufeld paper with tau_s = 0, x_i=0
    # nom = np.sqrt(6)/24.0 * xlst**2 / atau
    # denom = np.cosh( np.sqrt(np.pi**4.0/54.0) * (np.abs(xlst**3) / atau))
    # myJ = nom / denom

    # Eq. (3.51) in Laursen (2010)
    nom = np.sqrt(6)/24.0 * xlst**2 / (atau * np.sqrt(np.pi))
    denom = np.cosh( np.sqrt(np.pi**3.0/54.0) * (np.abs(xlst**3 - xi**3) / atau))
    myJ = nom / denom    
    return myJ


def neufeld_solution_spherical(tau0, xlst, T):
    """Neufeld spherical solution as computed by Dijkstra+06 (?)
    """
    atau = a(T) * tau0

    nom = np.sqrt(np.pi) * xlst**2
    denom = np.sqrt(24) * atau * (1 + np.cosh(np.sqrt(2 * np.pi**3/27.0) * 
                                              np.abs(xlst**3)/atau))

    return nom / denom



def pdf_upar(upar, T, x):
    """u_parallel pdf given a frequency `x` and a temperature `T`
    """
    const = np.pi**1.5/a(T) * voigt(x, 1/np.sqrt(2), a(T))
    other = np.exp(-upar**2)/((x - upar)**2 + a(T)**2)
    return other/const
