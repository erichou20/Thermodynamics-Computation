import numpy as np
from ResidualProperties import B_0, B_1

"""
This file gives funcitons to compute fugacity calculations using generalized correlations.
This includes pitzer correlations for pure species and mixtures
Author: Eric Hou
Version: 1.0
"""


def pitzer_Fugacity(P, T, Pc, Tc, w, coeff=False):
    """
    This function computes fugacity or fugacity coefficients for a pure species using the
    generalized pitzer correlation.
    :param P: Pressure of substance in desired units
    :param T: Temperature of substance in desired units
    :param Pc: Critical pressure of substance in desired units
    :param Tc: Critical temperature of substance in desired units
    :param w: Acentric factor of the substance
    :param coeff: Parameter to return the fugacity coefficient (Default: false)
    :return: Fugacity of the substance in given units of pressure
    """
    Pr = P / Pc
    Tr = T / Tc
    B0 = B_0(Tr)
    B1 = B_1(Tr)

    phi = np.exp((Pr / Tr) * (B0 + (B1 * w)))
    if (coeff):
        return phi
    return phi * P


def wij(wi, wj):
    """
    This functions computes the effective acentric factor of a mixture of pure gases
    to calculate the fugacity using the virial equation
    :param wi: Acentric factor of species i
    :param wj: Acentric factor of species j
    :return: Effective acentric factor
    """
    return (wi + wj) / 2


def Tcij(Tci, Tcj, kij=0):
    """
    This function computes the effective critical temperature of a mixture of pure gases
    to calculate the fugacity using the virial equation
    :param Tci: Critical temperature of species i in desired units
    :param Tcj: Critical temperature of species j in desired units
    :param kij: Empirical interaction parater (0 for like species)
    :return: Effective critical temperature in desired units
    """
    return np.sqrt(Tci*Tcj)*(1-kij)

def Zcij(Zci, Zcj):
    """
    This function computes the effective critical compressibility factor of a mixture of pure
    gases to calculate the fugacity using the virial equation
    :param Zci: Critical compressibility factor of species i
    :param Zcj: Critical compressibility factor of species j
    :return: Effective critical compressibility factor
    """
    return (Zci+Zcj)/2


def Vcij(Vci, Vcj):
    """
    This function computes the effective critical volume of a mixture of pure gases
    to calculate the fugacity using the virial equation
    :param Vci: Critical volume of species i in desired units
    :param Vcj: Critical volume of species j in desired units
    :return: Effective critical volume in desired units
    """
    return ((Vci**(1/3)+Vcj**(1/3))/2)**3
