import numpy as np

"""
This file provides functions that estimates the latent heats of vaporization using
tools such as :
Author: Eric Hou
Version: 1.0
"""


#TODO: Clapeyron Derivation


def troutons(Tn, R=8.314):
    """
    This function uses Trouton's to estimate the heat of vaporization of a species
    given the normal boiling point. This is a rough estimation
    :param Tn: Normal boiling point in K
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Heat of vaporization in desired units
    """
    return R * Tn * 10


def riedel(Tc, Pc, Tn, R=8.314):
    """
    This function uses Riedel's Equation to estimate the heat of vaporization of a species
    :param Tc: Critical Temperature in K
    :param Pc: Critical Pressure in bar
    :param Tn: Normal boiling point in K
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Heat of vaporization in desired units
    """
    Trn = Tn / Tc
    return R * Tn * (1.092 * (np.log(Pc) - 1.013) / (0.93 - Trn))


def watson(deltaH1, T1, T2, Tc):
    """
    This function uses Watson's method to estimate the heat of vaporization of any liquid
    at any temperature given another known value at a known temperature.
    :param deltaH1: Known heat of vaporization at T1 in desired units
    :param T1: Known temperature in K
    :param T2: Desired temperature in K
    :param Tc: Critical temperature of species in K
    :return: Latent heat of vaporization at T2 in desired units
    """
    Tr1 = T1/Tc
    Tr2 = T2/Tc
    if Tr1 > Tc or Tr2 > Tc:
        raise TypeError("Temperature cannot be above critical temperature")
    return deltaH1*((1-Tr2)/(1-Tr1))**0.38

