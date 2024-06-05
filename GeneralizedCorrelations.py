import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import root
from scipy.integrate import quad

"""
This file provides functions that solve for the compressibility factor of gases 
and liquids under various conditions using some generalized correlations. These include the 
Pitzer correlation for the second and third virial coefficients, the rackett equation,
and the Lyderson-Greenkorn-Hougen correlation.
Author: Eric Hou
Version: 1.0
"""


def pitz_Cor_2_Virial(P, T, Pc, Tc, w, flag=False):
    """
    This function solves for the compressibility factor of a gas
    using the Pitzer correlations for the second Virial coefficient
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of gas
    :param flag: if set to true, returns the B_hat value
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    B_hat = (0.083 - (.422 / (Tr ** 1.6))) + w * (0.139 - (.172 / (Tr ** 4.2)))

    if flag:
        return B_hat
    return 1 + B_hat * Pr / Tr


def pitz_Cor_3_Virial(P, T, Pc, Tc, w):
    """
    This function solves for the compressibility factor of a gas
    using the Pitzer correlations for the third Virial coefficient
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of gas
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    B_hat = pitz_Cor_2_Virial(P, T, Pc, Tc, w, True)
    C_hat = 0.01407 + 0.02432 / Tr + .00313 * Tr ** -10.5 + w * (
                -0.02676 + (0.05539 * Tr ** -2.7) - 0.00242 * Tr ** -10.5)

    def solver(Z):
        return 1 + B_hat * Pr / (Tr * Z) + C_hat * (Pr / (Tr * Z)) ** 2 - Z

    sol = root(solver, [1])
    return sol.x[0]


def rackett(Vc, Zc, T, Tc):
    """
    This function estimates the molar volume of a saturated liquid using
    the Rackett equation
    :param Vc: Critical Volume of gas in desired units
    :param Zc: Critical compressibility factor of gas
    :param T: Temperature of gas in K
    :param Tc: Critical temperature of gas in K
    :return: Molar volume of saturated liquid in desired units
    """
    Tr = T / Tc
    return Vc * Zc ** (1 - Tr) ** (2 / 7)


def LGH(rho1, rho2, V1):
    """
    This function caclulates the molar volume of a liquid given the two parameters
    using the Lyderson-Greenkorn-Hougen generalized correlation
    :param rho2:
    :param V1:
    :return:
    """
    return V1 * rho1 / rho2
