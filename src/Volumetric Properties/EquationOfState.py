import numpy as np
from scipy.optimize import root

"""
This file provides functions that solve for the compressibility factor 
of gases and liquids under various conditions using common equations of state.
These equations include van der Waals, Peng-Robinson, Redlich-Kwong, and Soave-Redlich-Kwong
and are computed using a generic cubic equation of state. In addition, functions which solve  
for the same properties using the virial equations are given.
 
Author: Eric Hou
Version: 1.0
"""

# TODO: add Volume func given P,R,T,Z
# TODO: add set flag for each EOS so they can solve directly for V

def gen_Cubic_EOS_V(Z, b, q, o, e):
    """
    A helper method used for generic cubic equations of state in vapor form
    """
    return 1 + b - q * b * (Z - b) / ((Z + e * b) * (Z + o * b)) - Z


def gen_Cubic_EOS_L(Z, b, q, o, e):
    """
    A helper method used for generic cubic equations of state in liquid form
    """
    return b + (Z + e * b) * (Z + o * b) * (1 + b - Z) / (q * b) - Z


params = {
    #[alpha(Tr), sigma, epsilon, omega, phi, Zc]
    "PR": [(1 + np.sqrt(2)), (1 - np.sqrt(2)), 0.07780, 0.45724, 0.30740],
    "RK": [1, 0, 0.08664, 0.42748, (1 / 3)],
    "SRK": [1, 0, 0.08664, 0.42748, (1 / 3)],
    "vdW": [0, 0, (1 / 8), (27 / 64), (3 / 8)],
}


def PR_gas(P, T, Pc, Tc, w):
    """
    This function uses the parameters of the Peng-Robin equation to calculate the
    compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of the gas
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["PR"][2] * Pr / Tr
    q = params["PR"][3] * a / (params["PR"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["PR"][0], params["PR"][1]))
    return sol.x[0]


def PR_liquid(P, T, Pc, Tc, w):
    """"
    This function uses the parameters of the Peng-Robin equation to calculate the
    compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of liquid in bar
    :param T: Temperature of liquid in K
    :param Pc: Critical pressure of liquid in bar
    :param Tc: Critical temperature of liquid in K
    :param w: Acentric factor of the liquid
    :return: Compressibility factor of liquid (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["PR"][2] * Pr / Tr
    q = params["PR"][3] * a / (params["PR"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["PR"][0], params["PR"][1]))
    return sol.x[0]


def RK_gas(P, T, Pc, Tc, w):
    """
    This function uses the parameters of the Redlich-Kwong equation to calculate the
    compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of the gas
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = Tr ** (-0.5)
    b = params["RK"][2] * Pr / Tr
    q = params["RK"][3] * a / (params["RK"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["RK"][0], params["RK"][1]))
    return sol.x[0]


def RK_liquid(P, T, Pc, Tc, w):
    """"
    This function uses the parameters of the Redlich-Kwong equation to calculate the
    compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of liquid in bar
    :param T: Temperature of liquid in K
    :param Pc: Critical pressure of liquid in bar
    :param Tc: Critical temperature of liquid in K
    :param w: Acentric factor of the liquid
    :return: Compressibility factor of liquid (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = Tr ** (-0.5)
    b = params["RK"][2] * Pr / Tr
    q = params["RK"][3] * a / (params["RK"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["RK"][0], params["RK"][1]))
    return sol.x[0]


def SRK_gas(P, T, Pc, Tc, w):
    """
    This function uses the parameters of the Soave-Redlich-Kwong equation to calculate
    the compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of the gas
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.480 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["SRK"][2] * Pr / Tr
    q = params["SRK"][3] * a / (params["SRK"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["SRK"][0], params["SRK"][1]))
    return sol.x[0]


def SRK_liquid(P, T, Pc, Tc, w):
    """"
    This function uses the parameters of the Soave-Redlich-Kwong equation to calculate
    the compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of liquid in bar
    :param T: Temperature of liquid in K
    :param Pc: Critical pressure of liquid in bar
    :param Tc: Critical temperature of liquid in K
    :param w: Acentric factor of the liquid
    :return: Compressibility factor of liquid (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.480 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["SRK"][2] * Pr / Tr
    q = params["SRK"][3] * a / (params["SRK"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["SRK"][0], params["SRK"][1]))
    return sol.x[0]


def vdW_gas(P, T, Pc, Tc, w):
    """
    This function uses the parameters of the van der Waals equation to calculate
    the compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param Pc: Critical pressure of gas in bar
    :param Tc: Critical temperature of gas in K
    :param w: Acentric factor of the gas (unnecessary for vdW)
    :return: Compressibility factor of gas (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = 1
    b = params["vdW"][2] * Pr / Tr
    q = params["vdW"][3] * a / (params["vdW"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["vdW"][0], params["vdW"][1]))
    return sol.x[0]


def vdW_liquid(P, T, Pc, Tc, w):
    """"
    This function uses the parameters of the van der Waals equation to calculate
    the compressibilty factor of the gas in a generic cubic equation of state
    :param P: Pressure of liquid in bar
    :param T: Temperature of liquid in K
    :param Pc: Critical pressure of liquid in bar
    :param Tc: Critical temperature of liquid in K
    :param w: Acentric factor of the liquid (unnecessary for vdW)
    :return: Compressibility factor of liquid (Z)
    """
    Tr = T / Tc
    Pr = P / Pc
    a = 1
    b = params["vdW"][2] * Pr / Tr
    q = params["vdW"][3] * a / (params["vdW"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["vdW"][0], params["vdW"][1]))
    return sol.x[0]


def virial_Z(Vig, B, C=0, R=83.14):
    """
    This function calculates the compressibilty factor of the gas using the
    Virial equations given the coeffecients B, and C if necessary
    :param Vig: Volume of gas in ideal state
    :param B: Virial coefficient B
    :param C: Virial coefficient C
    :param R: Universal gas constant (83.14 bar cm^3 / mol K)
    :return: Compressibility factor of gas (Z)
    """
    Z = 1 + (B / Vig) + (C / Vig ** 2)
    return Z


def virial_V(P, T, B, C=0, R=83.14):
    """
    This function solves for the molar volume of the gas using the virial equations
    given the coeffecients B, and C if necessary
    :param P: Pressure of gas in bar
    :param T: Temperature of gas in K
    :param B: Virial coefficient B
    :param C: Virial coefficient C
    :param R: Universal gas constant (83.14 bar cm^3 / mol K)
    :return: Volume of gas in units given by R
    """
    def solver(V):
        return (R * T * (1 + (B / V) + (C / V ** 2)) / P) - V

    V0 = R * T / P
    sol = root(solver, [V0])
    return sol.x[0]
