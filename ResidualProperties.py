import numpy as np
from scipy.misc import derivative

"""
Currently, this file holds function that use the second-virial coefficient
correlation to calculate the residual properties of a species
Author: Eric Hou
Version: 1.0
"""


def B_0(Tr):
    return 0.083 - (0.422 / (Tr ** 1.6))


def B_1(Tr):
    return 0.139 - (0.172 / (Tr ** 4.2))


def dB_0(Tr):
    return 0.675 / (Tr ** 2.6)


def dB_1(Tr):
    return 0.722 / (Tr ** 5.2)


def res_Enthalpy(P, T, Pc, Tc, w, R=8.314):
    """
    This function uses the generalized second-virial-coefficient correlation
    to calculate the residual enthalpy of a compound at a certain state
    :param P: Pressure in desired units
    :param T: Temperature in absolute units
    :param Pc: Critical pressure in desired units
    :param Tc: Critical temperature in absolute units
    :param w: Acentric factor
    :param R: Universal gas constant (8.314 J / mol K)
    :return: Residual enthalpy in desired units
    """
    Pr = P / Pc
    Tr = T / Tc
    temp1 = B_0(Tr) - (Tr*dB_0(Tr))
    temp2 = B_1(Tr) - (Tr*dB_1(Tr))
    print(temp1, temp2)

    return R * Tc * Pr * (temp1 + (w*temp2))


def res_Entropy(P, T, Pc, Tc, w, R=8.314):
    """
    This function uses the generalized second-virial-coefficient correlation
    to calculate the residual entropy of a compound at a certain state
    :param P: Pressure in desired units
    :param T: Temperature in absolute units
    :param Pc: Critical pressure in desired units
    :param Tc: Critical temperature in absolute units
    :param w: Acentric factor
    :param R: Universal gas constant (8.314 J / mol K)
    :return: Residual entropy in desired units
    """
    Tr = T / Tc
    dB0 = dB_0(Tr)
    dB1 = dB_1(Tr)
    return R * -(P / Pc) * (dB0 + w * dB1)