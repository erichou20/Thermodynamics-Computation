from scipy.optimize import root

"""
This file gives functions to calculate heat changes assuming a non-constant heat capacity
Author: Eric Hou
Version: 1.0
"""


#TODO: Add entropy change functions

def Cp_ig(T, A, B, C, D, R=8.314):
    """
    This function calculates Cp_ig of a certain species at a specified temperature
    :param T: Temperature in Kelvin
    :param A: Constant A
    :param B: Constant B
    :param C: Constant C
    :param D: Constant D
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Cp_ig in desired units
    """
    if T == 0:
        raise ZeroDivisionError("Temperature cannot be zero")
    return R * (A + B * T + C * T ** 2 + D * T ** -2)


def Cv_ig(Cp_ig, R=8.314):
    """
    This function calculates Cv_ig of a certain species given Cp_ig
    :param Cp_ig: Cp_ig in desired units
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Cv_ig in desired units
    """
    return Cp_ig - R


def deltaH(T, T0, A, B, C, D, R=8.314):
    """
    This function calculates the change in enthalpy of a certain species from
    one temperature to another
    :param T: Final temperature in Kelvin
    :param T0: Initial temperature in Kelvin
    :param A: Constant A
    :param B: Constant B
    :param C: Constant C
    :param D: Constant D
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Change in enthalpy in desired units
    """
    return R * (A * (T - T0) + (B / 2) * (T ** 2 - T0 ** 2) +
                (C / 3) * (T ** 3 - T0 ** 3) + D * ((T - T0) / (T * T0)))


def T_Given_deltaH(dH, T0, A, B, C, D, R=8.314):
    """
    This function calculates the final temperature of a process given the change
    in enthalpy of a certain species.
    :param dH: Change in enthalpy in desired units
    :param T0: Initial temperature in Kelvin
    :param A: Constant A
    :param B: Constant B
    :param C: Constant C
    :param D: Constant D
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Final temperature in Kelvin
    """

    def solver(T):
        return deltaH(T, T0, A, B, C, D, R) - dH

    sol = root(solver, [T0 + 100])
    return sol.x[0]
