import numpy as np
from scipy.optimize import root

"""
Author: Eric Hou
Version: 1.0
"""


# TODO: Carnot pumps, possibly P_given_dS, work lost, and reservoirs

def eff_Carnot(Tc, Th):
    """

    :param Tc: Temperature of cold reservoir in K
    :param Th: Temperature of hot reservoir in K
    :return: Efficiency of carnot cycle
    """
    if Th == 0:
        raise ZeroDivisionError('Temperature cannot be 0')
    elif Tc > Th:
        raise ValueError('Tc cannot be larger than Th')
    return 1 - (Tc / Th)


def deltaS_ig(T, T0, A, B, C, D, P=1, P0=1, R=8.314):
    """
    This function calculates the change in entropy of a thermodynamic
    process assuming a non-constant heat capacity
    :param T: Final temperature in K
    :param T0: Initial temperature in K
    :param A: Constant A
    :param B: Constant B
    :param C: Constant C
    :param D: Constant D
    :param P: Final pressure in desired units
    :param P0: Initial pressure in desired units
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Change in entropy in desired units
    """
    if T == T0:
        return 0
    elif T == 0 or T0 == 0:
        raise ZeroDivisionError("Temperature Cannot be Zero")
    ds = A * np.log(T / T0) + (B + (C + D / (T0 ** 2 * T ** 2)) * ((T + T0) / 2)) * (T - T0)
    ds = (ds - np.log(P / P0)) * R
    return ds


def T_given_deltaS_ig(dS, T0, A, B, C, D, P=1, P0=1, R=8.314):
    """
    This function calculates the final temperature of a process given the change
    in entropy of a certain species.
    :param dS: Change in entropy in desired units
    :param T0: Initial temperature in K
    :param A: Constant A
    :param B: Constant B
    :param C: Constant C
    :param D: Constant D
    :param P: Final pressure in desired units
    :param P0: Initial pressure in desired units
    :param R: Universal gas constant (8.314 J/mol K)
    :return: Final Temperature in K
    """

    def T_solver(T):
        return deltaS_ig(T, T0, A, B, C, D, P, P0, R) - dS

    # Shifts initial guess depending on dS
    if dS > 0:
        x0 = [T0 + 100]
    elif dS < 0:
        x0 = [T0 - 100]
    elif dS == 0 and P == P0:
        return T0
    else:
        x0 = [T0+1]
    sol = root(T_solver, [x0])
    return sol.x[0]


def P_given_deltaS_ig(dS, T, T0, A, B, C, D, P0, R=8.314):
    def P_solver(P):
        if dS == 0:
            return deltaS_ig(T, T0, A, B, C, D, P, P0, R)
        return deltaS_ig(T, T0, A, B, C, D, P, P0, R) - dS

    if dS < 0:
        x0 = [P0*2]
    elif dS > 0:
        x0 = [P0/2]
    elif dS == 0 and T == T0:
        return P0
    else:
        x0 = [100]
    sol = root(P_solver, [x0], method = 'krylov')
    return sol.x[0]
