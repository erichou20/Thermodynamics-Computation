import numpy as np
from ResidualProperties import B_0, B_1
from GeneralizedCorrelationsForFugacity import Tcij, wij, Vcij, Zcij

""""
This file contains the a function calculate fugacity of species in solutions
Author: erichou20
Version: 1.0
"""


# MAY NEED TO HAVE HELPER FUNCTIONS

def fugacity_Solution(P, t, Pc_arr, Tc_arr, w_arr, zc_arr, Vc_arr, y_arr,
                      coeff=False, R=83.14, k_arr=None):
    """
    This function calculates the fugacity and fugacity coefficients for multi-component
    gas mixtures of n species by using a generalized procedure including virial coefficients.
    :param P: Pressure of the system in desired units
    :param t: Temperature of the system in desired units
    :param Pc_arr: List or array of respective critical pressures in desired units
    :param Tc_arr: List or array of respective critical temperatures in desired units
    :param w_arr: List or array of respective acentric values
    :param zc_arr: List or array of respective compressibility factors
    :param Vc_arr: List or array of respective critical volumes in desired units
    :param y_arr: List or array of respective gaseous composition
    :param coeff: Parameter to return the fugacity coefficient (Default: false)
    :param R: Universal gas constant (83.14 cm^3 bar / mol K)
    :param k_arr: Array of k values for unlike solutions
    :return: List of fugacities of the species in given units of pressure
    """
    # Checking array inputs (may need for k_arr)
    if not (len(Pc_arr) == len(Tc_arr) == len(w_arr) == len(y_arr) == len(zc_arr) == len(Vc_arr)):
        raise IndexError('All arrays must be the same length')
    elif len(Pc_arr) < 2:
        raise IndexError('Solution must have at least 2 species')
    elif np.sum(y_arr) != 1:
        raise ValueError('Gas composition must add to 1')

    n = len(Pc_arr)
    B_arr = np.empty([n, n])
    del_arr = np.empty([n, n])
    phi_arr = np.empty(n)
    if k_arr is None:
        k_arr = np.zeros_like(B_arr)

    # Building 2d array of B values
    for j in range(n):
        for k in range(n):
            if j == k:
                Tcjk, Pcjk, wjk = Tc_arr[j], Pc_arr[j], w_arr[j]
                Trjk = t / Tcjk
                B_arr[j][k] = (R * Tcjk / Pcjk) * (B_0(Trjk) + (B_1(Trjk) * wjk))

            elif j < k:
                Tcjk = Tcij(Tc_arr[j], Tc_arr[k], k_arr[j][k])
                Trjk = t / Tcjk
                wjk = wij(w_arr[j], w_arr[k])
                Zcjk = Zcij(zc_arr[j], zc_arr[k])
                Vcjk = Vcij(Vc_arr[j], Vc_arr[k])
                Pcjk = R * Zcjk * Tcjk / Vcjk

                B_arr[j][k] = R * Tcjk * (B_0(Trjk) + B_1(Trjk) * wjk) / Pcjk
            else:
                B_arr[j][k] = B_arr[k][j]

    # Building 2d array of delta values
    for i in range(n):
        for j in range(n):
            if i == j:
                del_arr[i][j] = 0
            elif i < j:
                del_arr[i][j] = 2 * B_arr[i][j] - B_arr[i][i] - B_arr[j][j]
            else:
                del_arr[i][j] = del_arr[j][i]

    # Generalized phi calculations
    for i in range(n):
        sum = 0
        for j in range(n):
            for k in range(n):
                sum += (y_arr[j] * y_arr[k] * (2 * del_arr[j][i] - del_arr[j][k]))
        phi_arr[i] = np.exp(P * (B_arr[i][i] + (sum / 2)) / (R * t))

    # Outputting phi values or f values
    if coeff:
        return phi_arr
    f_arr = np.empty(n)
    for i in range(n):
        f_arr[i] = phi_arr[i] * y_arr[i] * P
    return f_arr
