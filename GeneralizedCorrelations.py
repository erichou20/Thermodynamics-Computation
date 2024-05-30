import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import root
from scipy.integrate import quad


def pitz_Cor_2_Virial(P, T, Pc, Tc, w, flag=False):
    Tr = T / Tc
    Pr = P / Pc
    B_hat = (0.083 - (.422 / (Tr ** 1.6))) + w * (0.139 - (.172 / (Tr ** 4.2)))

    if flag:
        return B_hat
    return 1 + B_hat * Pr / Tr


def pitz_Cor_3_Virial(P, T, Pc, Tc, w):
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
    Tr = T / Tc
    return Vc * Zc ** (1 - Tr) ** (2 / 7)


def LGH(rho1, rho2, V1):
    return V1 * rho1 / rho2
