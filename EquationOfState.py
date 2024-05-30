import numpy as np
from scipy.optimize import root


def gen_Cubic_EOS_V(Z, b, q, o, e):
    return 1 + b - q * b * (Z - b) / ((Z + e * b) * (Z + o * b)) - Z


def gen_Cubic_EOS_L(Z, b, q, o, e):
    return b + (Z + e * b) * (Z + o * b) * (1 + b - Z) / (q * b) - Z


params = {
    "PR": [(1 + np.sqrt(2)), (1 - np.sqrt(2)), 0.07780, 0.45724, 0.30740],
    "RK": [1, 0, 0.08664, 0.42748, (1 / 3)],
    "SRK": [1, 0, 0.08664, 0.42748, (1 / 3)],
    "vdW": [0, 0, (1 / 8), (27 / 64), (3 / 8)],
}


def PR_gas(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["PR"][2] * Pr / Tr
    q = params["PR"][3] * a / (params["PR"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["PR"][0], params["PR"][1]))
    return sol.x[0]


def PR_liquid(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.37464 + 1.54226 * w - 0.26992 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["PR"][2] * Pr / Tr
    q = params["PR"][3] * a / (params["PR"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["PR"][0], params["PR"][1]))
    return sol.x[0]


def RK_gas(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = Tr ** (-0.5)
    b = params["RK"][2] * Pr / Tr
    q = params["RK"][3] * a / (params["RK"][2] * Tr)
    print(q, b)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["RK"][0], params["RK"][1]))
    return sol.x[0]


def RK_liquid(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = Tr ** (-0.5)
    b = params["RK"][2] * Pr / Tr
    q = params["RK"][3] * a / (params["RK"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["RK"][0], params["RK"][1]))
    return sol.x[0]


def SRK_gas(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.480 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["SRK"][2] * Pr / Tr
    q = params["SRK"][3] * a / (params["SRK"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["SRK"][0], params["SRK"][1]))
    return sol.x[0]


def SRK_liquid(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = (1 + (0.480 + 1.574 * w - 0.176 * w ** 2) * (1 - np.sqrt(Tr))) ** 2
    b = params["SRK"][2] * Pr / Tr
    q = params["SRK"][3] * a / (params["SRK"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["SRK"][0], params["SRK"][1]))
    return sol.x[0]


def vdW_gas(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = 1
    b = params["vdW"][2] * Pr / Tr
    q = params["vdW"][3] * a / (params["vdW"][2] * Tr)

    sol = root(gen_Cubic_EOS_V, [1], args=(b, q, params["vdW"][0], params["vdW"][1]))
    return sol.x[0]


def vdW_liquid(P, T, Pc, Tc, w):
    Tr = T / Tc
    Pr = P / Pc
    a = 1
    b = params["vdW"][2] * Pr / Tr
    q = params["vdW"][3] * a / (params["vdW"][2] * Tr)

    sol = root(gen_Cubic_EOS_L, [0.01], args=(b, q, params["vdW"][0], params["vdW"][1]))
    return sol.x[0]


def virial_Z(P, V, T, B, C=0, R=83.14):
    Z = 1 + (B / V) + (C / V ** 2)
    return Z


def virial_V(P, T, B, C=0, R=83.14):
    def solver(V):
        return (R * T * (1 + (B / V) + (C / V ** 2)) / P) - V

    V0 = R * T / P
    sol = root(solver, [V0])
    return sol.x[0]
