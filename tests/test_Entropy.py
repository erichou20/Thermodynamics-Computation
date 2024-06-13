import unittest
import Entropy as e
from scipy.integrate import quad
import numpy as np

# Methane from 550 K 5 bar to 411.34 K and 1 bar
params = [411.34, 550, 1.702, 9.081E-3, -2.164E-6, 0, 1, 5]
T, T0, A, B, C, D, P, P0 = params

# Finds delta S using numeric integration
def numeric_deltaS(T, A, B, C, D):
    return (A + B * T + C * T ** 2 + D * T ** -2) / T

class TestEntropy(unittest.TestCase):

    def test_eff_Carnot(self):
        self.assertRaises(ZeroDivisionError, e.eff_Carnot, 273, 0)

    def test_eff_Carnot_2(self):
        self.assertRaises(ValueError, e.eff_Carnot, 373,273)

    def test_eff_Carnot_3(self):
        actual = e.eff_Carnot(273,373)
        self.assertEqual(round(actual,2), 0.27, "Incorrect effeciency")

    def test_deltaS(self):
        actual = e.deltaS_ig(*params)
        self.assertEqual(round(actual), 0, "Incorrect deltaS value")
        actual = e.deltaS_ig(T0, T, A, B, C, D, P0, P)
        self.assertEqual(round(actual), 0, "Incorrect deltaS value")

    def test_deltaS_2(self):
        actual = e.deltaS_ig(T, T0, A, B, C, D, P, P)
        self.assertEqual(round(actual,3), -13.380, "Incorrect deltaS value")

        ds = quad(numeric_deltaS, T0, T, args=(A, B, C, D))[0] * 8.314
        self.assertEqual(round(ds,3), -13.380, "Incorrect deltaS value")

    def test_deltaS_3(self):
        actual = e.deltaS_ig(T, T0, A, B, C, D, 1, 10)
        self.assertEqual(round(actual,3), 5.763, "Incorrect deltaS value")

        ds = quad(numeric_deltaS, T0, T, args=(A, B, C, D))[0] - np.log(1/10)
        ds = ds * 8.314
        self.assertEqual(round(ds, 3), 5.763, "Incorrect deltaS value")

    def test_deltaS_4(self):
        self.assertRaises(ZeroDivisionError, e.deltaS_ig, 0, T, A, B, C, D, P, P0)
        self.assertRaises(ZeroDivisionError, e.deltaS_ig, T0, 0, A, B, C, D, P, P0)

    def test_T_given_deltaS(self):
        actual = e.T_given_deltaS_ig(0, T0, A, B, C, D, P, P0)
        self.assertEqual(round(actual,1), round(T,1), "Incorrect T value")

    def test_T_given_deltaS_2(self):
        actual = e.T_given_deltaS_ig(-13.380, T0, A, B, C, D)
        self.assertEqual(round(actual,2), T, "Incorrect T value")

    def test_T_given_deltaS_3(self):
        # R in cal/mol K
        dS = 5.763 * 1.98 / 8.314
        actual = e.T_given_deltaS_ig(dS, T0, A, B, C, D, 1, 10, 1.98)
        self.assertEqual(round(actual, 2), T, "Incorrect T value")

    #def test_P_given_deltaS(self):
    #    actual = e.P_given_deltaS_ig(0, T, T0, A, B, C, D, P0)
    #    self.assertEqual(round(actual,2), P, "Incorrect P value")


    if __name__ == '__main__':
        unittest.main()



