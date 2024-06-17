import unittest
import ResidualProperties as RP
import numpy as np

# CO2 at 343.15 K and 150 bar
params1 = [150, 343.15, 73.83, 304.03, 0.224]
P, T, Pc, Tc, w = params1

# Propane at 450 K and 140 bar
params2 = [140, 450, 42.48, 369.8, 0.152]


class TestResidualProperties(unittest.TestCase):

    def test_res_Enthalpy(self):
        actual = RP.res_Enthalpy(*params1)
        self.assertEqual(round(actual), -4674, "Incorrect Hr")

    def test_res_Enthalpy_2(self):
        actual = RP.res_Enthalpy(*params2)
        self.assertEqual(round(actual), -7668, "Incorrect Hr")

    def test_res_Entropy(self):
        actual = RP.res_Entropy(*params1)
        self.assertEqual(round(actual,1), -9.8, "Incorrect Sr")

    def test_res_Entropy_2(self):
        actual = RP.res_Entropy(*params2)
        self.assertEqual(round(actual,1), -12.2, "Incorrect Sr")

    if __name__ == '__main__':
        unittest.main()