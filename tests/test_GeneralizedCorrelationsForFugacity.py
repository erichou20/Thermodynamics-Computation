from src.Fugacity import GeneralizedCorrelationsForFugacity as GCF
import numpy as np
from scipy.stats import gmean
import unittest

# test 1 for virial coefficient data: mixture of ethylene(1)/propylene(2) as gases
params = np.array([[282.3, 365.6],
                   [50.4, 46.64],
                   [131.0, 188.4],
                   [0.281, 0.289],
                   [0.087, 0.140]])


class TestGeneralizedCorrelationsForFugacity(unittest.TestCase):

    def test_pitzer_Fugacity_1(self):
        actual = GCF.pitzer_Fugacity(5.28, 313.15, 36.48, 408.1, 0.181)
        self.assertEqual(round(actual, 2), 4.69, "Incorrrect Fugactiy")

    def test_pitzer_Fugacity_2(self):
        actual = GCF.pitzer_Fugacity(5.28, 313.15, 36.48, 408.1, 0.181, True)
        self.assertEqual(round(actual, 3), 0.888, "Incorrrect Fugactiy")

    def test_pitzer_Fugacity_3(self):
        return

    def test_pitzer_Fugacity_4(self):
        return

    def test_wij_1(self):
        actual = GCF.wij(params[4][0], params[4][1])
        self.assertEqual(round(actual, 3), 0.114, "Incorrrect omega coefficient")

    def test_wij_2(self):
        for i in range(10):
            for j in range(10):
                actual = GCF.wij(i, j)
                self.assertEqual(actual, np.average([i, j]), "Incorrect omega calculation")

    def test_Tcij_1(self):
        actual = GCF.Tcij(params[0][0], params[0][1], 0)
        self.assertEqual(round(actual, 1), 321.3, "Incorrrect Tcij coefficient")

    def test_Tcij_2(self):
        for i in range(10):
            for j in range(10):
                actual = GCF.Tcij(i, j)
                self.assertEqual(round(actual,10), round(gmean([i, j]), 10),
                                 "Incorrrect Tcij calculation")

    def test_Zcij_1(self):
        actual = GCF.Zcij(params[3][0], params[3][1])
        self.assertEqual(round(actual, 3), 0.285, "Incorrrect Zcij coefficient")

    def test_Zcij_2(self):
        for i in range(10):
            for j in range(10):
                actual = GCF.wij(i, j)
                self.assertEqual(actual, np.average([i, j]), "Incorrect Zcij calculation")
                self.assertEqual(actual, GCF.wij(i,j), "Incorrect Zcij calculation")
    def test_Vcij_1(self):
        actual = GCF.Vcij(params[2][0], params[2][1])
        self.assertEqual(round(actual, 2), 157.97, "Incorrrect Vcij coefficient")

    def test_Vcij_2(self):
        for i in range(10):
            for j in range(10):
                actual = GCF.Vcij(i, j)**(1/3)
                self.assertEqual(actual, (i**(1/3) + j**(1/3))/2, "Incorrrect Vcij calculation")

    if __name__ == '__main__':
        unittest.main()
