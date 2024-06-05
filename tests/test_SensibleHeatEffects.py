import unittest
import SensibleHeatEffects as SHE
from scipy.integrate import quad

# Methane from 250C - 600C
params = [873.15, 533.15, 1.702, 9.081E-3, -2.164E-6, 0]
T, T0, A, B, C, D = params
dH = 19777.52


class TestSensibleHeatEffects(unittest.TestCase):

    def test_Cp_ig(self):
        actual = SHE.Cp_ig(T, A, B, C, D)
        self.assertEqual(round(actual, 3), 66.356, 'Incorrect Cp')

    def test_Cp_ig_2(self):
        self.assertRaises(ZeroDivisionError, SHE.Cp_ig, 0, A, B, C, D)

    def test_Cv_ig(self):
        actual = SHE.Cv_ig((7 * 8.314 / 2))
        self.assertEqual(actual, (5 * 8.314 / 2), 'Incorrect Cv')

    def test_deltaH(self):
        actual = SHE.deltaH(*params)
        self.assertEqual(round(actual), 19778, "Incorrect enthalpy change")

        dH = quad(SHE.Cp_ig, T0, T, args=(A, B, C, D))[0]
        self.assertEqual(actual, dH, "Does not match integrated Cp")

    def test_deltaH_2(self):
        params[1] = T
        self.assertEqual(SHE.deltaH(*params), 0, "Cannot handle T = T0")
        params[1] = T0

    def test_T_Given_deltaH(self):
        params[0] = dH
        actual = SHE.T_Given_deltaH(*params)
        self.assertEqual(round(actual, 2), T, 'Incorrect T given deltaH')
        params[0] = T

    def test_T_Given_deltaH_2(self):
        actual = SHE.T_Given_deltaH(-dH, T, A, B, C, D)
        self.assertEqual(round(actual, 2), T0, 'Incorrect T given -deltaH')

    def test_T_Given_deltaH_3(self):
        actual = SHE.T_Given_deltaH(0, T0, A, B, C, D)
        self.assertEqual(actual, T0, 'Incorrect T given deltaH = 0')

    def test_T_Given_deltaH_4(self):
        # R in cal/mol K
        delH = 19777.52 * 1.98 / 8.314
        actual = SHE.T_Given_deltaH(delH, T0, A, B, C, D, 1.98)
        self.assertEqual(round(actual, 2), T, 'Incorrect T given deltaH = 0')

    if __name__ == '__main__':
        unittest.main()
