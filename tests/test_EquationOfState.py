import unittest
import EquationOfState as EOS
# test 1: n-Butane at 9.4573 bar, 350K
params1 = [9.4573, 350, 37.96, 425.1, 0.2]
# test 2: Benzene at 10 bar, 373.15K
params2 = [10, 373.15, 48.98, 562.2, 0.21]

# Isopropanol at 200C, 10 bar
T = 473.15
B = -388
C = -26000
P = 10
Vig = 83.14*T/P

class TestEquationOfState(unittest.TestCase):
    def test_PR_liquid(self):
        actual = EOS.PR_liquid(*params1)
        self.assertEqual(round(actual, 4), .0366, "Incorrect Z value")

    def test_PR_liquid_2(self):
        actual = EOS.PR_liquid(*params2)
        self.assertEqual(round(actual, 4), .0307, "Incorrect Z value")

    def test_PR_gas(self):
        actual = EOS.PR_gas(*params1)
        self.assertEqual(round(actual, 4), .8081, "Incorrect Z value")

    def test_PR_gas_2(self):
        actual = EOS.PR_gas(*params2)
        self.assertEqual(round(actual, 4), .6537, "Incorrect Z value")

    def test_RK_liquid(self):
        actual = EOS.RK_liquid(*params1)
        self.assertEqual(round(actual, 4), .0433, "Incorrect Z value")

    def test_RK_liquid_2(self):
        actual = EOS.RK_liquid(*params2)
        self.assertEqual(round(actual, 4), .0357, "Incorrect Z value")

    def test_RK_gas(self):
        actual = EOS.RK_gas(*params1)
        self.assertEqual(round(actual, 4), .8305, "Incorrect Z value")

    def test_RK_gas_2(self):
        actual = EOS.RK_gas(*params2)
        self.assertEqual(round(actual, 4), .7082, "Incorrect Z value")

    def test_SRK_liquid(self):
        actual = EOS.SRK_liquid(*params1)
        self.assertEqual(round(actual, 4), .0415, "Incorrect Z value")

    def test_SRK_liquid_2(self):
        actual = EOS.SRK_liquid(*params2)
        self.assertEqual(round(actual, 4), .0347, "Incorrect Z value")

    def test_SRK_gas(self):
        actual = EOS.SRK_gas(*params1)
        self.assertEqual(round(actual, 4), .8191, "Incorrect Z value")

    def test_SRK_gas_2(self):
        actual = EOS.SRK_gas(*params2)
        self.assertEqual(round(actual, 4), .6626, "Incorrect Z value")

    def test_vdW_liquid(self):
        actual = EOS.vdW_liquid(*params1)
        self.assertEqual(round(actual, 4), .0621, "Incorrect Z value")

    def test_vdW_liquid_2(self):
        actual = EOS.vdW_liquid(*params2)
        self.assertEqual(round(actual, 4), .0522, "Incorrect Z value")

    def test_vdW_gas(self):
        actual = EOS.vdW_gas(*params1)
        self.assertEqual(round(actual, 4), .8667, "Incorrect Z value")

    def test_vdW_gas_2(self):
        actual = EOS.vdW_gas(*params2)
        self.assertEqual(round(actual, 4), .8080, "Incorrect Z value")

    def test_virial_Z(self):
        actual = EOS.virial_Z(Vig, B)
        self.assertEqual(round(actual, 4), .9014, "Incorrect Z value")

    def test_virial_Z_2(self):
        actual = EOS.virial_Z(Vig, B, C)
        self.assertEqual(round(actual, 4), .8997,"Incorrect Z value")

    def test_virial_V(self):
        actual = EOS.virial_V(P, T, B)
        self.assertEqual(round(actual), 3497, "Incorrect V value")

    def test_virial_V_2(self):
        actual = EOS.virial_V(P, T, B, C)
        self.assertEqual(round(actual), 3488, "Incorrect V value")

    if __name__ == '__main__':
        unittest.main()