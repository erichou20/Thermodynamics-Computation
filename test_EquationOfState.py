import unittest
import EquationOfState as EOS


class TestEquationOfState(unittest.TestCase):

    # add Volume func given P,R,T,Z
    #explain test description
    def test_PR_gas(self):
        actual = EOS.PR_gas(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 3), .808)

    def test_PR_gas_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_PR_liquid(self):
        actual = EOS.PR_liquid(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .0366)

    def test_PR_liquid_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_RK_liquid(self):
        actual = EOS.RK_liquid(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .0433)

    def test_RK_liquid_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_RK_gas(self):
        actual = EOS.RK_gas(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .8305)

    def test_RK_gas_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_SRK_liquid(self):
        actual = EOS.SRK_liquid(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .0415)

    def test_SRK_liquid_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_SRK_gas(self):
        actual = EOS.SRK_gas(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .8191)

    def test_SRK_gas_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_vdW_liquid(self):
        actual = EOS.vdW_liquid(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .0621)

    def test_vdW_liquid_2(self):
        # Benzene, 10 bar, 373K
        return

    def test_vdW_gas(self):
        actual = EOS.vdW_gas(9.4573, 350, 37.96, 425.1, 0.2)
        self.assertEqual(round(actual, 4), .8667)

    def test_vdW_gas_2(self):
        # Benzene, 10 bar, 373K
        return


    def test_virial_gas(self):
        return


