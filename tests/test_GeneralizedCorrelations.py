import unittest
import GeneralizedCorrelations as GC

# test 1: n-Butane at 25 bar, 510K
params1 = [25, 510, 37.96, 425.1, 0.2]
# test 2: Benzene at 10 bar, 373.15K
params2 = [10, 373.15, 48.98, 562.2, 0.21]

class TestGeneralizedCorrelations(unittest.TestCase):

    def test_Pitz_Cor_2_Virial(self):
        actual = GC.pitz_Cor_2_Virial(*params1)
        self.assertEqual(round(actual, 4), .8789, "Incorrect Z value")

    def test_Pitz_Cor_2_Virial_2(self):
        actual = GC.pitz_Cor_2_Virial(*params2)
        self.assertEqual(round(actual, 4), .7223, "Incorrect Z value")

    def test_Pitz_Cor_3_Virial(self):
        actual = GC.pitz_Cor_3_Virial(*params1)
        self.assertEqual(round(actual, 4), .8760, "Incorrect Z value")

    def test_Pitz_Cor_3_Virial_2(self):
        actual = GC.pitz_Cor_3_Virial(*params2)
        self.assertEqual(round(actual, 4), .6191, "Incorrect Z value")

    # test 1: Ammonia at 310K
    # test 2: Water at 373K
    def test_rackett(self):
        actual = GC.rackett(72.47, .242, 310, 405.7)
        self.assertEqual(round(actual, 4), 28.3350, "Incorrect V value")

    def test_rackett_2(self):
        return

    def test_LGH(self):
        actual = GC.LGH(2.34,2.38,29.14)
        self.assertEqual(round(actual, 2), 28.65, "Incorrect V value")

    def test_LGH_2(self):
        return