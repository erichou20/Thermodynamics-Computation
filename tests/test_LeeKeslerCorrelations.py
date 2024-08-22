import unittest
from src.Misc import LeeKeslerCorrelations as LKC


class TestLeeKeslerCorrelations(unittest.TestCase):

    # test 1: Lee-Kesler correlations for Z given Z0, Z1
    # test 2: Lee-Kesler correlations for Bhat given B0, B1

    def test_lee_Kesler(self):
        actual = LKC.lee_Kesler(.865,.038,.2)
        self.assertEqual(round(actual,3), .873)

    def test_lee_Kesler_2(self):
        actual = LKC.lee_Kesler(-0.232,.059,.2)
        self.assertEqual(round(actual,3), -.220)

    if __name__ == '__main__':
        unittest.main()