import unittest

import LatentHeats as LH

#test 1: Water
class TestLatentHeats(unittest.TestCase):

    def test_troutons(self):
        self.assertEqual(round(LH.troutons(373)), 31011, "Incorrect trouton estimation")

    def test_trounts_2(self):
        self.assertEqual(round(LH.troutons(373, 83.14)), 310112,
                         "Incorrect dH estimation")

    def test_Riedel(self):
        actual = LH.riedel(647.1, 220.55, 373.15)
        self.assertEqual(round(actual), 42024, "Incorrect riedel estimation")

    def test_Riedel_2(self):
        return

    # Water estimated at 300C with known dH at 100C
    def test_watson(self):
        actual = LH.watson(2257, 373.15, 573.15, 647.1)
        self.assertEqual(round(actual), 1372, "Incorrect watson estimation")

    def test_watson_2(self):
        actual = LH.watson(2257, 373.15, 373.15, 647.1)
        self.assertEqual(actual,2257,"Incorrect watson estimation")

    def test_watson_3(self):
        self.assertRaises(TypeError, LH.watson, 373.15, 673.15, 373.15)

    if __name__ == '__main__':
        unittest.main()

#dH = watson(2257,373.15,573.15,647.1)
#print(dH)
