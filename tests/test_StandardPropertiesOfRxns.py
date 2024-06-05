import unittest
import StandardPropertiesOfRxns as SPR

class TestStandardPropertiesOfRxns(unittest.TestCase):

    # test 1: Heat of combustion for 4HCl + O2 -> 2H2O + 2Cl2
    def test_standard_Prop(self):
        vals = [-92307, 0, -241818,0]
        vs = [-4, -1, 2, 2]
        self.assertEqual(SPR.standard_Prop(vals, vs), -114408, "Incorrect dH")

    def test_standard_Prop2(self):
        vals = [-92307, 0, -241818,0]
        vs = [-4, -1, 2]
        self.assertRaises(IndexError, SPR.standard_Prop, vals, vs)

    if __name__ == '__main__':
        unittest.main()
