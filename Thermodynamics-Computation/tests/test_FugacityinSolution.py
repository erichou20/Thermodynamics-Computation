import numpy as np
import unittest
from src.Fugacity import FugacityinSolution as FS

# test 1: System of methane(1)/ethane(2)/propane(3)
P = 35  # Bar
t = 373  # K
params = [[45.99, 48.72, 42.48],  # Crit. pressures in bar
          [190.6, 305.3, 369.8],  # Crit. temps in K
          [0.012, 0.1, 0.152],  # Acentric factors
          [0.286, 0.279, 0.276],  # Crit. compressibility factors
          [98.6, 145.5, 200.0]]  # Crit. volumes in cm^3/mol
y_arr = [0.21, 0.43, 0.36]

# test 2: mixture of ethylene(1)/propylene(2)
P2 = 30  # Bar
t2 = 423.15  # K
params2 = [[50.4, 46.65],  # Crit. pressures in bar
           [282.3, 365.6],  # Crit. temps in K
           [0.087, 0.140],  # Acentric factors
           [0.281, 0.289],  # Crit. compressibility factors
           [131, 188.4],  # Crit. volumes in cm^3/mol
           [0.35, 0.65]]  # gaseous compositions


class TestFugacityinSolution(unittest.TestCase):
    def test_fugacity_In_Solution_1(self):
        self.assertRaises(IndexError, FS.fugacity_Solution, P, t, *params, y_arr=[0, 1])
        self.assertRaises(IndexError, FS.fugacity_Solution, P, t, *params, y_arr=[0, 1, 1, 1])
        self.assertRaises(IndexError, FS.fugacity_Solution, P, t, [1], [1], [1], [1], [1], [1])
        self.assertRaises(ValueError, FS.fugacity_Solution, P, t, *params, y_arr=[0.33, 0.33, 0.33])
        self.assertRaises(ValueError, FS.fugacity_Solution, P, t, *params, y_arr=[0.33, 0.33, 0.35])

    def test_fugacity_In_Solution_2(self):
        actual = np.round(FS.fugacity_Solution(P, t, *params, y_arr), 2)
        self.assertEqual(actual[0], 7.49, "Incorrect fugacity of species 1")
        self.assertEqual(actual[1], 13.25,  "Incorrect fugacity of species 2")
        self.assertEqual(actual[2], 9.76, "Incorrect fugacity of species 3")
        actual2 = np.round(FS.fugacity_Solution(P, t, *params, y_arr, True), 3)
        self.assertEqual(actual2[0], 1.019, "Incorrect fugacity coefficient of species 1")
        self.assertEqual(actual2[1], 0.881, "Incorrect fugacity coefficient of species 2")
        self.assertEqual(actual2[2], 0.775, "Incorrect fugacity coefficient of species 3")

    def test_fugacity_In_Solution_3(self):
        actual = np.round(FS.fugacity_Solution(P2, t2, *params2), 2)
        self.assertEqual(actual[0], 10.05, "incorrect fugacity of species 1")
        self.assertEqual(actual[1], 17.06, "incorrect fugacity of species 2")
        actual2 = np.round(FS.fugacity_Solution(P2, t2, *params2, coeff=True), 3)
        self.assertEqual(actual2[0], 0.957, "incorrect fugacity coefficient of species 1")
        self.assertEqual(actual2[1], 0.875, "incorrect fugacity coefficient of species 2")

    if __name__ == '__main__':
        unittest.main()
