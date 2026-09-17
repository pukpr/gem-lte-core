import math
import unittest

import numpy as np

from fit_column3_by_column4 import plateau_frequency


class PlateauFrequencyTests(unittest.TestCase):
    def test_uses_monthly_plateau_increment_instead_of_reset(self):
        time = 2000.0 + np.arange(8) / 12.0
        x = np.array([0.0, 0.25, 0.5, 0.75, 5.0, 5.25, 5.5, 5.75])

        slope, frequency = plateau_frequency(time[::-1], x[::-1])

        self.assertAlmostEqual(slope, 0.25)
        self.assertAlmostEqual(frequency, 2.0 * math.pi / 0.25)


if __name__ == "__main__":
    unittest.main()
