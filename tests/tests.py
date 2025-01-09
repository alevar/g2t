# test_core.py
import unittest
from g2t.core import G2T

class TestG2T(unittest.TestCase):
    def test_summary(self):
        g2t = G2T('data/sample_data.csv')
        self.assertEqual(g2t.summary(), "Data summary: count=3, sum=6.0")

if __name__ == "__main__":
    unittest.main()
