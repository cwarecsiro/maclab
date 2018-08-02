from unittest import TestCase
import maclab
from maclab import colours

class TestCol(TestCase):
    def test_is_list(self):
        s = colours.col()
        self.assertTrue(isinstance(s, list))