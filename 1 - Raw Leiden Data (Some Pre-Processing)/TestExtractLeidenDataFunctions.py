import unittest  # Testing framework
from ExtractLeidenDataFunctions import *  # Module under test

class TestExtractLeidenData(unittest.TestCase):
    def test_match_column_order(self):
        # Simple example with two swaps
        header_expected = ['a', 'b', 'c', 'd']
        header_input = ['b', 'a', 'd', 'c']
        data = [[2, 1, 4, 3], [2, 1, 4, 3], [2, 1, 4, 3]]
        result = [[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]
        self.assertEquals(match_column_order(header_expected, header_input, data), result)

        # Headers do not contain the same elements
        header_expected = ['a', 'b', 'c', 'd']
        header_input = ['b', 'a', 'd', 'f']
        data = [[2, 1, 4, 3], [2, 1, 4, 3], [2, 1, 4, 3]]
        self.assertRaises(ValueError, match_column_order, header_expected, header_input, data)

        # Headers are not the same length
        header_expected = ['a', 'b', 'c', 'd', 'e']
        header_input = ['b', 'a', 'd', 'f']
        data = [[2, 1, 4, 3], [2, 1, 4, 3], [2, 1, 4, 3]]
        self.assertRaises(ValueError, match_column_order, header_expected, header_input, data)

if __name__ == '__main__':
    unittest.main()

