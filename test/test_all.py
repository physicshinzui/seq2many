import unittest 
from mutantsGenetator import *

class TestAll(unittest.TestCase):

    def test_ref2mut(self):
        seq = "ALKKL"
        ans = ref2mut(seq, 0, "K")
        expected = "KLKKL"
        self.assertEqual(ans, expected)

    def test_hamming_distance(self):
        seq1, seq2 = "AAL", "ALK"
        hd = hamming_distance(seq1, seq2)
        expected = 2
        self.assertEqual(hd, expected)

    def test_multiple_mutation1(self):
        seq = "AAAAA"
        mutations_dict = {'K':2, 'L':4}
        ans = multiple_mutation(seq, mutations_dict)
        expected = "AAKAL"
        self.assertEqual(ans, expected)
        print(ans, expected)

    def test_multiple_mutation2(self):
        seq = "AAAAA"
        mutations_dict = {'K':2, 'L':5}
        ans = multiple_mutation(seq, mutations_dict)
        expected = "AAKAA"
        self.assertNotEqual(ans, expected)
        print(ans, expected)

    def test_deep_mutational_scanning1(self):
        seq = "AAA"
        seqs = deep_mutational_scanning(seq, 0, 3)
        n_seqs = len(seqs)
        expected = 20*3
        self.assertEqual(n_seqs, expected)
        print(n_seqs, expected)

    def test_deep_mutational_scanning2(self):
        seq = "AAA"
        seqs = deep_mutational_scanning(seq, 0, 5)
        n_seqs = len(seqs)
        expected = 20*5
        self.assertEqual(n_seqs, expected)
        print(n_seqs, expected)


