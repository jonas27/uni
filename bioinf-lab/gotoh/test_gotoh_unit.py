import unittest
from unittest.mock import patch
import gotoh
import gotoh_helpers as helpers


class Test_GotohUnitTest(unittest.TestCase):
    
    def test_complete_d_p_q_computation(self):
        seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
        d,p,q = gotoh.complete_d_p_q_computation(seq1,seq2,-3,-1,gotoh.dna_sub)
        self.assertEqual(d[1], [-4, 1, -3, -4, -5]) 
        self.assertEqual(p[2], [-10000000000, -3, -7, -8, -9])
        self.assertEqual(q[3], [-1000000, -10, -8, -8, -3])

    def test_compute_tracebacks(self):
        seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
        d,p,q = gotoh.complete_d_p_q_computation(seq1,seq2,-3,-1,gotoh.dna_sub)
        p = gotoh.compute_tracebacks(seq1, seq2, d,p,q,-3,-1,gotoh.dna_sub)
        self.assertIn([((0, 0), 'd', 4), ((0, 1), 'd', 3), ((1, 2), 'd', 2), ((2, 3), 'd', 1), ((3, 4), 'd', 0)],p)

    def test_initd(self):
        tests = [
            ([" ","A","B"],[" ","A","B","C","D"],[-5, None, None, None, None]),
            ([" ","A","B"],[" ","A","B"],[-4, None, None])
        ]
        for test in tests:
            result = gotoh.initd(test[0],test[1],-3,-1)
            self.assertIn(test[2], result)

    def test_initp(self):
        tests = [
            ([" ","A","B"],[" ","A","B","C","D"],[0, -1000000, -1000000, -1000000, -1000000]),
            ([" ","A","B"],[" ","A","B"],[-10000000000, -10000000000,-10000000000])
        ]
        for test in tests:
            result = gotoh.initp(test[0],test[1],)
            self.assertIn(test[2], result)

    def test_initq(self):
        tests = [
            ([" ","A","B"],[" ","A","B","C","D"],[0, -10000000000, -10000000000, -10000000000, -10000000000]),
            ([" ","A","B"],[" ","A","B"],[-1000000, -10000000000, -10000000000])
        ]
        for test in tests:
            result = gotoh.initq(test[0],test[1])
            self.assertIn(test[2], result)

    def test_find_previous(self):
        seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
        d,p,q = gotoh.complete_d_p_q_computation(seq1,seq2,-3,-1,gotoh.dna_sub)
        il, jl = gotoh.size(d)
        p = gotoh.find_previous(((il-1,jl-1), "d", 0), seq1, seq2, d,p,q,-3,-1,[],gotoh.dna_sub)
        self.assertEqual(p[0][1],'d')

    def test_read_fasta_file(self):
        tests = [
            (helpers.read_fasta_file("data/s1.fasta"), ' ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN')
        ]
        for t in tests:
            self.assertEqual(''.join(t[0]),t[1])
        
    def test_read_substitution_matrix(self):
        scores = helpers.read_substitution_matrix("data/pam250.txt")
        self.assertEqual(scores.get(('*','*')),1)
        self.assertEqual(scores.get(('R','A')),-2)

    def test_size(self):
        tests = [
            ([" ","A","B"],[" ","A","B","C","D"],(3,5)),
            ([" ","A","B"],[" ","A","B","C","D","E","A","B","C"],(3,9))
        ]
        for test in tests:
            p = gotoh.initp(test[0],test[1])
            self.assertEqual(gotoh.size(p),test[2])

    @patch('builtins.print')
    def test_show(self, mock_print):
        tests = [
            ([" ","A","B"],[" ","A","B","C","D"],[-10000000000, -10000000000, -10000000000, -10000000000, -10000000000]),
            ([" ","A","B"],[" ","A","B"],[-10000000000, -10000000000, -10000000000])
        ]
        for test in tests:
            p = gotoh.initp(test[0],test[1])
            helpers.show(p)
            mock_print.assert_called_with(test[2])



if __name__ == "__main__":
    unittest.main()

