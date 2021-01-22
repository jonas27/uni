import unittest
import gotoh
from gotoh_helpers import read_fasta_file

class Test_GotohUnitTest(unittest.TestCase):

    def test_gotoh1(self):
        tests = [
            # (helpers.read_substitution_matrix('data/pam250.txt')
            #  Edge case testing
            ([' ','A','T'], [' ','T','A'],  -3,-1,None,
            (1,-2,
            'AT',
            'TA',
            '||',)),
            ([' ','A','T','T','A'], [' ','T','A'],  -3,-1,None,
            (1,-3,
            'ATTA',
            '__TA',
            '  **',)),
            ([' ','T','A'],[' ','A','T','T','A'],  -3,-1,None,
            (1,-3,
            '__TA',
            'ATTA',
            '  **',)),
            # Real tests
            (read_fasta_file('data/test1_1.fasta'), read_fasta_file('data/test1_2.fasta'),
            -3,-1,None,
            (10,-10,
            'AATCCGCTAGAAACCCTTAG',
            '_ATC______TTACCTTCAT',
            ' ***      ||***|*|*|',)),
            (read_fasta_file('data/test2_1.fasta'), read_fasta_file('data/test2_2.fasta'),
            -3,-1,None,
            (1,2,
            'ATTCACTAGATTACCTTCACT',
            '________GATTACCTTCACT',
            '        *************',)),
        ]
        for t in tests:
            score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3])
            should = t[5]
            self.assertEqual(should[0],len(alignments))
            self.assertEqual(should[1],score)
            self.assertEqual(should[2],''.join(alignments[0][0]))
            self.assertEqual(should[3],''.join(alignments[0][1]))
            self.assertEqual(should[4],''.join(alignments[0][2]))

if __name__ == "__main__":
    unittest.main()

