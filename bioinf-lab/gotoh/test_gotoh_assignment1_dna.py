import unittest
import gotoh
from gotoh_helpers import read_fasta_file

class Test_GotohUnitTest(unittest.TestCase):

    def test_gotoh_edgecases(self):
        tests = [
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
        ]
        for t in tests:
            score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3])
            should = t[5]
            self.assertEqual(should[0],len(alignments))
            self.assertEqual(should[1],score)
            self.assertEqual(should[2],''.join(alignments[0][0]))
            self.assertEqual(should[3],''.join(alignments[0][1]))
            self.assertEqual(should[4],''.join(alignments[0][2]))
        g = gotoh.Gotoh()
        for t in tests:
            score, alignments = g.run(t[0],t[1],t[2],t[3])
            should = t[5]
            self.assertEqual(should[0],len(alignments))
            self.assertEqual(should[1],score)
            self.assertEqual(should[2],''.join(alignments[0][0]))
            self.assertEqual(should[3],''.join(alignments[0][1]))
            self.assertEqual(should[4],''.join(alignments[0][2]))

    def test_goto_assignment(self):
        tests =[
            (read_fasta_file('data/test1_1.fasta'), read_fasta_file('data/test1_2.fasta'),
            -3,-1,None,
            (10,-10,
            'AATCCGCTAGAAACCCTTAG',
            '_ATC______TTACCTTCAT',
            ' ***      ||***|*|*|',),
            'results/assignment1_1.txt'),
            (read_fasta_file('data/test2_1.fasta'), read_fasta_file('data/test2_2.fasta'),
            -3,-1,None,
            (1,2,
            'ATTCACTAGATTACCTTCACT',
            '________GATTACCTTCACT',
            '        *************',),
            'results/assignment1_2.txt',),
        ]
        g = gotoh.Gotoh()
        gotoh_funcs = [gotoh.gotoh, g.run]
        for func in gotoh_funcs:
            for t in tests:
                score, alignments = func(t[0],t[1],t[2],t[3])
                should = t[5]
                self.assertEqual(should[0],len(alignments))
                self.assertEqual(should[1],score)
                self.assertEqual(should[2],''.join(alignments[0][0]))
                self.assertEqual(should[3],''.join(alignments[0][1]))
                self.assertEqual(should[4],''.join(alignments[0][2]))
        f = None
        for t in tests:
            score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3])
            # open(t[6], "w").close()
            f = open(t[6], "w")
            score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3],t[4])
            f.write('Score ' + str(score) + ' and ' + str(len(alignments)) + ' optimal alignments. \n')
            for a in alignments:
                f.write(''.join(a[0])+'\n')
                f.write(''.join(a[2])+'\n')
                f.write(''.join(a[1])+'\n\n')
            f.write('\n\n\n') 
            f.close()
        
if __name__ == "__main__":
    unittest.main()

