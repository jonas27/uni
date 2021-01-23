import itertools
import unittest

import gotoh
from gotoh_helpers import read_fasta_file, read_substitution_matrix

class Test_GotohUnitTest(unittest.TestCase):

    def test_gotoh1(self):
        tests = []
            # Real tests
        for combination in list(itertools.combinations(range(1,5), 2)):
            tests.append([read_fasta_file('data/s'+str(combination[0])+'.fasta'), 
            read_fasta_file('data/s'+str(combination[1])+'.fasta'),
            -11,-1, 
            read_substitution_matrix('data/pam250.txt'),
            'Tested sequences: ' + str(combination)])
        for i in range(2):
            f = None
            if i == 0:
                open("scores_blossum.csv", "w").close()
                f = open("scores_blossum.csv", "a")
            elif i == 1:
                open("scores_pam.csv", "w").close()
                f = open("scores_pam.csv", "a")
            for t in tests:
                if i == 0:
                    t[4] = read_substitution_matrix('data/pam250.txt') 
                elif i == 1:
                    t[4] = read_substitution_matrix('data/blosum62.txt')
                score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3],t[4])
                f.write(t[5]+'\nWith score: ' + str(score) + ' and ' + str(len(alignments)) + ' possible alignments. \n')
                for a in alignments:
                    f.write(''.join(a[0])+'\n')
                    f.write(''.join(a[2])+'\n')
                    f.write(''.join(a[1])+'\n\n')
                f.write('\n\n\n') 
            f.close()


if __name__ == "__main__":
    unittest.main()

