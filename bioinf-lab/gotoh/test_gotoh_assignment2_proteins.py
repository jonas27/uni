"""
test_gotoh_assignment2_proteins fullfils the requirments for A2.
It outputs the files './results/scores_pam250.csv' and './results/scores_blosum62.csv'.
It then checks the output against the expected file.
Expected is here a output I considered correct.
"""

import itertools
import unittest

import gotoh
from gotoh_helpers import read_fasta_file, read_substitution_matrix

class TestGotohAssignment2(unittest.TestCase):
    """
    TestGotohAssignment2 is a convience class to make A2 as test
    """

    def test_gotoh_assognment2(self):
        """
        test_gotoh_assognment2 manipulates data in results folder.
        If results is missing please create, if not tests will fail.
        """
        tests = []
        for combination in list(itertools.combinations(range(1,5), 2)):
            tests.append([read_fasta_file('data/s'+str(combination[0])+'.fasta'),
            read_fasta_file('data/s'+str(combination[1])+'.fasta'),
            -11,-1,
            read_substitution_matrix('data/pam250.txt'),
            'Tested sequences: ' + str(combination)])
        for i in range(2):
            f = None
            if i == 0:
                f = open("results/scores_blosum62.csv", "w")
            elif i == 1:
                f = open("results/scores_pam250.csv", "w")
            for t in tests:
                if i == 0:
                    t[4] = read_substitution_matrix('data/blosum62.txt')
                elif i == 1:
                    t[4] = read_substitution_matrix('data/pam250.txt')
                score, alignments = gotoh.gotoh(t[0],t[1],t[2],t[3],t[4])
                f.write(t[5]+'\nWith score ' + str(score) + ' and '
                    + str(len(alignments)) + ' optimal alignments. \n')
                for a in alignments:
                    f.write(''.join(a[0])+'\n')
                    f.write(''.join(a[2])+'\n')
                    f.write(''.join(a[1])+'\n\n')
                f.write('\n\n\n')
            f.close()
        iss = ["results/scores_pam250.csv","results/scores_blosum62.csv"]
        shoulds = ["results/scores_pam250_should.csv","results/scores_blosum62_should.csv"]
        for count, should in enumerate(shoulds):
            file_is = open(iss[count], "r")
            lines_is = file_is.readlines()
            file_should = open(should, "r")
            lines_should = file_should.readlines()
            for i in range(len(lines_is)):
                self.assertEqual(lines_is[count], lines_should[count])
            file_is.close()
            file_should.close()

if __name__ == "__main__":
    unittest.main()
