"""

"""
import unittest
from . import msa as msa
from . import helpers as h

class Test_MSNUnitTest(unittest.TestCase):
    """
    TestGotohAssignment1 mainly tests A2 but also some edge cases
    This is a convience wrapper to run the assignment as a test.
    """
    def test_msn(self):
        tests = [
            [
                [list(" ACCAT"), list(" ACGGAT"), list(" AACCAT")],
                [0,-1],
                [
                    ['#', 'A', 'C', '#', 'C', 'A', 'T'],
                    ['#', 'A', 'C', 'G', 'G', 'A', 'T'],
                    ['A', 'A', 'C', '#', 'C', 'A', 'T'],
                ]
            ],
            [
                [list(" ACCAT"), list(" ACGGAT"), list(" AACCAT")],
                [-1,-1],
                [
                    ['#', 'A', 'C', '#', 'C', 'A', 'T'],
                    ['#', 'A', 'C', 'G', 'G', 'A', 'T'],
                    ['A', 'A', 'C', '#', 'C', 'A', 'T'],
                ]
            ],
            # [
            #     [list(" A"), list(" ACGGAT"), list(" AACCATAACCAT")],
            #     [-1,-1],
            #     [
            #         ['#', 'A', '#', '#', '#', '#', '#', '#', '#', '#', '#', '#'],
            #         ['#', 'A', 'C', '#', '#', '#', '#', '#', 'G', 'G', 'A', 'T'],
            #         ['A', 'A', 'C', 'C', 'A', 'T', 'A', 'A', 'C', 'C', 'A', 'T'],
            #     ]
            # ],
            [
                [list(" AAAAATTTACATTTAA"), list(" ATTTTACATTTGGG"), list(" ATTTTACATTTGGGG")],
                [-1,-1],
                [
                    ['A', 'A', 'A', 'A', 'A', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', '#', 'A', 'A'], ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', 'G', 'G', 'G'], ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', 'G', 'G', 'G', 'G'],
                ]
            ],
        ]
        for t in tests:
            aligns = msa.run(t[0][0],t[0][1],t[0][2],t[1][0],t[1][1])
            self.assertEqual(aligns, t[2])

if __name__ == "__main__":
    unittest.main()


        # s1 = h.read_fasta_file("data/s1.fasta")
        # s2= h.read_fasta_file("data/s2.fasta")
        # s3 = h.read_fasta_file("data/s3.fasta")