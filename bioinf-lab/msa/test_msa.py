"""
Tests for the msa implementation
"""
import unittest
from msa import run

class Test_MSNTest(unittest.TestCase):
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
            [
                [list(" AAAAATTTACATTTAA"), list(" ATTTTACATTTGGG"), list(" ATTTTACATTTGGGG")],
                [-1,-1],
                [
                    ['A', 'A', 'A', 'A', 'A', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', '#', 'A', 'A'],
                    ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', 'G', 'G', 'G'],
                    ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', 'G', 'G', 'G', 'G'],
                ]
            ],
            [
                [list(" AAAAATTTACATTTAA"), list(" ATTTTACATTTGGG"), list(" ATTTTACATTTGGGG")],
                [-3,-1],
                [
                    ['A', 'A', 'A', 'A', 'A', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', '#', 'A', 'A'],
                    ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', '#', 'G', 'G', 'G'],
                    ['#', '#', '#', 'A', 'T', 'T', 'T', 'T', 'A', 'C', 'A', 'T', 'T', 'T', 'G', 'G', 'G', 'G'],
                ]
            ],
        ]
        for t in tests:
            aligns = run(t[0][0],t[0][1],t[0][2],t[1][0],t[1][1])
            self.assertEqual(aligns, t[2])
        t = tests.pop()
        aligns = run(t[0][0],t[0][1],t[0][2],t[1][0],t[1][1])
        f = open("result.txt", "w")
        f.write('Alignment:\n')
        f.write(''.join(aligns[0])+'\n')
        f.write(''.join(aligns[1])+'\n')
        f.write(''.join(aligns[2])+'\n')
        f.close()

if __name__ == "__main__":
    unittest.main()
