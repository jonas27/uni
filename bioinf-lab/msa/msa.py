"""

"""

from . import gotoh as g

def run(seq1, seq2, seq3, gap_open, gap_extend, substitions=None):
    if substitions == None :
        substitions = g.dna_sub
        pass
    seq, align1, align2 = bestPairwiseAlign(seq1,seq2,seq3, gap_open, gap_extend, substitions)
    a1 = (align1[0],["#" if x == "_" else x for x in align1[1]])
    a2 = (align2[0],["#" if x == "_" else x for x in align2[1]])
    a1[1].insert(0," ")
    a2[1].insert(0," ")
    g1 = a1[0], g.gotoh(seq[1], a1[1], gap_open, gap_extend,substitions)
    g2 = a2[0], g.gotoh(seq[1], a2[1], gap_open, gap_extend,substitions)
    if g1[1][0] > g2[1][0]:
        return align_final((g1[0], g1[1][1][0][1]), a2, (seq[0], g1[1][1][0][0]))
    else:
        return align_final((g2[0], g2[1][1][0][1]), a1, (seq[0], g2[1][1][0][0]))

def align_final(correct_align, expand_align, aligned):
    expand_align[1].remove(" ")
    correct_aligned = (aligned[0],['#' if x == '_' else x for x in aligned[1]])
    for i, a in enumerate(correct_align[1]):
            if a == "_":
                expand_align[1].insert(i,"_")
    alignments = [correct_align,expand_align,aligned]
    alignments.sort(key=lambda x: x[0])
    aligns = []
    for align in alignments:
        aligns.append(['#' if x == '_' else x for x in align[1]])
    return aligns

def bestPairwiseAlign(seq1,seq2,seq3,gap_open,gap_extend, substitions):
    """
    return tuple of seq number and alignment to keep order.
    """
    ab = g.gotoh(seq1,seq2,gap_open, gap_extend, substitions)
    ac = g.gotoh(seq1,seq3,gap_open, gap_extend, substitions)
    bc = g.gotoh(seq2,seq3,gap_open, gap_extend, substitions)
    # sort here and then return guide alignment
    # alignments = [ab,ac,bc]
    # alignments.sort(key=lambda x: -x[0])



    m = max([ab[0],ac[0],bc[0]])
    if m == ab[0]:
        # aligns = []
        # in1 = ab[1][0][0]
        # in1.insert(0," ")
        # for a in bc[1]:
        #     in2 = a[0]
        #     in2.insert(0," ")
        #     aligns.append(g.gotoh(in1, in2,gap_open,gap_extend, substitions))
        return (2,seq3), (0,ab[1][0][0]),(1,ab[1][0][1])
    elif m == ac[0]:
        return (1,seq2), (0,ac[1][0][0]),(2,ac[1][0][1])
    else:
        return (0,seq1), (1,bc[1][0][0]),(2,bc[1][0][1])