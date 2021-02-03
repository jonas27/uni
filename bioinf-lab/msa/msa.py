"""
msa is an implementation for a 3 sequences aligment.
it exports the run method.
"""

from gotoh import gotoh, dna_sub

def run(seq1, seq2, seq3, gap_open, gap_extend, substitions=None):
    """
    run is the main method, used to find a msa.
    """
    if substitions is None :
        substitions = dna_sub
    seq, align1, align2 = _bestPairwiseAlign(seq1,seq2,seq3, gap_open, gap_extend, substitions)
    a1 = (align1[0],["#" if x == "_" else x for x in align1[1]])
    a2 = (align2[0],["#" if x == "_" else x for x in align2[1]])
    a1[1].insert(0," ")
    a2[1].insert(0," ")
    g1 = a1[0], gotoh(seq[1], a1[1], gap_open, gap_extend,substitions)
    g2 = a2[0], gotoh(seq[1], a2[1], gap_open, gap_extend,substitions)
    if g1[1][0] > g2[1][0]:
        return _align_final((g1[0], g1[1][1][0][1]), a2, (seq[0], g1[1][1][0][0]))
    return _align_final((g2[0], g2[1][1][0][1]), a1, (seq[0], g2[1][1][0][0]))

def _align_final(correct_align, expand_align, aligned):
    expand_align[1].remove(" ")
    for i, a in enumerate(correct_align[1]):
        if a == "_":
            expand_align[1].insert(i,"_")
    alignments = [correct_align,expand_align,aligned]
    alignments.sort(key=lambda x: x[0])
    aligns = []
    for align in alignments:
        aligns.append(['#' if x == '_' else x for x in align[1]])
    return aligns

def _bestPairwiseAlign(seq1,seq2,seq3,gap_open,gap_extend, substitions):
    """
    return tuple of seq number and alignment to keep order.
    """
    ab = gotoh(seq1,seq2,gap_open, gap_extend, substitions)
    ac = gotoh(seq1,seq3,gap_open, gap_extend, substitions)
    bc = gotoh(seq2,seq3,gap_open, gap_extend, substitions)
    m = max([ab[0],ac[0],bc[0]])
    if m == ab[0]:
        return (2,seq3), (0,ab[1][0][0]),(1,ab[1][0][1])
    if m == ac[0]:
        return (1,seq2), (0,ac[1][0][0]),(2,ac[1][0][1])
    return (0,seq1), (1,bc[1][0][0]),(2,bc[1][0][1])
