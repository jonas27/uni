"""
gotoh has a main function 'gotoh' and a convience class 'Gotoh'.
The class only has one method which in turn just runs the gotoh method.
(I guess we should have this. I would rather delete the class.)
"""

def gotoh(seq1, seq2, cost_open, cost_extension, substition=None):
    """gotoh returns a score and the alignments for two sequences.
    seq1 and seq2: The sequences to be compared to.
          They need to be char lists with an empty starting character.
    cost_open: The costs to open a new gap as an integer.
    cost_extension: The costs to extend an existing gap as an integer.
    substition: substition is a function returning a comparison value
          of two nucleotides/amino acids.
    The function takes as an input a tuple of two chars, ie. substition(('E','E'))."""
    if substition is None:
        substition = dna_sub
    d,p,q = complete_d_p_q_computation(seq1,seq2, cost_open, cost_extension,substition)
    tracebacks = compute_tracebacks(seq1, seq2, d,p,q, cost_open,cost_extension,substition)
    alignments = []
    for t in tracebacks:
        alignments.append(build_alignment(seq1,seq2,t))
    score = score_of_alignment(alignments[0][0],alignments[0][1],
        cost_open,cost_extension,substition)
    return score, alignments

class Gotoh:
    """Gotoh is a convinience class with only the run method."""
    def run(self, seq1, seq2, cost_open, cost_extension, substition=None):
        """run only delegates to the gotoh() function."""
        return gotoh(seq1, seq2, cost_open, cost_extension, substition)

def score_of_alignment(align_seq1, align_seq2, cost_open,
                       cost_extension, substitution):
    """ score_of_alignment computes the score of a single alignment.
    Input example:
    __CG
    AACA
    """
    score = 0
    for i in range(max([len(align_seq1),len(align_seq2)])):
        if align_seq1[i] == '_' or align_seq2[i] == '_':
            if align_seq1[i-1] == '_' or align_seq2[i-1] == '_':
                score+=cost_extension
            else:
                score+=cost_open + cost_extension
        else:
            score+=substitution((align_seq1[i],align_seq2[i]))
    return score

def initd(seq1, seq2, cost_open, cost_extend):
    """initd initialized the D-Matrix of gotoh"""
    matrix, counter, row = [], 1, [0]
    for x in range(1,len(seq2)):
        if x != ' ':
            x = cost_open + counter * cost_extend
            row.append(x)
            counter+=1
    matrix.append(row)
    counter = 1
    for _ in range(1,len(seq1)):
        row = [None if x != ' ' else cost_open + counter * cost_extend for x in seq2]
        counter+=1
        matrix.append(row)
    return matrix

def initp(seq_1, seq_2):
    """initp initialized the P-Matrix of gotoh"""
    matrix = []
    col = [-10**6 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for _ in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**10 for x in seq_2]
        matrix.append(col)
    return matrix

def initq(seq_1, seq_2):
    """initq initialized the Q-Matrix of gotoh"""
    matrix = []
    col = [-10**10 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for _ in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**6 for x in seq_2]
        matrix.append(col)
    return matrix

def complete_d_p_q_computation(seq1, seq2, cost_open, cost_extend, substitution=None):
    """complete_d_p_q completes all three matrixes"""
    d = initd(seq1,seq2,cost_open, cost_extend)
    p = initp(seq1,seq2)
    q = initq(seq1,seq2)
    # all matrices have the same size, so iterate over entire matrix
    for i in range(1,size(d)[0]):
        for j in range(1,size(d)[1]):
            p[i][j] = max([d[i-1][j]+cost_open+cost_extend,p[i-1][j]+cost_extend])
            q[i][j] = max([d[i][j-1]+cost_open+cost_extend,q[i][j-1]+cost_extend])
            d[i][j] = max([d[i-1][j-1]+substitution((seq1[i],seq2[j])),p[i][j],q[i][j]])
    return d,p,q

def dna_sub(ins):
    """dna_sub is the dna substition function.
    ins: a tuple of Nucleoids and it returns 1 if equal, -1 else."""
    return 1 if ins[0] == ins[1] else -1

def compute_tracebacks(seq1, seq2, d, p, q,
                           cost_open, cost_extend, substitution):
    """compute_tracebacks computes all tracebacks."""
    il,jl=size(d)
    paths = []
    cells = find_previous(((il-1,jl-1), 'd',0),seq1, seq2,
        d, p, q,cost_open, cost_extend, [], substitution)
    path = [((il-1,jl-1), 'd',0)]
    for pos in range(1,len(cells)):
        if cells[pos-1][2] == cells[pos][2]-1:
            path.append(cells[pos])
        else:
            paths.append(path)
            path = path[:cells[pos][2]]
            path.append(cells[pos])
    paths.append(path)
    rev = [x[::-1] for x in paths]
    return rev

def find_previous(cell, seq1, seq2, d, p, q,
                   cost_open, cost_extend, history, substitution):
    """find_previous recersevly finds the a path from to (0,0) to (il,jl).
    It returns an array of all possible aligments with the form
        one step inside the trace has the form:
        ((row,col),matrix, step_counter)"""
    parent_cells = []
    i,j = cell[0]
    curr = None
    step_counter = cell[2]+1
    history.append(cell)
    if cell[0] == (0,1) or cell[0] == (1,0) or (cell[0] == (1,1) and cell[1] == 'd'):
        # if second last field we go to d in last field
        history.append(((0,0),'d', cell[2]+1))
        return history
    if cell[0][0] == 0 and cell[1] == 'd':
        # This catches the case where we open a gap in first position
        parent_cells.append(find_previous(((i,j-1),'d', step_counter),
        seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    if cell[0][1] == 0 and cell[1] == 'd':
        parent_cells.append(find_previous(((i-1,j),'d', step_counter),
        seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    # This builds
    if cell[1] == 'd':
        curr = d[i][j]
        if curr ==  d[i-1][j-1]+substitution((seq1[i],seq2[j])):
            parent_cells.append(find_previous(((i-1,j-1),'d', step_counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i][j]:
            parent_cells.append(find_previous(((i,j),'p', step_counter) ,
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j]:
            parent_cells.append(find_previous(((i,j),'q', step_counter) ,
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'p':
        curr = p[i][j]
        if curr ==  d[i-1][j]+cost_open+cost_extend:
            parent_cells.append(find_previous(((i-1,j),'d', step_counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i-1][j]+cost_extend:
            parent_cells.append(find_previous(((i-1,j),'p', step_counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'q':
        curr = q[i][j]
        if curr ==  d[i][j-1]+cost_open+cost_extend:
            parent_cells.append(find_previous(((i,j-1),'d', step_counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j-1]+cost_extend:
            parent_cells.append(find_previous(((i,j-1),'q', step_counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    return history

def build_alignment(seq1, seq2, traceback):
    """build_alignments builds the alignment for two sequences and a particular traceback.
    it returns the two resulting alignment sequences and the signs.
    To see the signs between the sequences use the print-order (s1,signs,s2)."""
    s1, s2, signs = [], [], []
    traceback = traceback[1:]
    c1,c2 = 0,0
    for cell in traceback:
        if cell[0][0] == c1 and cell[0][1] == c2:
            pass
        elif cell[0][0] == c1:
            s1.append('_')
            s2.append(seq2[cell[0][1]])
            signs.append(' ')
        elif cell[0][1] == c2:
            s2.append('_')
            s1.append(seq1[cell[0][0]])
            signs.append(' ')
        else:
            s1.append(seq1[cell[0][0]])
            s2.append(seq2[cell[0][1]])
            if seq1[cell[0][0]] ==seq2[cell[0][1]]:
                signs.append('*')
            else:
                signs.append('|')
        c1,c2 = cell[0]
    return s1,s2,signs

def size(matrix):
    """size is a convience method to get the size of the matrix
    it returns (len(rows), len(cols)).
    Prints Error: if matrix is not recatangular, i.e. len(row[i]) != len(row[j])"""
    rows = len(matrix)
    row = len(matrix[0][:])
    # check that all rows have same length
    for r in matrix:
        if len(r) != row:
            raise Exception('Matrix not rectangular.')
    return rows, row
    