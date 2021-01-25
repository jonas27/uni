# gotoh returns a score and the alignments for two sequences.
# seq1 and seq2: The sequences to be compared to. They need to be char lists with an empty starting character.
# cost_open: The costs to open a new gap as an integer.
# cost_extension: The costs to extend an existing gap as an integer.
# substition: substition is a function returning a comparison value of two nucleotides/amino acids.
# The function takes as an input a tuple of two chars, ie. substition(('E','E')).
def gotoh(seq1, seq2, cost_open, cost_extension, substition=None):
    # seq1, seq2 = read_fasta_file(fasta_file_1), read_fasta_file(fasta_file_2)
    if substition == None:
        substition = dna_sub  
    d,p,q = complete_d_p_q_computation(seq1,seq2, cost_open, cost_extension,substition)
    tracebacks = compute_tracebacks(seq1, seq2, d,p,q, cost_open,cost_extension,substition)
    alignments = []
    for t in tracebacks:
        alignments.append(build_alignment(seq1,seq2,t))
    score = score_of_alignment(alignments[0][0],alignments[0][1],cost_open,cost_extension,substition)
    return score, alignments

# Gotoh is a convinience class.
# Its only function is run, which returns the gotoh() function.
class Gotoh:
    def run(self, seq1, seq2, cost_open, cost_extension, substition=None):
        return gotoh(seq1, seq2, cost_open, cost_extension, substition=None)

# score_of_alignment computes the score of a single alignment.
# 
def score_of_alignment(align_seq1, align_seq2, cost_open, 
                       cost_extension, substitution):
    """
    A nice helper function which computes the score of the given alignment.
    This is only used for the self checks.
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

# initd initialized the D-Matrix of gotoh
def initd(seq1, seq2, cost_open, cost_extend):
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

# initp initialized the P-Matrix of gotoh
def initp(seq_1, seq_2):
    matrix = []
    col = [-10**6 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**10 for x in seq_2]
        matrix.append(col)
    return matrix

# initq initialized the Q-Matrix of gotoh
def initq(seq_1, seq_2):
    matrix = []
    col = [-10**10 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**6 for x in seq_2]
        matrix.append(col)
    return matrix

def complete_d_p_q_computation(seq1, seq2, cost_open, cost_extend, substitution=None):
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

# dna_sub is the dna substition function.
# ins: a tuple of Nucleoids and it returns 1 if equal, -1 else.
def dna_sub(ins):
    return 1 if ins[0] == ins[1] else -1

# compute_tracebacks 
def compute_tracebacks(seq1, seq2, d, p, q,
                           cost_open, cost_extend, substitution):
    il,jl=size(d)
    paths = []
    cells = find_previous(((il-1,jl-1), 'd',0),seq1, seq2, d, p, q,cost_open, cost_extend, [], substitution)
    path = [((il-1,jl-1), 'd',0)]
    for p in range(1,len(cells)):
        if cells[p-1][2] == cells[p][2]-1:
            path.append(cells[p])
        else:
            paths.append(path)
            path = path[:cells[p][2]]
            path.append(cells[p])
    paths.append(path)
    rev = [x[::-1] for x in paths]
    return rev

# find_previous recersevly finds the a path from to (0,0) to (il,jl).
# It returns an array of all possible aligments with the form
# return[0] 
def find_previous(cell, seq1, seq2, d, p, q,
                   cost_open, cost_extend, history, substitution):
    parent_cells = []
    """
    Implement a search for all possible previous cells.
    """
    # if second last field we go to d in last field 
    i,j = cell[0]
    curr = None
    counter = cell[2]+1
    history.append(cell)
    if cell[0] == (0,1) or cell[0] == (1,0) or (cell[0] == (1,1) and cell[1] == 'd'):
        history.append(((0,0),'d', cell[2]+1))
        return history
    # This catches the case where we open a gap in first position
    elif cell[0][0] == 0 and cell[1] == 'd':
        parent_cells.append(find_previous(((i,j-1),'d', counter),
        seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    elif cell[0][1] == 0 and cell[1] == 'd':
        parent_cells.append(find_previous(((i-1,j),'d', counter),
        seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    # This builds 
    if cell[1] == 'd':
        curr = d[i][j]
        if curr ==  d[i-1][j-1]+substitution((seq1[i],seq2[j])): 
            parent_cells.append(find_previous(((i-1,j-1),'d', counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i][j]:
            parent_cells.append(find_previous(((i,j),'p', counter) ,
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j]:
            parent_cells.append(find_previous(((i,j),'q', counter) ,
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'p':
        curr = p[i][j]
        if curr ==  d[i-1][j]+cost_open+cost_extend: 
            parent_cells.append(find_previous(((i-1,j),'d', counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i-1][j]+cost_extend:
            parent_cells.append(find_previous(((i-1,j),'p', counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'q':
        curr = q[i][j]
        if curr ==  d[i][j-1]+cost_open+cost_extend: 
            parent_cells.append(find_previous(((i,j-1),'d', counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j-1]+cost_extend:
            parent_cells.append(find_previous(((i,j-1),'q', counter),
            seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    return history

# build_alignments builds the alignment for two sequences and a particular traceback.
# it returns the two resulting alignment sequences and the signs.
# To see the signs between the sequences use the print-order (s1,signs,s2).
def build_alignment(seq1, seq2, traceback):
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
    row = len(matrix)
    col = len(matrix[0][:])
    # check that all rows have same length
    [print("Error: matrix is not rectangular") for r in matrix if col != len(r)]
    return row, col