def gotoh(seq1, seq2, cost_open, cost_extension, file_substitution_matrix=None):
    # seq1, seq2 = read_fasta_file(fasta_file_1), read_fasta_file(fasta_file_2)
    if file_substitution_matrix == None:
        file_substitution_matrix = dna_sub  
    d,p,q = complete_d_p_q_computation(seq1,seq2, cost_open, cost_extension,file_substitution_matrix)
    tracebacks = compute_tracebacks(seq1, seq2, d,p,q, cost_open,cost_extension,file_substitution_matrix)
    alignments = []
    for t in tracebacks:
        alignments.append(build_alignment(seq1,seq2,t))
    score = score_of_alignment(alignments[0][0],alignments[0][1],cost_open,cost_extension,file_substitution_matrix)
    return score, alignments

class Gotoh:
    def run(fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix=None):

        alignment_score, alignments = [],[]
        return alignment_score, alignments


def score_of_alignment(align_seq1, align_seq2, cost_open, 
                       cost_extension, substitution=None):
    """
    A nice helper function which computes the score of the given alignment.
    This is only used for the self check.
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
            score+=substitution(align_seq1[i],align_seq2[i])
    return score

def initd(seq1, seq2, cost_open, cost_extend):
    """
    Implement initialization of the matrix D
    """
    matrix = []
    counter=1
    row = [0]
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
    matrix = []
    # Assumption: we never reach a score of -1,000,000 --> is infinity!
    col = [-10**6 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**10 for x in seq_2]
        matrix.append(col)
    return matrix

def initq(seq_1, seq_2):
    matrix = []
    # Assumption: we never reach a score of -1,000,000 --> is infinity!
    col = [-10**10 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**6 for x in seq_2]
        matrix.append(col)
    return matrix

def complete_d_p_q_computation(seq1, seq2, cost_open, cost_extend, substitution=None):
    """
    Implement the recursive computation of matrices D, P and Q
    """
    d = initd(seq1,seq2,cost_open, cost_extend)
    p = initp(seq1,seq2)
    q = initq(seq1,seq2)
    # d = [[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15]]
    # all matrices have the same size, so iterate over entire matrix
    for i in range(1,size(d)[0]):
        for j in range(1,size(d)[1]):
            p[i][j] = max([d[i-1][j]+cost_open+cost_extend,p[i-1][j]+cost_extend])
            q[i][j] = max([d[i][j-1]+cost_open+cost_extend,q[i][j-1]+cost_extend])
            # cost = 1 if seq1[i] == seq2[j] else -1
            d[i][j] = max([d[i-1][j-1]+substitution(seq1[i],seq2[j]),p[i][j],q[i][j]])
    return d,p,q

def dna_sub(a,b):
    return 1 if a == b else -1


"""
You are working with 3 matrices simultaneously.
You can store your path as a list of cells.
A cell can be a tuple: coordinates, matrix_name.
And coordinates is a tuple of indexex i, j.

Cell example: ((0, 2), "d")
Path example: [((2, 4), 'd'), ((2, 4), 'q'), ((2, 3), 'q'), ((2, 2), 'd'), ((1, 1), 'd'), ((0, 0), 'd')]

"""
def compute_tracebacks(seq1, seq2, d, p, q,
                           cost_open, cost_extend, substitution=None):
    """
    Implement a search for all possible paths from the bottom right corner to the top left.
    Implement 'find_all_previous' and check_complete first.
   
    """
    il,jl=size(d)
    cells = find_previous(((il-1,jl-1), 'd',0),seq1, seq2, d, p, q,cost_open, cost_extend, [], substitution)
    paths = []
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
    
def find_previous(cell, seq1, seq2, d, p, q,
                   cost_open, cost_extend, history, substitution=None):
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
        parent_cells.append(find_previous(((i,j-1),'d', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    elif cell[0][1] == 0 and cell[1] == 'd':
        parent_cells.append(find_previous(((i-1,j),'d', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        return history
    
    if cell[1] == 'd':
        curr = d[i][j]
        if curr ==  d[i-1][j-1]+substitution(seq1[i],seq2[j]): 
            parent_cells.append(find_previous(((i-1,j-1),'d', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i][j]:
            parent_cells.append(find_previous(((i,j),'p', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j]:
            parent_cells.append(find_previous(((i,j),'q', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'p':
        curr = p[i][j]
        if curr ==  d[i-1][j]+cost_open+cost_extend: 
            parent_cells.append(find_previous(((i-1,j),'d', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == p[i-1][j]+cost_extend:
            parent_cells.append(find_previous(((i-1,j),'p', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    if cell[1] == 'q':
        curr = q[i][j]
        if curr ==  d[i][j-1]+cost_open+cost_extend: 
            parent_cells.append(find_previous(((i,j-1),'d', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
        if curr == q[i][j-1]+cost_extend:
            parent_cells.append(find_previous(((i,j-1),'q', counter) ,seq1, seq2, d, p, q,cost_open, cost_extend, history,substitution))
    return history

def build_alignment(seq1, seq2, traceback):
    s1, s2, signs = [], [], []
    traceback = traceback[1:]
    c1,c2 = 0,0
    for cell in traceback:
        # we never stay in 0,0 so this isn't needed. Left as a reminder. 
        # TODO: delete soon, for visuals
        if cell[0][0] == c1 and cell[0][1] == c2:
            # s1.append('_')
            # s2.append('_')
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
        


def visualize_traceback_test():
    seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
    traceback = [((0, 0), 'd', 4), ((0, 1), 'd', 3), ((1, 2), 'd', 2), ((2, 3), 'd', 1), ((3, 4), 'd', 0)]
    s1,s2,signs = visualize_traceback(seq1, seq2, traceback)
    print(s1)
    print(signs)
    print(s2)
    traceback = [((0, 0), 'd', 5), ((1, 1), 'd', 4), ((1, 2), 'q', 3), ((1, 2), 'd', 2), ((2, 3), 'd', 1), ((3, 4), 'd', 0)]
    s1,s2,signs = visualize_traceback(seq1, seq2, traceback)
    print(s1)
    print(signs)
    print(s2)
    traceback = [((0, 0), 'd', 5), ((1, 1), 'd', 4), ((2, 2), 'd', 3), ((3, 3), 'd', 2), ((3, 4), 'q', 1), ((3, 4), 'd', 0)]
    s1,s2,signs = visualize_traceback(seq1, seq2, traceback)
    print(s1)
    print(signs)
    print(s2)

def size(matrix):
    row = len(matrix)
    col = len(matrix[0][:])
    # check that all rows have same length
    [print("Error: matrix is not rectangular") for r in matrix if col != len(r)]
    return row, col