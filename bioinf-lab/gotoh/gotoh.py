# Read and append empty sapce in beginning
def read_fasta_file(fasta_file):
    line = open(fasta_file).readlines()[1].replace("\n","")
    sequence = list(" " + line)
    return sequence

def read_fasta_file_test():
    sequence = read_fasta_file("data/s1.fasta")
    if ''.join(sequence) != ' ILDMDVVEGSAARFDCKVEGYPDPEVMWFKDDNPVKESRHFQIDYDEEGN':
        print("Error reading fasta file.")

def score_of_alignment(align_seq1, align_seq2, cost_gap_open, 
                       cost_gap_extension, substitutions=None):
    """
    A nice helper function which computes the score of the given alignment.
    This is only used for the self check.
    Input example:
    --CG
    AACA
    """
    return score

def read_substitution_matrix(file_substitution_matrix):
    """
    Implement reading the scores file.
    It can be stored as a dictionary of example: scores[("A", "R")] = -1
    """
    lines = open(file_substitution_matrix).readlines()
    matrix = []
    # Convert to matrix
    firstline = True
    for l in lines:
        if l[0]=='#':
            continue
        if firstline :
            l = l.replace("   ","-  ").replace(" ","").replace("\n","")
            matrix.append(list(l))
            firstline=False
        else:
            line=[]
            line.append(list(l)[0])
            istr = [x.strip('') for x in ''.join(list(l)[3:]).replace("\n","").replace("  ", " ").split(" ")]
            line.extend(istr)
            matrix.append(line)

    scores = dict()
    for row in range(1,len(matrix)):
        for col in range(1,len(matrix[0])):
            scores[(matrix[0][col],matrix[row][0])] = int(matrix[row][col])
    return scores

def read_substitution_matrix_test():
    scores = read_substitution_matrix("data/pam250.txt")
    if scores.get(('*','*')) != 1 or scores.get(('R','A')) != -2:
        print("Error in converting to dict")


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
    for i in range(1,len(seq1)):
        row = [None if x != ' ' else cost_open + counter * cost_extend for x in seq2]
        counter+=1
        matrix.append(row)

    return matrix
def initd_test():
    m = initd([" ","C","G"],[" ","C","C","G","A"], -3, -1)
    show(m)
    size(m)



def initp(seq_1, seq_2):
    matrix = []
    # Assumption: we never reach a score of -1,000,000 --> is infinity!
    col = [-10**6 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else "-" for x in seq_2]
        matrix.append(col)
    return matrix

def initp_test():
    # seq1, seq2 = read_fasta_file("data/s1.fasta"),read_fasta_file("data/s2.fasta")
    # result = init_matrix_p(seq1,seq2)
    seq1, seq2 = [" ","A","B"],[" ","A","B","C","D","E"]
    result = initp(seq1,seq2)
    show(result)


def initq(seq_1, seq_2):
    matrix = []
    # Assumption: we never reach a score of -1,000,000 --> is infinity!
    col = [-10**10 if x != ' ' else 0 for x in seq_2]
    matrix.append(col)
    for i in range(1,len(seq_1)):
        col = [-10**10 if x != ' ' else -10**6 for x in seq_2]
        matrix.append(col)
    return matrix

def initq_test():
    seq1, seq2 = [" ","A","B"],[" ","A","B","C","4"]
    result = initq(seq1,seq2)
    show(result)
    size(result)


def show(matrix):
    for row in matrix:
        print(row)

def size(matrix):
    row = len(matrix)
    col = len(matrix[0][:])
    for c in matrix:
        if len(c) != col:
            print(matrix)
            raise Exception('Matrix not rectangular')
    return row, col

def show_size_test():
    seq1, seq2 = [" ","A","B"],[" ","A","B","C","4"]
    result = initp(seq1,seq2)
    show(result)
    size(result)


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

def com(a,b):
    return 1 if a == b else -1

def complete_d_p_q_computation_test():
    seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
    d,p,q = complete_d_p_q_computation(seq1,seq2,-3,-1,com)
    if d[1] != [-4, 1, -3, -4, -5] or p[2] != ['-', -3, -7, -8, -9] or q[3] != [-1000000, -10, -8, -8, -3]:
        raise Exception('matrix computation is wrong!') 


"""
You are working with 3 matrices simultaneously.
You can store your path as a list of cells.
A cell can be a tuple: coordinates, matrix_name.
And coordinates is a tuple of indexex i, j.

Cell example: ((0, 2), "d")
Path example: [((2, 4), 'd'), ((2, 4), 'q'), ((2, 3), 'q'), ((2, 2), 'd'), ((1, 1), 'd'), ((0, 0), 'd')]

"""
# seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
# d,p,q = complete_d_p_q_computation(seq1,seq2,-3,-1,com)
def compute_tracebacks(seq1, seq2, d, p, q,
                           cost_open, cost_extend, substitution=None):
    """
    Implement a search for all possible paths from the bottom right corner to the top left.
    Implement 'find_all_previous' and check_complete first.
   
    """
    il,jl=size(d)
    previous = find_previous(((il-1,jl-1), 'd'),seq1, seq2, d, p, q,cost_open, cost_extend, substitution)
    paths = [previous]
    for x in range(1000):
        li = []
        # print("preve")
        # print(previous)
        # print(li)
        for pr in previous:
            # print(pr)
            print("pr")
            print(pr)
            pre = find_previous(pr ,seq1, seq2, d, p, q,cost_open, cost_extend, substitution)
            for fi in pre:
                if fi[0] != (0,0):
                    li.append(pre[0])
        # print("li")
        # print(li)
        if len(li) == 0:
            break
        previous = li
        # print("preve2")
        # print(previous)
        paths.append(li)

    return paths

def compute_tracebacks_test():
    seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
    d,p,q = complete_d_p_q_computation(seq1,seq2,-3,-1,com)
    p = compute_tracebacks(seq1, seq2, d,p,q,-3,-1,com)
    print("tracebacks")
    print(p)
    
    
def find_previous(cell, seq1, seq2, d, p, q,
                   cost_open, cost_extend, substitution=None):
    parent_cells = []
    """
    Implement a search for all possible previous cells.
    """
    i,j, = cell[0]
    curr = None
    if cell[1] == 'd':
        curr = d[i][j]
    if cell[1] == 'p':
        curr = p[i][j]
    if cell[1] == 'q':
        curr = q[i][j]
    # print(max([d[i-1][j-1]+substitution(seq1[i],seq2[j])]))
    # if second last field we go to d in last field 
    if cell[0] == (0,1) or cell[0] == (1,0):
        return [((0,0),'d')]
    
    if curr ==  max([d[i-1][j-1]+substitution(seq1[i],seq2[j]),p[i][j],q[i][j]]):
        parent_cells.append(((i-1,j-1),'d'))
    if curr == max([d[i-1][j]+cost_open+cost_extend,p[i-1][j]+cost_extend]):
        parent_cells.append(((i-1,j-1),'p'))
    if curr == max([d[i][j-1]+cost_open+cost_extend,q[i][j-1]+cost_extend]):
        parent_cells.append(((i,j),'q'))
    # print("parent cell")
    # print(parent_cells)
    return parent_cells

def find_previous_test():
    seq1, seq2 = [" ","C","G","G"],[" ","C","C","G","A"]
    d,p,q = complete_d_p_q_computation(seq1,seq2,-3,-1,com)
    p = find_previous(((3,4), "d"), seq1, seq2, d,p,q,-3,-1,com)
    if p[0][1] != 'd':
        print(p[0][1])
        raise Exception('find previous is wrong!') 

def check_complete(path):
    """
    Implement a function which checks if the traceback path is complete.
    """


if __name__ == "__main__":
    # find_previous_test()
    compute_tracebacks_test()
    # compute_all_tracebacks(seq1,seq2,d,p,q,-3,-1)