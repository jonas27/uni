
# Read and append empty sapce in beginning
def read_fasta_file(fasta_file):
    file = open(fasta_file)
    line = file.readlines()[1].replace("\n","")
    sequence = list(" " + line)
    # Close file to prevent mem leak. Sadly python has no defer.
    file.close()
    return sequence

def read_substitution_matrix(file_substitution_matrix):
    """
    Implement reading the scores file.
    It can be stored as a dictionary of example: scores[("A", "R")] = -1
    """
    file = open(file_substitution_matrix)
    lines = file.readlines()
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
    # Close file to prevent mem leak
    file.close()
    return scores

def show(matrix):
    for row in matrix:
        print(row)