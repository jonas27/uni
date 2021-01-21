
def show(matrix):
    for row in matrix:
        print(row)

def size(matrix):
    row = len(matrix)
    col = len(matrix[0][:])
    # check that all rows have same length
    [print("Error: matrix is not rectangular") for r in matrix if col != len(r)]
    return row, col