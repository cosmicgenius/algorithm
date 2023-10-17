# given an n by m matrix, find a subset of rows that sum to 0 in F_2
# outputs list of linearly independent subsets (hopefully the maximum amount)
# Algorithm from: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
def generate_subsets(matrix):
    n = len(matrix)
    if n == 0:
        return None

    m = len(matrix[0])
    if m == 0:
        return None
    
    for r in range(n):
        if len(matrix[r]) != m:
            return None

    marked = [False] * n
    row_marked = []
    for c in range(m):
        pivot = None
        for r in range(n):
            if matrix[r][c] % 2 == 1:
                pivot = r
                break
         
        row_marked.append(pivot)
        if pivot != None:
            marked[pivot] = True
            for c2 in range(m):
                if matrix[pivot][c2] % 2 == 1 and c2 != c:
                    for r in range(n):
                        matrix[r][c2] = (matrix[r][c2] + matrix[r][c]) % 2

    all_subsets = []
    for r in range(n):
        if not marked[r]:
            needed_rows = [r]
            for c in range(m):
                if matrix[r][c] % 2 == 1:
                    needed_rows.append(row_marked[c])
            all_subsets.append(needed_rows)
    return all_subsets

if __name__ == '__main__':
    m = [[1, 1, 0, 0],
         [1, 1, 0, 1],
         [0, 4, 1, 1],
         [2, 0, 1, 0],
         [0, 2, 0, 1]]
    print(generate_subsets(m))
