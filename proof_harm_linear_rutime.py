import sys

def print_matrix(matrix, seeds_total):
    print()
    print("matrix for", seeds_total, "seeds:")
    for row in matrix:
        for a, b in row:
            print("| s", seeds_total - a, ", ", b, sep="", end="\t")
        print("|")

def compute_tree(stack, num_seeds):
    def add_matrices(m1, m2):
        m = []
        for x, row in enumerate(m1):
            m.append([])
            for y, ele in enumerate(row):
                a, b = ele
                a_, b_ = m2[x][y]
                m[x].append( (a + a_, b + b_) )
        return m
    print(".", end="")
    sys.stdout.flush()
    matrix = [[(0,0)]*num_seeds]*num_seeds
    # no contradictions
    stack_1 = stack + num_seeds
    for x, row in enumerate(compute_tree(stack_1, num_seeds-1)):
        for y, ele in enumerate(row):
            a_, b_ = matrix[x][y]
            a, b = ele
            matrix[x][y] = (a + a_, b + b_)

    # contradictions

    #return value
    return matrix


num_seeds = 5
print_matrix(compute_tree([], num_seeds), num_seeds)
