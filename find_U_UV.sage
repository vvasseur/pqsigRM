import pqsigrm
import random


# Copied from
# https://groups.google.com/a/list.nist.gov/g/pqc-forum/c/3l4UmEJFi6k
def dual(mat):
    mat1 = mat.rref()
    T1 = []
    for j in range(0, mat1.nrows()):
        if mat1.row(j) == 0:
            T1.append(j)
    mat2 = mat1.delete_rows(T1)
    T3 = []
    T2 = []
    mr = mat2.rank()
    for j in range(0, mat2.ncols()):
        if j < mr and j not in mat2.pivots():
            T3.append(j)
        if j > mr - 1 and j in mat2.pivots():
            T2.append(j)
    for j in range(0, len(T2)):
        mat2.swap_columns(T3[j], T2[j])
    mat2 = mat2.rref()
    mat3 = mat2.submatrix(0, mr, mr, mat2.ncols() - mr)
    mat4 = mat3.transpose()
    i5 = matrix.identity(GF(2), mat4.nrows())
    mat5 = mat4.augment(i5)
    for j in range(0, len(T2)):
        mat5.swap_columns(T3[j], T2[j])
    return mat5


def check_structure(G, A, U, V, permutation):
    G2 = matrix.block(
        [
            [A],
            [matrix.block([[U, U], [0, V]])],
        ]
    )

    return (
        A.dimensions()[0] == pqsigrm.K_APP
        and U.dimensions()[0] == pqsigrm.DIM_U
        and V.dimensions()[0] == pqsigrm.DIM_V
        and G.rref() == G2[:, permutation].rref()
    )


def find_swaps(HP):
    GP = dual(HP)

    K, N = GP.dimensions()
    R = N // 2

    J = matrix.block(
        [
            [matrix.identity(GF(2), R), matrix.identity(GF(2), R)],
            [matrix.identity(GF(2), R), 0],
        ]
    )

    swaps = [0 for _ in range(R)]

    GPJ = GP * J
    GPJ.echelonize()

    # The left side of GPJ has a rank equal to the dimension of V + A.
    # (A is the span of the appended rows.)
    dim_VA = GPJ[:, :R].rank()

    dim_U = GPJ[dim_VA:, R:].rank()
    pivots_U = GPJ[dim_VA:, R:].pivots()
    supp_U2 = [R + i for i in range(R) if i not in pivots_U]
    len_U2 = R - dim_U

    # The permutation of two matched pairs in GP is in fact a linear operation
    # on GP * J. In this loop we precompute all the differences.
    GPJ2 = copy(GPJ)
    equations_row = [matrix(GF(2), R, len_U2) for i in range(dim_VA)]
    for i in range(R):
        if i % ((R + 19) // 20) == 0:
            print(".", file=sys.stderr, end="")

        GPJ2[:, i + R] += GPJ2[:, i]
        GPJ2.echelonize()
        diff = GPJ[:dim_VA, supp_U2] + GPJ2[:dim_VA, supp_U2]
        for j in range(dim_VA):
            equations_row[j][i] = diff[j]

        # Reverse previous operations
        GPJ2[:, i + R] += GPJ2[:, i]
    print("", file=sys.stderr)

    # This heuristic finds a permutation while handling the appended rows. In
    # the end, the submatrix in the upper right corner of GPJ should have a
    # rank equal to `K_APP`.
    # For each row, a linear system can be solved to find suitable column
    # swapping that cancels that row if its component on A is zero. If not, we
    # append the row to our system, hoping that it is a vector of a basis of A.
    rank = GPJ[:dim_VA, R:].rank()
    while rank > pqsigrm.K_APP:
        unsolved = []
        for j in range(dim_VA):
            if j % ((dim_VA + 19) // 20) == 0:
                print(".", file=sys.stderr, end="")

            A = equations_row[j]
            if unsolved:
                A = A.stack(GPJ[unsolved, supp_U2])

            try:
                sol = A.solve_left(vector(GPJ[j, supp_U2]))
            except Exception:
                unsolved.append(j)
                pass
            else:
                if sol[:R] != 0:
                    for i, pi in enumerate(sol[:R]):
                        if pi == 1:
                            GPJ[:, i + R] += GPJ[:, i]
                            swaps[i] ^^= 1
                    GPJ.echelonize()

        rank = GPJ[:dim_VA, R:].rank()

        if rank > pqsigrm.K_APP:
            for i in range(R):
                GPJ[:, i + R] += GPJ[:, i]
                GPJ.echelonize()
                this_rank = GPJ[:dim_VA, R:].rank()
                if this_rank < rank:
                    rank = this_rank
                    swaps[i] ^^= 1
                    break
                GPJ[:, i + R] += GPJ[:, i]

    dim_A = rank
    dim_V = dim_VA - dim_A
    for i, swap in enumerate(swaps):
        if swap:
            left[i], right[i] = right[i], left[i]

    permutation = left + right
    permutation_inv = sorted(range(pqsigrm.CODE_N), key=permutation.__getitem__)

    # Remove the span of the appended rows.
    S = matrix.block(
        [[GPJ[:dim_VA, R:], matrix.identity(GF(2), dim_VA)]]
    ).echelon_form()[:, R:]
    S = matrix.block([[S, 0], [0, matrix.identity(GF(2), K - dim_VA)]])
    SGPJ = S * GPJ

    # Extract the different codes.
    A = matrix.block([[SGPJ[:dim_A, R:], SGPJ[:dim_A, :R] + SGPJ[:dim_A, R:]]])
    V = SGPJ[dim_A : dim_A + dim_V, :R]
    U = SGPJ[dim_A + dim_V :, R:]
    check = check_structure(dual(H), A, U, V, permutation_inv)

    return check, permutation_inv


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {} <pk> <matched_pairs> <permutation>".format(sys.argv[0]))
        exit(1)

    filename_pk = sys.argv[1]
    filename_pairs = sys.argv[2]
    filename_permutation = sys.argv[3]

    # Import matched pairs from signatures.
    pairs = list(pqsigrm.read_matched(filename_pairs))
    pairs.sort(key=min)
    left, right = map(list, zip(*pairs))

    T = matrix(
        GF(2),
        pqsigrm.read_matrix(filename_pk, pqsigrm.DIM_T, pqsigrm.LEN_T),
        sparse=False,
    )
    H = matrix.block([[matrix.identity(GF(2), pqsigrm.LEN_H - pqsigrm.LEN_T), T]])

    # Rearrange the columns so that the matched pairs are necessarily at a
    # distance `CODE_K` from each other.
    HP = matrix(GF(2), pqsigrm.DIM_H, pqsigrm.LEN_H, sparse=False)
    for i, (a, b) in enumerate(pairs):
        HP[:, i] = H[:, a]
        HP[:, i + pqsigrm.CODE_K] = H[:, b]

    is_U_UV, permutation_inv = find_swaps(HP)
    print(is_U_UV)

    # Output the permutation in a format similar to the private key.
    pqsigrm.write_uint16s(filename_permutation, permutation_inv)
