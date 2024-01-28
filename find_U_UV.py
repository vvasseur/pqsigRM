from itertools import chain

import logging

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

import pqsigrm


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


def compute_equations(GPJ):
    K, N = GPJ.dimensions()
    R = N // 2

    dim_VA = GPJ[:, :R].rank()
    dim_U = GPJ[dim_VA:, R:].rank()
    len_U2 = R - dim_U
    pivots_U = GPJ[dim_VA:, R:].pivots()
    supp_U2 = [R + i for i in range(R) if i not in pivots_U]

    # The permutation of two matched pairs in GP is in fact a linear operation
    # on GP * J. In this loop we precompute all the differences.
    equations_row = [matrix(GF(2), R, len_U2) for i in range(dim_VA)]
    for i in range(R):
        column_left = GPJ[:dim_VA, i]
        support_indices = [j for j, cj in enumerate(column_left) if cj[0] == 1]
        if i in pivots_U:
            r = pivots_U.index(i)
            row_right = GPJ[dim_VA + r, supp_U2]
            for j in support_indices:
                equations_row[j][i] = row_right
        else:
            i2 = supp_U2.index(R + i)
            for j in support_indices:
                equations_row[j][i, i2] = 1

    return equations_row, supp_U2


# Generate all the strings with exactly `u` 'U' and `v` 'V'
def generate_UV_strings(u, v, current_string="", current_u=0, current_v=0):
    if current_u == u and current_v == v:
        # If the string contains the right number of 'U's and 'V's, yield it
        yield current_string
        return

    if current_u < u:
        # If we can still add 'U's, do it
        yield from generate_UV_strings(
            u, v, current_string + "U", current_u + 1, current_v
        )

    if current_v < v:
        # If we can still add 'V's, do it
        yield from generate_UV_strings(
            u, v, current_string + "V", current_u, current_v + 1
        )


def get_pearson(sigs, start_paths=[""], u=0, v=0):
    def common_size(strings):
        if not strings:
            return None

        size = len(strings[0])
        for string in strings:
            if len(string) != size:
                return None

        return size

    if (l := common_size(start_paths)) is None:
        return None
    l += u + v

    correlation = pqsigrm.Correlation(pqsigrm.CODE_N // 2**l)
    for start_path in start_paths:
        for path in generate_UV_strings(u, v, start_path):
            sigs_a = sigs.clone()

            r = pqsigrm.CODE_R
            m = pqsigrm.CODE_M

            for dir in path:
                if dir == "U":
                    sigs_a.apply_and_half(len=2**m)
                    m -= 1
                if dir == "V":
                    sigs_a.apply_xor_half(len=2**m)
                    r -= 1
                    m -= 1
            correlation.update(sigs_a)
    pearson = correlation.finish()
    return pearson


def get_pairs_from_pearson(pearson):
    matching = min_weight_full_bipartite_matching(csr_matrix(pearson), maximize=True)
    pairs = {(min(a, b), max(a, b)) for a, b in zip(*matching)}

    return pairs


def get_permutation_from_pairs(pairs):
    permutation = list(chain.from_iterable(zip(*pairs)))

    return permutation


def apply_swaps(permutation, swaps):
    permutation_swapped = permutation[:]
    k = len(swaps)
    for i in range(0, len(permutation), 2 * k):
        for j in range(k):
            if swaps[j]:
                permutation_swapped[i + j], permutation_swapped[i + j + k] = (
                    permutation_swapped[i + j + k],
                    permutation_swapped[i + j],
                )
    return permutation_swapped


def apply_permutation_blockwise(L, permutation):
    k = len(permutation)
    permuted_L = []

    for i in range(0, len(L), k):
        block = L[i : i + k]
        permuted_block = [block[j] for j in permutation]
        permuted_L.extend(permuted_block)

    return permuted_L


def find_swaps(GP, dimA=0):
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

    equations_row, supp_U2 = compute_equations(GPJ)

    # This heuristic finds a permutation while handling the appended rows. In
    # the end, the submatrix in the upper right corner of GPJ should have a
    # rank equal to `K_APP`.
    # For each row, a linear system can be solved to find suitable column
    # swapping that cancels that row if its component on A is zero. If not, we
    # append the row to our system, hoping that it is a vector of a basis of A.
    rank = GPJ[:dim_VA, R:].rank()
    while rank > dimA:
        unsolved = []
        for j in range(dim_VA):
            if vector(GPJ[j, supp_U2]) == 0:
                continue

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
                            swaps[i] ^= 1
                    GPJ.echelonize()

        rank = GPJ[:dim_VA, R:].rank()

        if rank > dimA:
            for i in range(R):
                GPJ[:, i + R] += GPJ[:, i]
                GPJ.echelonize()
                this_rank = GPJ[:dim_VA, R:].rank()
                if this_rank < rank:
                    rank = this_rank
                    swaps[i] ^= 1
                    break
                GPJ[:, i + R] += GPJ[:, i]

    return swaps


def extract_codes(GP):
    K, N = GP.dimensions()
    R = N // 2

    J = matrix.block(
        [
            [matrix.identity(GF(2), R), matrix.identity(GF(2), R)],
            [matrix.identity(GF(2), R), 0],
        ]
    )

    GPJ = GP * J
    GPJ.echelonize()

    # The left side of GPJ has a rank equal to the dimension of V + A.
    # (A is the span of the appended rows.)
    dim_VA = GPJ[:, :R].rank()
    dim_A = GPJ[:dim_VA, R:].rank()
    dim_V = dim_VA - dim_A

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

    return A, U, V


def depth_1(sigs, G):
    sigs_current = sigs.clone()

    logging.info("Computing Pearson coefficients")
    pearson = get_pearson(sigs_current)

    logging.info("Computing matched pairs")
    pairs = get_pairs_from_pearson(pearson)
    permutation = get_permutation_from_pairs(pairs)

    logging.info("Reordering the matched pairs")
    swaps = find_swaps(G[:, permutation], pqsigrm.K_APP)
    permutation = apply_swaps(permutation, swaps)

    return permutation


def depth_2(sigs, G, permutation):
    sigs_current = sigs.clone()
    sigs_current.apply_bitwise_permutation(permutation)

    logging.info("Computing Pearson coefficients")
    pearson = get_pearson(sigs_current, ["U", "V"])

    logging.info("Computing matched pairs")
    pairs = get_pairs_from_pearson(pearson)
    permutation_small = get_permutation_from_pairs(pairs)
    permutation = apply_permutation_blockwise(permutation, permutation_small)

    logging.info("Reordering the matched pairs")
    A, G_6_12, G_5_12 = extract_codes(G[:, permutation])
    swaps = find_swaps(G_6_12)

    permutation = apply_swaps(permutation, swaps)

    return permutation


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        print("Usage: {} <pk> <sigs> <permutation>".format(sys.argv[0]))
        exit(1)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(funcName)s - %(message)s"
    )

    filename_pk = sys.argv[1]
    filename_sigs = sys.argv[2]
    filename_permutation = sys.argv[3]

    # Compute the generator matrix
    T = matrix(
        GF(2),
        pqsigrm.read_matrix(filename_pk, pqsigrm.DIM_T, pqsigrm.LEN_T),
        sparse=False,
    )
    H = matrix.block([[matrix.identity(GF(2), pqsigrm.LEN_H - pqsigrm.LEN_T), T]])
    G = dual(H)

    # Load signatures
    sigs = pqsigrm.Signatures.from_file(filename_sigs)

    # Find a global permutation that uncovers RM(5, 11)
    permutation = depth_1(sigs, G)
    permutation = depth_2(sigs, G, permutation)

    permutation_inv = sorted(range(pqsigrm.CODE_N), key=permutation.__getitem__)

    # Output the permutation in a format similar to the private key.
    pqsigrm.write_uint16s(filename_permutation, permutation_inv)
