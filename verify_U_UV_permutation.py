import pqsigrm
import sys


def check_permutation(Q, P, depth):
    delta = 8192 // 2**depth
    tests = [i * delta for i in range(2**depth)]

    Q_Pinv = [Q[P.index(i)] for i in range(pqsigrm.CODE_N)]

    return all(
        [[Q_Pinv[i + j] - Q_Pinv[i] for j in tests] == tests for i in range(delta)]
    )


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <secret key> <permutation>".format(sys.argv[0]))
        exit(1)

    filename_sk = sys.argv[1]
    filename_perm = sys.argv[2]

    Q = pqsigrm.get_uint16s(filename_sk, pqsigrm.CODE_N)
    P = pqsigrm.get_uint16s(filename_perm, pqsigrm.CODE_N)

    print(check_permutation(Q, P, 2))
