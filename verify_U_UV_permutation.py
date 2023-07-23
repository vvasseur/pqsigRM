import pqsigrm
import sys


def check_permutation(Q, P):
    return (
        (
            set([P[Q.index(i)] for i in range(pqsigrm.CODE_N // 2)])
            == set(range(pqsigrm.CODE_N // 2))
        )
        and (
            set([P[Q.index(i)] for i in range(pqsigrm.CODE_N // 2, pqsigrm.CODE_N)])
            == set(range(pqsigrm.CODE_N // 2, pqsigrm.CODE_N))
        )
        and all(
            [
                P[Q.index(i + pqsigrm.CODE_N // 2)]
                == P[Q.index(i)] + pqsigrm.CODE_N // 2
                for i in range(pqsigrm.CODE_N // 2)
            ]
        )
    )


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <secret key> <permutation>".format(sys.argv[0]))
        exit(1)

    filename_sk = sys.argv[1]
    filename_perm = sys.argv[2]

    Q = pqsigrm.get_uint16s(filename_sk, pqsigrm.CODE_N)
    P = pqsigrm.get_uint16s(filename_perm, pqsigrm.CODE_N)

    print(check_permutation(Q, P))
