import pqsigrm
import sys


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} <secret key> <matched_pairs>".format(sys.argv[0]))
        exit(1)

    filename_sk = sys.argv[1]
    filename_pairs = sys.argv[2]

    Q = pqsigrm.get_uint16s(filename_sk, pqsigrm.CODE_N)

    real_matched_pairs = set()
    for i in range(pqsigrm.CODE_K):
        a = Q.index(i)
        b = Q.index(i + pqsigrm.CODE_K)
        real_matched_pairs.add(frozenset({a, b}))

    matched_pairs = pqsigrm.read_matched(filename_pairs)
    matched_pairs = {frozenset(pair) for pair in matched_pairs}

    print(
        "Found {} matched pairs out of {}".format(
            len(matched_pairs & real_matched_pairs), len(real_matched_pairs)
        )
    )
