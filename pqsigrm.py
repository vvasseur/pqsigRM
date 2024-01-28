from pqsigrm_ctypes import Correlation, Signatures

import struct

CODE_K = 4096
CODE_N = 8192

DIM_H = CODE_N - CODE_K - 1
LEN_H = CODE_N

# T is the non-systematic part of H.
DIM_T = CODE_N - CODE_K - 1
LEN_T = CODE_K + 1

K_APP = 2
DIM_U = 2508
DIM_V = 1587

CODE_R = 6
CODE_M = 13


def read_matrix(filename, nrows, ncols):
    with open(filename, "rb") as f:
        byte_data = f.read()

    colsize = (ncols + 63) // 64 * 64

    bit_matrix = [
        [
            (byte_data[i // 8] >> (i % 8)) & 1
            for i in range(j * colsize, j * colsize + ncols)
        ]
        for j in range(nrows)
    ]

    return bit_matrix


def read_matched(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    tuples = [tuple(int(x) for x in line.split()) for line in lines]
    return set((a, b) if a < b else (b, a) for (a, b) in tuples)


def get_uint16s(filename, len):
    with open(filename, "rb") as file:
        data = file.read(len * 2)
        numbers = struct.unpack("H" * len, data)
        return numbers


def write_uint16s(filename, values):
    with open(filename, "wb") as file:
        data = struct.pack("H" * len(values), *values)
        file.write(data)
