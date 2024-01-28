from ctypes import (
    CDLL,
    POINTER,
    Structure,
    addressof,
    c_double,
    c_int64,
    c_size_t,
    c_ubyte,
    c_uint64,
    c_uint8,
    c_void_p,
    cast,
)
from numpy.ctypeslib import as_array


CODE_N = 8192
MLEN = 32
CRYPTO_BYTES = (64 + 8192 + 64) // 8

double = c_double
int64_t = c_int64
size_t = c_size_t
uint64_t = c_uint64
uint8_t = c_ubyte
uint8_t = c_uint8
void_p = c_void_p

Array2D = (double * CODE_N) * CODE_N
Sigs_Array = uint8_t * (MLEN + CRYPTO_BYTES)

_lib = CDLL("./pqsigrm613/signatures.so")


class Tables(Structure):
    _fields_ = [
        ("n", uint64_t),
        ("sum_X", POINTER(int64_t)),
        ("sum_XY", POINTER(int64_t)),
    ]


class Signatures:
    def __init__(self, raw, ptr, n):
        self.raw = raw
        self.ptr = ptr
        self.n = n

    def clone(self):
        new_raw = (type(self.raw))(*self.raw)

        offset = addressof(self.ptr.contents) - addressof(self.raw)
        new_ptr_address = addressof(new_raw) + offset
        new_ptr = cast(new_ptr_address, POINTER(type(self.ptr.contents)))

        return Signatures(new_raw, new_ptr, self.n)

    @classmethod
    def from_file(cls, filename, MLEN=32):
        with open(filename, "rb") as f:
            file_data = f.read()

        ArrayType = uint8_t * len(file_data)
        data_array = ArrayType.from_buffer_copy(file_data)

        # Get pointer to the ctypes array
        data_ptr = cast(data_array, POINTER(uint8_t))

        # Skip the length of the message and the message
        shifted_ptr_address = addressof(data_ptr.contents) + 8 + MLEN
        shifted_ptr = cast(void_p(shifted_ptr_address), POINTER(Sigs_Array))

        n = len(file_data) // (MLEN + CRYPTO_BYTES)

        return cls(data_array, shifted_ptr, n)

    def apply_bitwise_permutation(self, permutation, len=CODE_N):
        permutation_array = (size_t * len)(*permutation)

        _lib.apply_bitwise_permutation(self.ptr, permutation_array, len, self.n)

    def apply_and_half(self, len):
        _lib.apply_and_half(self.ptr, len, self.n)

    def apply_xor_half(self, len):
        _lib.apply_xor_half(self.ptr, len, self.n)


class Correlation:
    def __init__(self, length):
        self.length = length
        self.table_ptr = _lib.init_tables(length)

    def update(self, sigs):
        _lib.update_correlation(self.table_ptr, sigs.ptr, self.length, sigs.n)

    def finish(self):
        pearson = Array2D()
        pearson_ptr = POINTER(Array2D)(pearson)

        _lib.finish_correlation(pearson_ptr, self.table_ptr, self.length)
        pearson = as_array(pearson)[: self.length, : self.length]
        return pearson

    def __del__(self):
        _lib.free_tables(self.table_ptr)


_lib.apply_bitwise_permutation.argtypes = [
    POINTER(Sigs_Array),
    POINTER(size_t),
    size_t,
    size_t,
]

_lib.apply_xor_half.argtypes = [
    POINTER(Sigs_Array),
    size_t,
    size_t,
]

_lib.apply_and_half.argtypes = [
    POINTER(Sigs_Array),
    size_t,
    size_t,
]

_lib.init_tables.argtypes = [size_t]
_lib.init_tables.restype = POINTER(Tables)

_lib.free_tables.argtypes = [POINTER(Tables)]
_lib.free_tables.restype = None

_lib.update_correlation.argtypes = [
    POINTER(Tables),
    POINTER(Sigs_Array),
    size_t,
    size_t,
]
_lib.update_correlation.restype = None

_lib.finish_correlation.argtypes = [POINTER(Array2D), POINTER(Tables), size_t]
_lib.finish_correlation.restype = None
