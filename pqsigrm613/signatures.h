#ifndef __SIGNATURES_H
#define __SIGNATURES_H
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include "src/api.h"
#include "src/parm.h"

#define MLEN 32

typedef struct tables {
    uint64_t n;
    int64_t *sum_X;
    int64_t *sum_XY;
} tables;

tables *init_tables(size_t len);
void free_tables(tables *t);
void update_correlation(tables *t, uint8_t (*e)[MLEN + CRYPTO_BYTES],
                        size_t len, size_t n);
void finish_correlation(double (*pearson)[CODE_N], tables *t, size_t len);

void apply_bitwise_permutation(unsigned char (*e)[MLEN + CRYPTO_BYTES],
                               size_t *permutation, size_t len, size_t n);
void apply_and_half(unsigned char (*e)[MLEN + CRYPTO_BYTES], size_t len,
                    size_t n);
void apply_xor_half(unsigned char (*e)[MLEN + CRYPTO_BYTES], size_t len,
                    size_t n);
#endif
