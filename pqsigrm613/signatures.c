#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "signatures.h"

static uint64_t expand(uint8_t b) {
    uint64_t shift = 0x0000040810204081ul; // bits set: 0, 7, 14, 21, 28, 35, 42
    uint64_t mask = 0x0001010101010101ul;  // bits set: 0, 8, 16, 24, 32, 40, 48
    return (uint64_t)(b & 0x7f) * shift & mask ^ (uint64_t)(b & 0x80) << 49;
}

tables *init_tables(size_t len) {
    int64_t *sum_X = calloc(len, sizeof(int64_t));
    if (!sum_X)
        return NULL;

    int64_t *sum_XY = calloc(len * (len - 1) / 2, sizeof(int64_t));
    if (!sum_XY)
        goto free_sum_X;

    tables *t = malloc(sizeof(tables));
    if (!t)
        goto free_sum_XY;

    t->sum_X = sum_X;
    t->sum_XY = sum_XY;
    t->n = 0;

    return t;

free_sum_XY:
    free(sum_XY);
free_sum_X:
    free(sum_X);
    return NULL;
}

void free_tables(tables *t) {
    if (t != NULL) {
        free(t->sum_X);
        free(t->sum_XY);
        free(t);
    }
}

void update_correlation(tables *t, uint8_t (*e)[MLEN + CRYPTO_BYTES],
                        size_t len, size_t n) {
#pragma omp parallel
    {
        tables *lt = init_tables(len);

#pragma omp for schedule(dynamic, 128)
        for (size_t i = 0; i < n; ++i) {
            uint8_t ei[len];

            for (size_t j = 0; j < (len + 7) / 8; ++j) {
                ((uint64_t *)ei)[j] = expand(e[i][j]);
            }

            size_t idx = 0;
            for (size_t j = 0; j < len; ++j) {
                uint8_t ej = ei[j];
                lt->sum_X[j] += 1 - 2 * ej;
                for (size_t k = j + 1; k < len; ++k) {
                    uint8_t ek = ei[k];
                    lt->sum_XY[idx] += 1 - 2 * (ej ^ ek);
                    ++idx;
                }
            }
        }
#pragma omp critical
        {
            for (size_t j = 0; j < len; ++j) {
                t->sum_X[j] += lt->sum_X[j];
            }
            for (size_t idx = 0; idx < len * (len - 1) / 2; ++idx) {
                t->sum_XY[idx] += lt->sum_XY[idx];
            }
        }
        free_tables(lt);
    }
    t->n += n;
}

void finish_correlation(double (*pearson)[CODE_N], tables *t, size_t len) {
    for (size_t j = 0; j < len; ++j) {
        pearson[j][j] = 0;
    }
    for (size_t j = 0; j < len; ++j) {
        size_t j_idx = j * (2 * len - j - 3) / 2 - 1;
        for (size_t k = j + 1; k < len; ++k) {
            double mean_XY = (double)t->sum_XY[j_idx + k] / t->n;
            double mean_X = (double)t->sum_X[j] / t->n;
            double mean_Y = (double)t->sum_X[k] / t->n;
            pearson[j][k] =
                (mean_XY - mean_X * mean_Y) /
                (sqrt((t->n - mean_X * mean_X) * (t->n - mean_Y * mean_Y)));
        }
    }
    // Correlations are symmetric so we only computed half of them in the loops
    // above.
    for (size_t j = 0; j < len; ++j) {
        for (size_t k = j + 1; k < len; ++k) {
            pearson[k][j] = pearson[j][k];
        }
    }
}

void apply_bitwise_permutation(unsigned char (*e)[MLEN + CRYPTO_BYTES],
                               size_t *permutation, size_t len, size_t n) {
    bool done[len];
    for (size_t i = 0; i < len; i++) {
        done[i] = false;
    }

    for (size_t i = 0; i < len; i++) {
        if (!done[i]) {
            size_t j = i;
            while (permutation[j] != i) {
                for (size_t k = 0; k < n; k++) {
                    unsigned char src_bit = (e[k][j / 8] >> (j % 8)) & 1;
                    unsigned char dest_bit =
                        (e[k][permutation[j] / 8] >> (permutation[j] % 8)) & 1;

                    // Swap the bits in place
                    e[k][j / 8] ^= (src_bit ^ dest_bit) << (j % 8);
                    e[k][permutation[j] / 8] ^= (src_bit ^ dest_bit)
                                                << (permutation[j] % 8);
                }

                done[j] = true;
                j = permutation[j];
            }
            done[j] = true;
        }
    }
}

void apply_xor_half(unsigned char (*e)[MLEN + CRYPTO_BYTES], size_t len,
                    size_t n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < len / 2; ++j) {
            size_t k = j + len / 2;
            uint64_t ek = e[i][k / 8] >> (k % 8) & 1;
            e[i][j / 8] ^= ek << (j % 8);
        }
    }
}

void apply_and_half(unsigned char (*e)[MLEN + CRYPTO_BYTES], size_t len,
                    size_t n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < len / 2; ++j) {
            size_t k = j + len / 2;
            uint64_t ek = e[i][k / 8] >> (k % 8) & 1;
            uint64_t ej = e[i][j / 8] >> (j % 8) & 1;
            if (ej && !ek)
                e[i][j / 8] ^= 1 << (j % 8);
        }
    }
}
