#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "src/api.h"
#include "src/parm.h"

#define MLEN 32

unsigned char *read_bin(const char *filename, size_t *size) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        printf("Failed to open file: %s\n", filename);
        return NULL;
    }

    // Get the size of the file.
    fseek(file, 0, SEEK_END);
    *size = ftell(file);
    rewind(file);

    // Allocate memory for the buffer.
    unsigned char *buffer = (unsigned char *)malloc(*size);
    if (buffer == NULL) {
        printf("Failed to allocate memory for the file buffer.\n");
        fclose(file);
        return NULL;
    }

    // Read the file into the buffer.
    size_t read_size = fread(buffer, 1, *size, file);
    if (read_size != *size) {
        printf("Failed to read file: %s\n", filename);
        free(buffer);
        fclose(file);
        return NULL;
    }

    fclose(file);

    return buffer;
}

typedef struct pair {
    uint64_t key;
    uint64_t value;
} pair;

static int pair_comp(const void *pa, const void *pb) {
    int64_t a = ((pair *)pa)->value;
    int64_t b = ((pair *)pb)->value;
    return (a > b) - (a < b);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    size_t size;

    unsigned char *sigs = read_bin(argv[1], &size);
    if (sigs == NULL) {
        printf("Failed to allocate memory for the signature buffer.\n");
        return -1;
    }

    uint64_t(*correlation)[CODE_N] = calloc(CODE_N * CODE_N, sizeof(uint64_t));
    if (correlation == NULL) {
        printf("Failed to allocate memory for the correlation table.\n");
        free(sigs);
        return -1;
    }

    size_t n = size / (MLEN + CRYPTO_BYTES);

    // Only consider the error part of the signatures.
    unsigned char *e = &sigs[8 + MLEN];

    for (size_t i = 0; i < n; i++) {
        if (i % (n / 20) == 0)
            fprintf(stderr, ".");

        for (size_t j = 0; j < CODE_N; ++j) {
            uint64_t ej = e[j / 8] >> (j % 8) & 1;
            for (size_t k = j + 1; k < CODE_N; ++k) {
                uint64_t ek = e[k / 8] >> (k % 8) & 1;
                correlation[j][k] += ej ^ ek;
            }
        }
        e += MLEN + CRYPTO_BYTES;
    }
    // Correlations are symmetric so we only computed half of them in the loops
    // above.
    for (size_t j = 0; j < CODE_N; ++j) {
        for (size_t k = j + 1; k < CODE_N; ++k) {
            correlation[k][j] = correlation[j][k];
        }
    }
    fprintf(stderr, "\n");

    // Output the matched pairs.
    for (size_t i = 0; i < CODE_N; ++i) {
        pair dict[CODE_N];
        for (size_t k = 0; k < CODE_N; ++k) {
            dict[k].key = k;
            dict[k].value = correlation[i][k];
        }
        qsort(dict, CODE_N, sizeof(pair), pair_comp);

        printf("%ld %ld\n", i, dict[1].key);
    }

    free(correlation);
    free(sigs);

    return 0;
}
