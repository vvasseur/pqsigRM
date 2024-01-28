#include <stdio.h>
#include "signatures.h"
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
    double value;
} pair;

static int pair_comp(const void *pa, const void *pb) {
    double a = ((pair *)pa)->value;
    double b = ((pair *)pb)->value;
    return (a < b) - (a > b);
}

int compute_matched_pairs(size_t *pairs, double pearson[CODE_N][CODE_N],
                          size_t len) {
    for (size_t i = 0; i < len; ++i) {
        pair dict[len];
        for (size_t k = 0; k < len; ++k) {
            dict[k].key = k;
            dict[k].value = pearson[i][k];
        }
        qsort(dict, len, sizeof(pair), pair_comp);
        pairs[i] = dict[0].key;
    }
    int ok = 1;
    for (size_t i = 0; i < len; ++i) {
        ok = ok && (pairs[pairs[i]] == i);
    }
    return ok;
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

    double(*pearson)[CODE_N] = calloc(CODE_N * CODE_N, sizeof(double));
    if (pearson == NULL) {
        printf("Failed to allocate memory for the correlation table.\n");
        free(sigs);
        return -1;
    }

    size_t n = size / (8 + MLEN + CRYPTO_BYTES);

    // Only consider the error part of the signatures.
    unsigned char(*e)[MLEN + CRYPTO_BYTES] =
        (unsigned char(*)[MLEN + CRYPTO_BYTES]) & sigs[8 + MLEN];
    size_t pairs[CODE_N];

    // Stage 1
    tables *t = init_tables(CODE_N);

    update_correlation(t, e, CODE_N, n);
    finish_correlation(pearson, t, CODE_N);
    fprintf(stderr, "Correlations done\n");

    // Output the matched pairs.
    int ok = compute_matched_pairs(pairs, pearson, CODE_N);
    fprintf(stderr, "Matched pairs %d\n", ok);
    for (size_t i = 0; i < CODE_N; ++i) {
        printf("%ld %ld\n", i, pairs[i]);
    }

    free_tables(t);
    free(pearson);

    free(sigs);

    return 0;
}
