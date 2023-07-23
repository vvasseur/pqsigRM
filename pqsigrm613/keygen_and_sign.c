#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "src/api.h"
#include "src/rng.h"

#define RUN_SUCCESS 0
#define RUN_FILE_OPEN_ERROR -1
#define RUN_DATA_ERROR -3
#define RUN_CRYPTO_FAILURE -4

#define MLEN 32

void write_bin(char *filename, unsigned char *data, size_t length) {
    FILE *file = fopen(filename, "wb");
    if (file != NULL) {
        fwrite(data, sizeof(unsigned char), length, file);
        fclose(file);
    }
    else {
        printf("Failed to open file: %s\n", filename);
    }
}

int keygen_and_sign(size_t n, char *basename) {
    // Add extensions to the basename.
    char filename_sigs[256], filename_sk[256], filename_pk[256];
    sprintf(filename_sigs, "%s%s", basename, "-sigs");
    sprintf(filename_sk, "%s%s", basename, "-sk");
    sprintf(filename_pk, "%s%s", basename, "-pk");

    unsigned char entropy_input[48];
    unsigned char *m, *sm;
    uint64_t smlen;
    unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    int ret_val;

    for (int i = 0; i < 48; i++)
        entropy_input[i] = i;

    randombytes_init(entropy_input, NULL, 256);

    // Generate the public/private keypair.
    if ((ret_val = crypto_sign_keypair(pk, sk)) != 0) {
        printf("crypto_sign_keypair returned <%d>\n", ret_val);
        return RUN_CRYPTO_FAILURE;
    }

    m = (unsigned char *)calloc(MLEN, sizeof(unsigned char));
    if (m == NULL) {
        printf("Failed to allocate memory for the message buffer.\n");
        return -1;
    }

    sm = (unsigned char *)calloc(n * (MLEN + CRYPTO_BYTES),
                                 sizeof(unsigned char));
    if (sm == NULL) {
        free(m);
        printf("Failed to allocate memory for the signature buffer.\n");
        return -1;
    }

    write_bin(filename_pk, pk, CRYPTO_PUBLICKEYBYTES);
    write_bin(filename_sk, sk, CRYPTO_SECRETKEYBYTES);

    // Generate `n` signatures.
    for (size_t i = 0; i < n; i++) {
        if (i % (n / 20) == 0)
            fprintf(stderr, ".");

        randombytes(m, MLEN);

        if ((ret_val =
                 crypto_sign(&sm[i * (MLEN + CRYPTO_BYTES)],
                             (unsigned long long *)&smlen, m, MLEN, sk)) != 0) {
            printf("crypto_sign returned <%d>\n", ret_val);
            return RUN_CRYPTO_FAILURE;
        }
    }
    fprintf(stderr, "\n");
    write_bin(filename_sigs, sm, n * (MLEN + CRYPTO_BYTES));

    free(m);
    free(sm);
    return RUN_SUCCESS;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <N> <basename>\n", argv[0]);
        return 1;
    }

    size_t n = strtoul(argv[1], NULL, 10); // Number of signatures
    char *basename = argv[2];              // Basename for outputs

    return keygen_and_sign(n, basename);
}
