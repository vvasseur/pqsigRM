CC = gcc

LDFLAGS =  -L/usr/local/lib
CFLAGS = -I/usr/local/include -Wunused-variable -Wunused-function -mavx2
LIBFLAGS = -lcrypto -lssl -lm
CFILES := $(shell find src -name '*.c' | sed -e 's/\.c/\.o/')

OBJS = ${CFILES}

ifeq ($(DEBUG), 1)
    DBG_FLAGS = -g -O0 -DDEBUG
else
    DBG_FLAGS = -O3 -march=native -fopenmp -DNDEBUG -Wunused-variable -Wunused-function
endif


%.o : %.c
	$(CC) $(CFLAGS) $(DBG_FLAGS) -o $@ -c $< 
	
sign: ${OBJS} ./PQCgenKAT_sign.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^ $(LIBFLAGS)

runtime: ${OBJS} ./runtime.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^ $(LIBFLAGS)

keygen_and_sign: ${OBJS} ./keygen_and_sign.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^ $(LIBFLAGS)

get_matched_pairs: ${OBJS} ./get_matched_pairs.c ./signatures.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^ $(LIBFLAGS)

matrix.o : matrix.h
rng.o : rng.h
api.o : api.h

clean:
	cd src; rm -f *.o; cd ..
	rm -f *.o
	rm -f sign runtime keygen_and_sign get_matched_pairs
