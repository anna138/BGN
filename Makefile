.PHONY: clean

CFLAGS  := -Wall -O3 -g -std=c++0x -D_GNU_SOURCE
CC      := g++
LDLIBS  := -fopenmp -lgmp -lsodium -lm -lpbc
INC     := -I /usr/local/include -I /usr/local/include/pbc

APPS    := runMain runMainopenmp

all: ${APPS}

multRun: ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp mainmulttest.cpp
	${CC} -o $@ $^ ${LDLIBS} ${CFLAGS} ${INC}

addRun:  ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp mainaddtest.cpp
	${CC} -o $@ $^ ${LDLIBS} ${CFLAGS} ${INC}

runMain: ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp main.cpp
	${CC} -o $@ $^ ${LDLIBS} ${CFLAGS} ${INC}

runMainopenmp: ciphertext.cpp bgn.cpp plaintext.cpp parser.cpp main.cpp
	${CC} -o $@ $^ ${LDLIBS} ${CFLAGS} -DPARAL=1 -fopenmp ${INC}

clean:
	rm -f *.o *.enc output.txt ${APPS}
