# your choice of compiler
CC = gcc
SEED = 42
# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g -fopenmp
LDLIBS = -lm -fopenmp
MATRIX = cfd1
all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

run :
	./cg --matrix ../Matrix/${MATRIX}.mtx --seed ${SEED} > x.txt

check:
	./checker --matrix ../Matrix/${MATRIX}.mtx --seed ${SEED} --solution x.txt
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
