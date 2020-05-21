# your choice of compiler
CC = clang

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g -fopenmp -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include
LDLIBS = -lm -fopenmp=libomp

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

run :
	./cg --matrix ../Matrix/cfd1.mtx > /dev/null
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
