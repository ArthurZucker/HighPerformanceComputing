# your choice of compiler
CC = gcc

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g
LDLIBS = -lm

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

run :
	./cg --matrix ../Matrix/cfd2.mtx > /dev/null
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
