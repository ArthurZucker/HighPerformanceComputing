# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -Wall -Wextra -g
LDLIBS = -lm

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

run :mpirun -n 4 -hostfile hostfile --map-by node ./cg --matrix /Infos/lmd/2019/master/ue/MU4IN903-2020fev/bcsstk13.mtx > /dev/null

.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
