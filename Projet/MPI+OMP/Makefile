# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -Wall -Wextra -g -fopenmp
LDLIBS = -lm -fopenmp

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h
#mpirun -n 4 -hostfile hostfile --map-by node zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/bcsstk13.mtx.gz | ./cg
run :
	mpirun -n 4 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd1.mtx > /dev/null
debug:
	mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/cfd1.mtx
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
