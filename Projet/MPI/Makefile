# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -Wall -Wextra -g
LDLIBS = -lm

all : cg_same_A

cg_same_A : cg_same_A.o mmio.o
mmio.o : mmio.c mmio.h
cg_same_A.o : cg_same_A.c mmio.h
#mpirun -n 4 -hostfile hostfile --map-by node zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/bcsstk13.mtx.gz | ./cg_same_A
run :
	mpirun -n 2 -hostfile hostfile --map-by node ./cg_same_A --matrix ../Matrix/cfd1.mtx > /dev/null
debug:
	mpirun -np 2 gdb -ex run --args ./cg_same_A --matrix ../Matrix/bcsstk13.mtx
.PHONY: clean
clean :
	rm -rf *.o cg_same_A *.receipt *.json
