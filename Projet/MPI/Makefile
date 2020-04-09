# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -Wall -Wextra -g
LDLIBS = -lm

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

run :
	mpirun -n 8 -hostfile nodes --map-by node sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/cfd1.mtx.gz | ./cg" > /dev/null
debug :
	mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/cfd1.mtx
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
