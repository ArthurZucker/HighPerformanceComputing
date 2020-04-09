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
	mpirun -n 20  -hostfile hostfile --map-by node sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/hood.mtx.gz  | ./cg " > /dev/null
#mpirun -n 2 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd2.mtx > ./x.txt
debug:
	mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/cfd1.mtx
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
