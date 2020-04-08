# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -Wall -Wextra -g
LDLIBS = -lm

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h
#mpirun -n 4 -hostfile hostfile --map-by node zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/bcsstk13.mtx.gz | ./cg
run :
	mpirun -n 2 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd1.mtx > /dev/null
debug:
<<<<<<< HEAD
	mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/cfd1.mtx
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
=======
	mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/bcsstk13.mtx
.PHONY: clean
clean :
	rm -rf *.o cg cg *.receipt *.json
>>>>>>> 9980e2f71e88196b364754da7c7a6bcccbf0d171
