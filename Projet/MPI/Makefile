# your choice of compiler
CC = mpicc
SEED = 42
# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g
LDLIBS = -lm
MATRIX = cfd2
all : cg_challenge
cg_challenge : cg_challenge.o mmio.o
mmio.o : mmio.c mmio.h
cg_challenge.o : cg_challenge.c mmio.h

run :
	mpirun -n 40 -mca bt1_tcp_if_include eth0 -hostfile hostfile305 --map-by node sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/${MATRIX}.mtx.gz  | ./cg  --seed ${SEED} "> x.txt
run2:
	mpirun -n 68 -hostfile hostfile305 ./cg_challenge > x.txt
debug:
	mpirun -np 40 gdb -ex run --args ./cg --matrix ../Matrix/${MATRIX}.mtx --hostfile hostfile305
check:
	./checker --matrix ../Matrix/${MATRIX}.mtx --seed ${SEED} --solution x.txt
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
