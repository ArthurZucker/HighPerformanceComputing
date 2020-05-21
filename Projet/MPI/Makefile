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
	#mpirun -n 4 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd1.mtx > ./x.txt
	#mpirun -n 20  -hostfile hostfile --map-by node sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/hood.mtx.gz  | ./cg " > /dev/null
	#mpirun -n 2 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd1.mtx > /dev/null
	#mpirun -np 2 gdb -ex run --args ./cg --matrix ../Matrix/cfd1.mtx
	#mpirun -n 4 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/cfd1.mtx > ./x.txt
	#mpirun -n 70 -mca bt1_tcp_if_include eth1 -hostfile hostfile408 --map-by node sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/nd24k.mtx.gz  | ./cg " > /dev/null
	mpirun -n 68 -hostfile hostfile --map-by node ./cg --matrix ../Matrix/${MATRIX}.mtx --seed ${SEED} > x.txt
run2:
	mpirun -n 136 -mca btl_tcp_if_include eth1 -hostfile hostfile406 ./cg_challenge > x.txt
debug:
	mpirun -np 40 gdb -ex run --args ./cg --matrix ../Matrix/${MATRIX}.mtx --hostfile hostfile305
check:
	./checker --matrix ../Matrix/${MATRIX}.mtx --seed ${SEED} --solution x.txt
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
