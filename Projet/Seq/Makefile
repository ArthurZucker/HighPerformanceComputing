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
<<<<<<< HEAD
	./cg --matrix ../Matrix/cfd2.mtx > /dev/null
=======
	sh -c "zcat /Infos/lmd/2019/master/ue/MU4IN903-2020fev/hood.mtx.gz | ./cg" > /dev/null
>>>>>>> c93c87ff2c6676390c9d5bfb8a3c723ea845ab69
.PHONY: clean
clean :
	rm -rf *.o cg *.receipt *.json
