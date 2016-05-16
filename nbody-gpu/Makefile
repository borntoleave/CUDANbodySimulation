#Define the compilers
NVCC=nvcc
#Define the compiler options
CFLAGS=-O3 -arch sm_20

#Define the dependencies
nbody: nbody.o
	$(NVCC) -o nbody nbody.o

nbody.o: nbody.cu
	$(NVCC) -c $(CFLAGS) nbody.cu

###
clean:
	rm -f nbody nbody.o
