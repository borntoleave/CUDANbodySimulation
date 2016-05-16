#Define the compilers
CC=g++
#Define the compiler options
CFLAGS=-O3

#Define the dependencies
nbody: nbody.o
	$(CC) -o nbody nbody.o

nbody.o: nbody.cpp
	$(CC) -c $(CFLAGS) nbody.cpp

###
clean:
	rm -f nbody nbody.o
