# type of compiler
CC = gcc
# compiler flags:
#  -g        adds debugging information to the executable file.
#  -c        compile or assemble the source files, but do not link.
#  -Wall     turns on most, but not all, compiler warnings.
#  -lm       linking C library not a part of the standard C libraries.
#  -fopenmp  activate the OpenMP extensions for thread parallelization.
CFLAGS1  = -g -c -Wall
CFLAGS2 = -lm -fopenmp

boucles: boucles.o file_management.o
	$(CC) -g -o boucles boucles.o file_management.o $(CFLAGS2)

boucles.o: boucles.c file_management.h
	$(CC) $(CFLAGS1) boucles.c $(CFLAGS2)

file_management.o: file_management.c file_management.h
	$(CC) $(CFLAGS1) file_management.c $(CFLAGS2)

clean: 
	rm -f boucles *.o