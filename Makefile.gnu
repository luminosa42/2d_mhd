.SUFFIXES : .f90 .c

.f90.o :
	gfortran -O3 -fdefault-real-8 -c $*.f90

#.c.o :
#	gcc -O3 -c $*.c

driver2d : driver2d.o init.o plm1d.o god1d.o mhd1d.o riemann.o
	gfortran -O3 -o driver2d driver2d.o init.o plm1d.o god1d.o mhd1d.o riemann.o

clean :
	rm driver2d *.o
