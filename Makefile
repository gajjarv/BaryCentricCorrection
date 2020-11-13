FC = gfortran
GC = gcc

BINDIR = ./bin
LIBDIR = ./include
SIGPROC="/Users/vishalgajjar/prog/bl_sigproc/src/"
GSL="/usr/local/Cellar/gsl/1.16/lib/"

CFLAGS1= -g -Wno-unused-result -O3 -DHAVE_CONFIG_H -I. -I$(LIBDIR) -I/usr/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS2= -L/usr/lib $(SIGPROC)/libsigproc.a -lm -L${GSL} -lgsl -lgslcblas -lz  

all: barycentre_seti.o barycentre_seti

barycentre_seti.o: barycentre_seti.c
	$(GC) $(CFLAGS1) -c barycentre_seti.c 

barycentre_seti: barycentre_seti.o
	$(GC) $(CFLAGS2) -o barycentre_seti barycentre_seti.o 

clean: 
	rm *.o	barycentre_seti


