FC = gfortran
GC = gcc

BINDIR = ./bin
LIBDIR = ./include
SIGPROC="/home/vgajjar/BaryCentricCorrection_sigproc/bl_sigproc/src/"
CFITSIO="/usr/lib/x86_64-linux-gnu/"
GSL="/usr/local/Cellar/gsl/1.16/lib/"

CFLAGS1= -g -Wno-unused-result -O3 -DHAVE_CONFIG_H -I. -I$(LIBDIR) -I/usr/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 
CFLAGS2= -L/usr/lib -L$(SIGPROC) -lsigproc -lm -L${GSL} -lgsl -lgslcblas -lz  -L${CFITSIO} -lcfitsio

all: barycentre_seti.o barycentre_seti 

barycentre_seti.o: barycentre_seti.c
	$(GC) -c barycentre_seti.c $(CFLAGS1) 

barycentre_seti: barycentre_seti.o
	$(GC) -o barycentre_seti barycentre_seti.o $(CFLAGS2) 

clean: 
	rm *.o	barycentre_seti


