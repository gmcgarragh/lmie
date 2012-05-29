include ../make.inc

BINARIES = calllmie calllmie_mpi

all: calllmie

mpi: calllmie_mpi

calllmie: calllmie.o ../src/liblmie.a ../misc/liblmie_misc.a
	$(CC) $(CCFLAGS) -o calllmie calllmie.o ../src/liblmie.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)

calllmie_mpi: calllmie_mpi.o ../src/liblmie_mpi.a ../misc/liblmie_misc.a
	$(CC) $(CCFLAGS) -o calllmie_mpi calllmie_mpi.o ../src/liblmie_mpi.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS) -L$(MPI_LIBDIR) $(MPI_LINK)

calllmie_mpi.o: calllmie.c
	$(CC) $(CCFLAGS) $(CCDEFINES) -c $(INCDIRS) -o calllmie_mpi.o calllmie.c -DUSE_MPI -I$(MPI_INCDIR)

include dep.inc
