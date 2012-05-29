include ../make.inc

OBJECTS     = lmie_core.o \
              lmie_interface.o \
              lmie_schedule.o

OBJECTS_MPI = lmie_core.o \
              lmie_interface_mpi.o \
              lmie_schedule_mpi.o

all: liblmie.a

mpi: liblmie_mpi.a

liblmie.a: $(OBJECTS)
	ar -rs liblmie.a $(OBJECTS)

liblmie_mpi.a: $(OBJECTS_MPI)
	ar -rs liblmie_mpi.a $(OBJECTS_MPI)

lmie_interface_mpi.o: lmie_interface.c
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o lmie_interface_mpi.o lmie_interface.c -DUSE_MPI -I$(MPI_INCDIR)
lmie_interface_mpi.c:

lmie_schedule.o: lmie_schedule.c
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o lmie_schedule.o lmie_schedule.c
lmie_schedule.c:

lmie_schedule_mpi.o: lmie_schedule.c
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o lmie_schedule_mpi.o lmie_schedule.c -DUSE_MPI -I$(MPI_INCDIR)
lmie_schedule_mpi.c:

include dep.inc
