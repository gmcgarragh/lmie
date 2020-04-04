include ../make.inc
include ../common.inc

OBJECTS = lmie_core.o \
          lmie_interface.o \
          lmie_schedule.o

OBJECTS_PTHREADS = lmie_core.o \
                   lmie_interface.o \
                   lmie_schedule_pthreads.o

all: liblmie.a liblmie_pthreads.a

liblmie.a: $(OBJECTS)
	ar -rs liblmie.a $(OBJECTS)

liblmie_pthreads.a: $(OBJECTS_PTHREADS)
	ar -rs liblmie_pthreads.a $(OBJECTS_PTHREADS)

lmie_schedule_pthreads.o: lmie_schedule.c
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o lmie_schedule_pthreads.o lmie_schedule.c -DUSE_PTHREADS

include dep.inc
