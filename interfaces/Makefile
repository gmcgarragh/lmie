include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS = lmie_int_f77.o \
          lmie_int_f90.o \
          lmie_int_util.o

all: liblmie_interfaces.a

liblmie_interfaces.a: $(OBJECTS)
	ar -rs liblmie_interfaces.a $(OBJECTS)

include dep.inc
