include ../make.inc
include ../common.inc

BINARIES = calllmie calllmie_pthreads

all: $(BINARIES)

calllmie: calllmie.o ../src/liblmie.a ../misc/liblmie_misc.a
	$(CC) $(CCFLAGS) -o calllmie calllmie.o ../src/liblmie.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)

calllmie_pthreads: calllmie.o ../src/liblmie_pthreads.a ../misc/liblmie_misc.a
	$(CC) $(CCFLAGS) -o calllmie_pthreads calllmie.o ../src/liblmie_pthreads.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)

include dep.inc
