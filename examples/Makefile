include ../make.inc
include ../common.inc

BINARIES = example_c example_f77 example_f90

all: $(BINARIES)

example_c: example_c.o ../src/liblmie.a ../misc/liblmie_misc.a
	$(CC) $(CCFLAGS)  -o example_c example_c.o ../src/liblmie.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)

example_f77: example_f77.o ../interfaces/liblmie_interfaces.a ../src/liblmie.a ../misc/liblmie_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_f77 example_f77.o ../interfaces/liblmie_interfaces.a ../src/liblmie.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)

example_f90: example_f90.o ../interfaces/liblmie_interfaces.a ../src/liblmie.a ../misc/liblmie_misc.a
	$(CFMIX_LD) $(CFMIX_LDFLAGS) -o example_f90 example_f90.o ../interfaces/liblmie_interfaces.a ../src/liblmie.a ../misc/liblmie_misc.a \
        $(LIBDIRS) $(LINKS)


include dep.inc
