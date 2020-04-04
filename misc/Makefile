include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS = gindex_name_value.o \
          gutil_alloc.o       \
          gutil_error.o       \
          gutil_parse.o       \
          rtutil_math.o       \
          rtutil_scat_int.o   \
          rtutil_scat_io.o    \
          rtutil_size_dist.o  \
          rtutil_support.o    \
          zeroin.o

all: liblmie_misc.a

liblmie_misc.a: $(OBJECTS)
	ar -rs liblmie_misc.a $(OBJECTS)

include dep.inc
