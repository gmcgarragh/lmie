#*******************************************************************************
# You should only need to modify the first two sections on this file
#*******************************************************************************


#*******************************************************************************
# Compiler and linker settings
#*******************************************************************************
CC            = gcc
F77           = gfortran
F90           = gfortran

CCFLAGS       = -O2 -Wall -Werror -Wno-uninitialized
F77FLAGS      = -O2 -Wall -Werror -Wno-uninitialized -fno-second-underscore \
                -ffixed-line-length-none
F90FLAGS      = -O2 -Wall -Werror -Wno-uninitialized -fno-second-underscore \
                -ffree-line-length-0

# Link program that should be used for a mixed C/Fortran 
CFMIX_LD      = $(F90)
CFMIX_LDFLAGS = # -lstdc++


#*******************************************************************************
# Uncomment the following lines if parallelization with OpenMP is desired.
#*******************************************************************************
#CCDEFINES     += -DUSE_OPENMP

#CCFLAGS       += -fopenmp
#F77FLAGS      += -fopenmp
#F90FLAGS      += -fopenmp

#CFMIX_LDFLAGS += -fopenmp


#*******************************************************************************
# If the MPI version of LMie is to be compiled set the MPI include and lib
# directorires and link flags.  If these are installed in standard places
# (/usr/include, /usr/lib) then the include and library variables may be left
# blank.
#*******************************************************************************
MPI_INCDIR    = ${HOME}/opt/mpich2/include
MPI_LIBDIR    = ${HOME}/opt/mpich2/lib
MPI_LINK      = -lmpich -lmpl -lpthread


#*******************************************************************************
# 
#*******************************************************************************
VERSION       = 0.9


#*******************************************************************************
# 
#*******************************************************************************
CCDEFINES    += -D'VERSION="$(VERSION)"'

INCDIRS       = -I../src -I../interfaces -I../misc
LIBDIRS       = 

LINKS         = -lm


#*******************************************************************************
# 
#*******************************************************************************
MAKE_SUBDIRS  = target=`echo $@ | sed s/-subdirs//`; \
                for subdir in $$list; do \
                     echo "Making $$target in $$subdir"; \
                     (cd $$subdir && $(MAKE) $$target) || exit 1; \
                done;


#*******************************************************************************
# 
#*******************************************************************************
all: all-subdirs

all-subdirs:
	+ @list='$(SUBDIRS)'; $(MAKE_SUBDIRS)


mpi: mpi-subdirs

mpi-subdirs:
	+ @list='$(SUBDIRS_MPI)'; $(MAKE_SUBDIRS)


clean: clean-subdirs
	rm -f *.a *.o *.mod $(BINARIES) $(PRODUCTS)

clean-subdirs:
	+ @list='$(SUBDIRS)'; $(MAKE_SUBDIRS)


#*******************************************************************************
# 
#*******************************************************************************
.SUFFIXES: .c .f .f90 .o

.c.o:
	$(CC) $(CCFLAGS) -c $(CCDEFINES) $(INCDIRS) -o $*.o $<

.f.o:
	$(F77) $(F77FLAGS) -c $(INCDIRS) -o $*.o $<

.f90.o:
	$(F90) $(F90FLAGS) -c $(F90DEFINES) $(INCDIRS) -o $*.o $<
