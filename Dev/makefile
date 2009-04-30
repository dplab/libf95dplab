# $Id: makefile,v 1.45 2001/08/22 18:00:07 balay Exp $ 

CFLAGS	         = -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE}
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/

default:  linear-solver.x

include ${PETSC_DIR}/bmake/common/base

#
MAIN.o = isaiah.o
MYOBJS = LinearSolverPetscModule.o $(OBJECTS)
MYLIB  = libf90.a 
#
#
#  My f90 options.
#
LIBDIR   = ../../f90-dplab/modules/BUILD
MODDIR   = $(LIBDIR)
LIBRARY  = $(LIBDIR)/libf90.a
LIBFLAGS = -L$(LIBDIR) -lf90 
MODFLAGS = -module $(MODDIR)
#
OTHER_OPTS = -error-limit 5 -check all -traceback -check noarg_temp_created
#
FFLAGS           = $(MODFLAGS) $(OTHER_OPTS)
LOCDIR		 = .

#FC_FLAGS =  -I. -fPIC -names lowercase
#FLINKER= mpif90 -I. -fPIC 
#PETSC_FORTRAN_LIB = -L$(PETSC_DIR)/lib/$(PETSC_ARCH)
#PETSC_LIBS = -lpetsc -lpetscksp -lpetscdm -lpetscmat-lpetscvec -lpetscts

FLINKER = mpif90 -fPIC -g

linear-solver.x:   $(MYOBJS) $(MAIN.o) 
	-${FLINKER} -o $@ $(MAIN.o) $(MYOBJS) $(LIBFLAGS) ${PETSC_KSP_LIB}


debug:
	echo Fortran linker:  ${FLINKER}
	echo Fortran lib:  ${PETSC_FORTRAN_LIB}
