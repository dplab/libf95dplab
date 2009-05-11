# $Id: makefile,v 1.45 2001/08/22 18:00:07 balay Exp $ 

FC = ifort

CFLAGS	         = -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE}
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =

default:  linear-solver.x

include ${PETSC_DIR}/bmake/common/base

#
MAIN.o = linear-solver.o
MYOBJS = LinearSolverPetscModule.o
MYLIB  = libf90.a 
#
#
#  My f90 options.
#
OTHER_OPTS = -error-limit 5 -check all -traceback -check noarg_temp_created
#
FFLAGS           = -g ${OTHER_OPTS}
LOCDIR		 = .

#FC_FLAGS =  -I. -fPIC -names lowercase
#FLINKER= mpif90 -I. -fPIC 
#PETSC_FORTRAN_LIB = -L$(PETSC_DIR)/lib/$(PETSC_ARCH)
#PETSC_LIBS = -lpetsc -lpetscksp -lpetscdm -lpetscmat-lpetscvec -lpetscts

FLINKER = /opt/openmpi-1.2.7/bin/mpif90

linear-solver.x:   $(MYOBJS) $(MAIN.o) 
	-${FLINKER} -o $@ $(MAIN.o) $(MYOBJS) $(LIBFLAGS) ${PETSC_KSP_LIB}


debug:
	echo Fortran linker:  ${FLINKER}
	echo Fortran lib:  ${PETSC_FORTRAN_LIB}
