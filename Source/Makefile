#    -*-mode: makefile; abbrev-mode: t -*-
#
#  f90dplab/modules/BUILD	
#
#  Makefile for building f90 modules library. [serial version]
#
LIBBASE = f95dplab
#
SOURCES=\
IntrinsicTypesModule.f95\
ConstantsModule.f95\
StringsModule.f95\
FilesModule.f95\
TimerModule.f95
#
F90 = ifort
F95 = ifort
F95FLAGS = -O  -module Modules
#
#--------------------*--------------------------------------------------
#
OBJECTS = $(SOURCES:.f95=.o)

LIBRARY = lib$(LIBBASE).a
#
default:  $(LIBRARY)
#
$(LIBRARY): Modules $(LIBRARY)($(OBJECTS))
	ranlib $@
Modules: 
	mkdir Modules
#
#  Dependencies:
#
$(LIBRARY)(ConstantsModule.o):    $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(FilesModule.o):        $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(StringsModule.o):      $(LIBRARY)(IntrinsicTypesModule.o)
#
MODMODULES = \
constantsmodule.mod\
crystaltypemodule.mod\
filesmodule.mod\
intrinsictypesmodule.mod\
tensor3dmodule.mod

DIRT=$(LIBRARY) Modules/*.mod
#
#  Suffix rules for fortran-90/95 files.
#
.SUFFIXES: .f90 .F90 .ff .INC .inc .f95 .F95
#
.f90.o:
	$(F90) -c $(F90FLAGS) $<
.f90:
	$(F90) -o $@ $(F90FLAGS) $<
.f95.o:
	$(F95) -c $(F95FLAGS) $<
.f95:
	$(F95) -o $@ $(F95FLAGS) $<
.F90.f90:
	/bin/rm -f $@
	$(CPP) $(CPPFLAGS) $< > $@
	chmod 444 $@
.INC.inc:
	/bin/rm -f $@
	$(CPP) $(CPPFLAGS) $< > $@
	chmod 444 $@
.ff.o:
	$(F90) -c $(FFLAGS) $<

#
#  Rules for `clean' target.
#
#  TARGETS:  clean
#    NEEDS:  DIRT
#  DEFINES:  REMOVE
#     USES:  force
#
REMOVE=/bin/rm -f
#
clean: force
	$(REMOVE) $(DIRT)
#
#
#  Rule to force execution of another rule.
#
#  TARGETS:  force
#
force:
