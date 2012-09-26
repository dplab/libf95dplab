#
# Suffix rules for fortran-90/95 files.
# 
# * includes "clean" and "force"
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
