#
# Defines sources, objects and library.
#
LIBBASE = f95dplab
LIBRARY = lib$(LIBBASE).a
#
SOURCES=\
IntrinsicTypesModule.f90\
ConstantsModule.f90\
StringsModule.f90\
FilesModule.f90\
TimerModule.f90\
CommandModule.f90\
OptimModule.f90\
LibF95.f90
#
OBJECTS = $(SOURCES:.f90=.o)
