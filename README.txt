
Libf95dplab
===========


This library has fortran utilities of general use.

USAGE

To access these routines, use the LibF95 module in
your programs, like this:

USE LibF95


BUILDING

To build, go to the Source directory.  There are sample
makefiles there for different compilers and platforms. 
Use one of those if possible; otherwise modify one to meet
your needs.  Build using the "-f" option to make. For example:

make -f Make.ifort-linux

will build for the ifort (intel) compiler on a linux system.

The build will produce the libary, "libf95dplab.a", and will
place module definitions in the "Modules" subdirectory of Source.
You will need to remember these locations for linking.


