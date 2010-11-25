##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Molecular Modeling Package  ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the TINKER executables
##   2. make rename           Move the executables to BINDIR
##   3. make remove_links     Create soft links in LINKDIR
##   4. make create_links     Remove soft links from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executables
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, Jan. 1998,
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##
##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Hard Copies of TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables
##

TINKERDIR = /Users/ponder/tinker
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible
##  Oper Sys: Fedora Core
##  Compiler: GNU gfortran
##

CPP = /usr/bin/g++
LIBS = -lm -fopenmp
CPPFLAGS = -c
OPTFLAGS = -O3
LIBFLAGS = -crusv
LINKFLAGS =

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible
##  Oper Sys: Fedora Core
##  Compiler: Portland Group (PGI) Fortran
##

#CPP = /opt/pgi/linux86/bin/pgcpp
#LIBS =
#CPPFLAGS = -c
#OPTFLAGS = -fast -O3
#LIBFLAGS = -crusv
#LINKFLAGS =


#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = clear.o \
       coordinates.o \
       mathplus.o \
       MC_BKL.o \
       MC_BKL_layer.o

EXEFILES = MC_BKL.x \
           MC_BKL_layer.x

.cpp.o:
	${CPP} ${CPPFLAGS} ${OPTFLAGS} $*.cpp

%.x: %.o libmonte.a
	${CPP} ${LINKFLAGS} -o $@ $^ ${LIBS}

all:	${EXEFILES}

clean:
	rm -f *.o *.a *.x

listing:
	cat *.i *.cpp > MC.txt

rename:
	mv  MC_BKL.x          $(BINDIR)/MC_BKL
	mv  MC_BKL_layer.x    $(BINDIR)/MC_BKL_layer

libmonte.a: ${OBJS} 
	ar ${LIBFLAGS} libmonte.a \
        clear.o \
        coordinates.o \
        mathplus.o

