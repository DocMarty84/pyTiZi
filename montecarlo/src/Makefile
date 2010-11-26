##
###################################################################
##                                                               ##
##         Makefile for Building the Monte Carlo Package         ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the executables
##   2. make rename           Move the executables to BINDIR
##   3. make listing          Concatenate source to MC.txt
##   4. make clean            Delete objects and executables
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##
##  MCDIR        MC Distribution Directory
##  BINDIR       Hard Copies of MC Executables
##

MCDIR = "."
BINDIR = $(MCDIR)/../bin

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible
##  Oper Sys: Ubuntu 10.04
##  Compiler: GNU g++
##

CPP = /usr/bin/g++
LIBS = -lm
OMPLIBS = -fopenmp
CPPFLAGS = -Wall -c
OPTFLAGS = -O3
LIBFLAGS = -crusv
LINKFLAGS =

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible
##  Oper Sys: CentOS 5.5
##  Compiler: Portland Group (PGI) Fortran
##  Notice:   Load the appropriate module if needed
##

#CPP = pgcpp
#LIBS =
#OMPLIBS = -mp
#CPPFLAGS = -c
#OPTFLAGS = -fast -O3
#LIBFLAGS = -crusv
#LINKFLAGS =


#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = buildgrid.o \
       clear.o \
       coordinates.o \
       deltae.o \
       deltav.o \
       dispatch.o \
       distance.o \
       electricfield.o \
       mathplus.o \
       mcbkl.o \
       mcfrm.o \
       MC_BKL.o \
       MC_BKL_layer.o \
       printsummary.o \
       read.o \
       transferrate.o

EXEFILES = MC_BKL.x \
           MC_BKL_layer.x

.cpp.o:
	${CPP} ${CPPFLAGS} $*.cpp ${OMPLIBS} ${OPTFLAGS}

%.x: %.o libmonte.a
	${CPP} ${LINKFLAGS} -o $@ $^ ${LIBS} ${OMPLIBS}

all:	${EXEFILES}

clean:
	rm -f *.o *.a *.x

listing:
	cat *.i *.cpp > MC.txt

rename:
	mkdir -p $(BINDIR)
	mv  MC_BKL.x          $(BINDIR)/MC_BKL
	mv  MC_BKL_layer.x    $(BINDIR)/MC_BKL_layer

libmonte.a: ${OBJS} 
	ar ${LIBFLAGS} libmonte.a \
	buildgrid.o \
	clear.o \
	coordinates.o \
	deltae.o \
	deltav.o \
	dispatch.o \
	distance.o \
	electricfield.o \
	mathplus.o \
	mcbkl.o \
	mcfrm.o \
	printsummary.o \
	read.o \
	transferrate.o

