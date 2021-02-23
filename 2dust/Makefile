#=======================================================================
#
#  MASTER MAKEFILE FOR 2-Dust CODE
#
#=======================================================================
#-------------------  macro definitions  -------------------------------
#=======================================================================
### Default 
#=======================================================================
#FC	= f90
#FL	= f90	
#OPTS    = -g -c
#OPTS	= -O -c
#=======================================================================
### MIPSpro f90 compiler
#=======================================================================
##### Single Processor
#=======================================================================
#FC	= f90
#FL	= f90
#OPTS	= -g3 -c -OPT:Olimit=0 -fullwarn -check_bounds
#OPTS	= -c -O -fast
#OPTS	= -c -O3 -OPT:Olimit=0
#OPTS   = -c -OPT:Olimit=0
#=======================================================================
##### Multiple Processor
#=======================================================================
#FC	= f90
#FL	= f90 -apo
#FL	= f90 -mp -apo
#OPTS	= -g3 -c -mp -apokeep -OPT:Olimit=0 -fullwarn -check_bounds
#OPTS	= -c -mp -apokeep -O3
#OPTS	= -c -mp -apokeep -O3 -OPT:Olimit=0
#OPTS	= -c -mp -apokeep -OPT:Olimit=0
#OPTS   = -c -mp -apokeep -O3 -OPT:Olimit=0
#OPTS	= -c -mp -apokeep -O3 -OPT:Olimit=0:IEEE_arithmetic=1:IEEE_NaN_inf=ON
#OPTS	= -c -apokeep -O3 -OPT:Olimit=0:IEEE_arithmetic=1
#=======================================================================
### Sun Workstation @ STScI
#=======================================================================
#FC	= /opt/SUNWspro-7.0/bin/f95 -L/opt/SUNspro-7.0/lib
#FL	= /opt/SUNWspro-7.0/bin/f95 -L/opt/SUNspro-7.0/lib
#OPTS	= -O -c
#=======================================================================
### Intel Fortran Compiler
#=======================================================================
#FC	= ifc
#FL	= ifc
#OPTS    = -g -c -e90
#OPTS	= -O3 -mp -c -e90
#=======================================================================
### Intel Fortran Compiler for Mac OSX
#=======================================================================
FC	= ifort
FL	= ifort
#OPTS	= -g -c
OPTS	= -O3 -mp -c
#=======================================================================
### Portland Group f90 compiler
#=======================================================================
#FC	= pgf90 -Bstatic
#FL	= pgf90 -Bstatic
#OPTS	= -O -c
#=======================================================================
### IBM XL Fortran for MacOSX
#=======================================================================
#FC	= /opt/ibmcmp/xlf/8.1/bin/f90 -L/opt/ibmcmp/xlf/8.1/lib
#FL	= /opt/ibmcmp/xlf/8.1/bin/f90 -L/opt/ibmcmp/xlf/8.1/lib
#OPTS	= -g -c 
#OPTS	= -O -c -qmaxmem=-1
#=======================================================================
### Absoft Fortran for MacOSX
#=======================================================================
#FC	= f90
#FL	= f90
#OPTS	= -z0 -en -Rb -Rc -Rs -Rp -g -c 
#OPTS	= -z0 -en -g -c 
#OPTS	= -O -c
#=======================================================================
#-------------------  object files  ------------------------------------
#=======================================================================
OBJ 	=	prcsns.o const.o defvars.o getrho.o phase.o dfunc.o\
		gauleg.o mie.o pah.o dustprep.o\
		locate.o mipcub.o gridgen.o mindens.o shmass.o bbtabgen.o\
		bbtabgen2.o temp.o extrap.o extrap2.o geom.o readsed.o\
		twodust.o
#
OBJd 	=	prcsnd.o const.o defvars.o getrho.o phase.o dfunc.o\
		gauleg.o mie.o pah.o dustprep.o\
		locate.o mipcub.o \
		gridgen.o mindens.o shmass.o bbtabgen.o bbtabgen2.o\
		temp.o extrap.o extrap2.o geom.o readsed.o twodust.o
#
OBJ1 	=	prcsns.o const.o defvars.o getrho.o phase.o dfunc.o\
		gauleg.o mie.o pah.o dustprep.o \
		locate.o mipcub.o \
		gridgen.o mindens.o shmass.o bbtabgen.o bbtabgen2.o\
		temp.o extrap.o extrap2.o geom.o agnprep.o twodust.o
#
OBJ1d 	=	prcsnd.o const.o defvars.o getrho.o phase.o dfunc.o\
		gauleg.o mie.o pah.o dustprep.o\
		locate.o mipcub.o \
		gridgen.o mindens.o shmass.o bbtabgen.o bbtabgen2.o    \
		temp.o extrap.o extrap2.o geom.o agnprep.o twodust.o
#
OBJ2	=	prcsnd.o const.o defvars.o dfunc.o getrho.o phase.o    \
		gauleg.o mie.o pah.o locate.o \
		gridgen.o mindens.o mapspec.o
#
OBJ3    =	prcsnd.o const.o fits_io.o mapgrid.o
#
.SUFFIXES: .f90 .o
#
.f90.o: 
	${FC} ${OPTS} $*.f90
#
2dust:	${OBJ}
	${FL} -o 2dust.exe ${OBJ}
#
2dustd:	${OBJd}
	${FL} -o 2dustd.exe ${OBJd}
#
agn: 	${OBJ1}
	${FL} -o 2dustagn.exe ${OBJ1}
#
agnd: 	${OBJ1d}
	${FL} -o 2dustagnd.exe ${OBJ1d}
#
mapspec:${OBJ2}
	${FL} -o mapspec.exe ${OBJ2}
#
mapgrid:${OBJ3}
	${FL} -o mapgrid.exe ${OBJ3}
#
clean:	
	\rm -f *.o *.M *~ *.list *.mod *.anl *.m *.stb
# 
2d:	clean 2dust
#
2dd:	clean 2dustd
#
ms: 	clean mapspec
#
mg:	clean mapgrid
#
2da: 	clean agn 
#
2dda: 	clean agnd
#
para:
	${FC} -c -O3 -OPT:Olimit=0 prcsns.f90
	${FC} -c -O3 -OPT:Olimit=0 const.f90
	${FC} -c -O3 -OPT:Olimit=0 defvars.f90
	${FC} -c -O3 -OPT:Olimit=0 getrho.f90
	${FC} -c -O3 -OPT:Olimit=0 phase.f90
	${FC} -c -O3 -OPT:Olimit=0 dfunc.f90
	${FC} -c -O3 -OPT:Olimit=0 gauleg.f90
	${FC} -c -O3 -OPT:Olimit=0 mie.f90
	${FC} -c -O3 -OPT:Olimit=0 dustprep.f90
	${FC} -c -O3 -OPT:Olimit=0 locate.f90
	${FC} -c -O3 -OPT:Olimit=0 mipsub.f90
	${FC} -c -O3 -OPT:Olimit=0 gridgen.f90
	${FC} -c -O3 -OPT:Olimit=0 mindens.f90
	${FC} -c -O3 -OPT:Olimit=0 shmass.f90
	${FC} -c -O3 -OPT:Olimit=0 bbtabgen.f90
	${FC} -c -O3 -OPT:Olimit=0 bbtabgen2.f90
	${FC} -c -O3 -OPT:Olimit=0 temp.f90
	${FC} -c -O3 -OPT:Olimit=0 extrap.f90
	${FC} -c -O3 -OPT:Olimit=0 extrap2.f90
	${FC} -c -O3 -OPT:Olimit=0 geom.f90
	${FC} -c -O3 -OPT:Olimit=0 agnprep.f90
	${FC} -c -O3 -OPT:Olimit=0 -apo -mp twodust.f90
	${FL} -o 2dustpara.exe ${OBJ}
#

