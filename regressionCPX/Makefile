SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#-----------------------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where EIGEN, CPLEX and CONCERT are installed.
#
#-----------------------------------------------------------------------------

CPLEXDIR   = /cluster/Software/Dcplex/cplex
CONCERTDIR = /cluster/Software/Dcplex/concert

# ----------------------------------------------------------------------------
# Compiler selection 
# ----------------------------------------------------------------------------

CCC = g++ -std=c++11

# ----------------------------------------------------------------------------
# Compiler options 
# ----------------------------------------------------------------------------

CCOPT = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -Wno-deprecated

# ----------------------------------------------------------------------------
# Link options and libraries
# ----------------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

# ----------------------------------------------------------------------------

PROGRAMS = regEx 

all:
	make $(PROGRAMS)
clean :
	/bin/rm -rf $(PROGRAMS) *.o *.txt

# ----------------------------------------------------------------------------

regEx: regEx.o regression.o knot.o
	$(CCC) $(CCOPT) regEx.o regression.o knot.o -o regEx $(CCLNFLAGS)
regEx.o: regEx.cpp 
	$(CCC) -c $(CCFLAGS) regEx.cpp -o regEx.o $(CCLNFLAGS)
regression.o: regression.cpp
	$(CCC) -c $(CCFLAGS) regression.cpp -o regression.o $(CCLNFLAGS)
knot.o: knot.cpp
	$(CCC) -c $(CCFLAGS) knot.cpp -o knot.o $(CCLNFLAGS)
