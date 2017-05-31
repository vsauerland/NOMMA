# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -std=c++11

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -Wno-deprecated

# ---------------------------------------------------------------------

PROGRAMS = regEx

all:
	make $(PROGRAMS)

# ---------------------------------------------------------------------

clean :
	/bin/rm -rf $(PROGRAMS) *.o *.txt

# ---------------------------------------------------------------------

regEx: regEx.o regression.o knot.o
	$(CCC) $(CCOPT) regEx.o regression.o knot.o -o regEx
regEx.o: regEx.cpp 
	$(CCC) -c $(CCOPT) regEx.cpp -o regEx.o
regression.o: regression.cpp
	$(CCC) -c $(CCOPT) regression.cpp -o regression.o
knot.o: knot.cpp
	$(CCC) -c $(CCOPT) knot.cpp -o knot.o
