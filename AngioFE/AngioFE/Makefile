
FECORE_INC = $(HOME)/FEBio2
FECORE_LIB = $(HOME)/FEBio2/build/lib

ANGIO_INC = ../../AngioLib
ANGIO_LIB = ../../AngioLib

INC = -I$(FECORE_INC) -I$(ANGIO_INC)

all:
	icpc -fPIC -c *.cpp $(INC) -g
	icpc -shared -Wl,-soname,libangiofe.so -o libangiofe.so *.o -L$(FECORE_LIB) -lfecore_lnx64 -L$(ANGIO_LIB) -langio
