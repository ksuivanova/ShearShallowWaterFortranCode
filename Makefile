# "make" => generation du code de calcul
# "make traite" => Utilitaire de traitement des fichiers resultats 
# "make clean"  => effacement de tous les objets et modules
# "make cleanall"  => effacement de tous les objets, modules et executables compiles
# "make all " => le tout (sauf cleanall)

# Nom du compilateur sequentiel
COMPSER  = f90

# Nom du compilateur parallele
COMPILE =  mpif90
# mpif90 #ifort
# Options de Compilation
# OPT = -C
OPT =  -O3 

# Chemin des sources
PATHSRC=./

# Executable
SRC = main_grad
#Or2Parallele
EXE = $(SRC)
OBJS = $(SRC).o

# Code de Calcul
$(EXE): $(OBJS) 
	$(COMPILE) -o $(EXE) $(OBJS) $(OPT)
	make clean
$(OBJS): $(PATHSRC)/$(SRC).f90 Makefile
	$(COMPILE) $(OPT) -c $<
	
# Utilitaires

clean:
	rm -f *.mod
	rm -f *.*~
	rm -f *~
	rm -f *.o

cleanall:
	rm -f *.mod
	rm -f *.*~
	rm -f *~
	rm -f *.o
	rm -f $(EXE)

all:
	make

