FC=gfortran
# FC=ifort

srcroot=main-nia-stl
srcfile=$(srcroot)$(version).f90
exefile=$(srcroot).x
OBJ=fufs.o math.o   mathstl.o

%.o: %.f90
	$(FC) $(FCOPT) -c $?
%.o: %.f
	$(FC) $(FCOPT) -c $?

cmp: $(OBJ)
	@echo
	@echo "#! Compiling $(srcfile) using $(FC) ... "
	@echo
	$(FC) $(srcfile)  $(OBJ) -o   $(exefile)
	@echo
	@echo "#! Compiling $(srcfile) done! \n"

run:
	mkdir -p data 
	./$(exefile)  stl/ascii-cube.stl

clean:
	rm -rf *.o *.mod *.x *.exc a.out NNB*.dat NNB*.csv 

all: clean cmp run 

