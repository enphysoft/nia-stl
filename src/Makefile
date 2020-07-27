
# FC=ifort
# CXX=icpc

FC=gfortran
CXX=g++ 

srcrootF90=main-nia-stl-ask
srcfileF90=$(srcrootF90).f90
ojbfileF90=$(srcrootF90).o
exefileF90=$(srcrootF90).x
OBJF90=fufs.o math.o mathstl.o $(ojbfileF90) 



srcrootCPP=main-nia-stl-ejj
srcfileCPP=$(srcrootCPP).cpp
ojbfileCPP=$(srcrootCPP).o
exefileCPP=$(srcrootCPP).x
OBJCPP=$(ojbfileCPP) 


%.o: %.f90
	$(FC) $(FCOPT) -c $?
%.o: %.f
	$(FC) $(FCOPT) -c $?
%.o: %.cpp
	$(CXX) $(CPPOPT) -c $?
%.o: %.C 
	$(CXX) $(CPPOPT) -c $?

ALL: allclean allcmp allrun allcompare
allcmp:  flnk clnk
allrun:  frun crun
allcompare: fcompare ccompare
allclean: fclean cclean stlclean


flnk: fcmp 
	@echo
	@echo "#! Compiling $(srcfileF90) using $(FC) ... "
	@echo
	$(FC) $(OBJF90) -o   $(exefileF90)
	@echo
	@echo "#! Compiling $(srcfileF90) done!"
	@echo "#! To run:\n   make frun \n"

clnk: ccmp 
	@echo
	@echo "#! Compiling $(srcfileCPP) using $(CXX) ... "
	@echo
	$(CXX) $(OBJCPP) -o   $(exefileCPP)
	@echo
	@echo "#! Compiling $(srcfileCPP) done!"
	@echo "#! To run:\n   make crun \n"


fcmp:  $(OBJF90)
ccmp:  $(OBJCPP)

frun:
	./$(exefileF90)  ascii-cube.stl
	./$(exefileF90)
	@echo "\n Also try: $ ./$(exefileF90)  mystl_nnbs.stl \n"

crun:
	./$(exefileCPP)  ascii-cube.stl
	./$(exefileCPP)  
	@echo "\n equivalent to : $ ./$(exefileCPP)  STL_INPUT.stl \n"

fclean:
	rm -f  $(OBJF90) $(exefileF90) *.mod NNB*.dat NNB*.csv
cclean:
	rm -f  $(OBJCPP) $(exefileCPP) 

stlclean:
	rm -f *nia*.stl



# all: clean cmp run 

fcompare:
	compare_nnb2.sh STL_INPUT_cnia1.stl STL_INPUT_fnia1.stl
ccompare:
	compare_nnb2.sh STL_INPUT_cnia0.stl STL_INPUT_fnia0.stl

