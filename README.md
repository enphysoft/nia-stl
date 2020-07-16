# nia-stl
The neighbor-index-added stereolithography (nia-stl)
# Stereolithography format with the neighbor-index added 
## Objective
This package aims to calculate and include the primary and secondary nearest neighbors of each facet within a duplicated STL (stereolithography) file.  
See the image below, where the red boxes indicate the appended information of three nearest neighbors, i.e., edge-sharing neibhgor facets.
![nnbd stl file, data structure](https://github.com/enphysoft/nia-stl/blob/master/image/cube-stl-org.png)
![nnbd stl file, data structure](https://github.com/enphysoft/nia-stl/blob/master/image/cube-stl-nia-boxed.png) 
 
 ## File list
- Source code
  - **main-nia-stl.f90** - the main f90 program 
- Library files
  - **fufs.f90** - includes frequently used format styles
  - **math.f90** - pure math functions including operator routines
  - **mathstl.f90** - basic math functions for STL format and data structures
- Make Utility 
  - **Makefile** - the _makefile_ used to compile the source codes and generate an executable file, e.g., append-stl-nnbs.x
- Input STL files
  - **stl/ascii-hinge.stl** - input stl-file used as the first argument, i.e., $ append-stl-nnbs.x hinge.stl, to generate hinge_nnbs.stl.
  - **stl/ascii-cube.stl** - default input stl-file if there is no argument, i.e., $ append-stl-nnbs.x, which will generate STL_INPUT_NNBS.stl.
- Output STL files: as explained above
  - **hinge_nnbs.stl**
  - **stl/ascii-cube_nia.stl**
- Output Analysis files
  - **data/NNBcheck.dat** - to compare two vertices, common to two edges of two contacting neighbors
  - **data/NNBfacet.dat** - to list three contacting facets (i.e., nearest neighbors) to a specific facet.
  - **data/NNBverix.csv** - to list the full nearest neighbors with vertices paired.  
  
## How to compile and run
### To clean 
```bash
$ make clean 
``` 
### To compile the source code *append-stl-nnbs.f90*
```bash
$ make 
```
### To run without an input file
```bash
$ make run 
```
equivalent to  
```bash
$ ./main-nia-stl.x stl/ascii-cube.stl 
```
In this case, main-nia-stl.x will open and read an input file **stl/ascii-cube.stl**, and generate  **ascii-cube_nia.stl**, which will include the all the nearest neighbor indices per facet. 


