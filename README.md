# nia-stl
The neighbor-index-added stereolithography (nia-stl).
A previous version can be found at https://github.com/enphysoft/append-stl-nnbs
# Stereolithography format with the neighbor-index added 
## Objective
This package aims to calculate and include the primary and secondary nearest neighbors of each facet within a duplicated STL (stereolithography) file. The left and right images below are the original and neighbor-index-added STL files. On the seventh line of the nia-stl file, the sequential index of facet, its three primary nearest neighbor indicies, and multiple secondary nearest neighbor indices are shown.  

**from**
![nnbd stl file, data structure](https://github.com/enphysoft/nia-stl/blob/master/image/cube-stl-org3.png) **to**
![nnbd stl file, data structure](https://github.com/enphysoft/nia-stl/blob/master/image/cube-stl-nia-boxed3.png)
 
 ## File list
- Source code
  - **main-nia-stl.f90** - the main f90 program 
- Library files
  - **fufs.f90** - includes frequently used format styles
  - **math.f90** - pure math functions including operator routines
  - **mathstl.f90** - basic math functions for STL format and data structures
- Make Utility 
  - **Makefile** - the _Makefile_ used to compile the source codes and generate an executable file, e.g., **main-nia-stl.x**
- Sample STL files in ./stl/ directory
  - **stl/ascii-hinge.stl**  
  - **stl/ascii-cube.stl** 
- Output STL files: 
  - Output files will have **nia** as postfix of the input STL file name.
- Output Analysis files
  - **data/NNBcheck.dat** - to compare two vertices, common to two edges of two contacting neighbors
  - **data/NNBfacet.dat** - to list three contacting facets (i.e., nearest neighbors) to a specific facet.
  - **data/NNBvertx.csv** - to list the full nearest neighbors: facet id, the primary nearest neighbors, and the secondary nearest neighbors.  
  
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


