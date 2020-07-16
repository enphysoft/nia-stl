# nia-stl
The neighbor-index-added stereolithography (nia-stl)
# Stereolithography format with the neighbor-index added 
## Objective
This package aims to calculate and include the primary and secondary nearest neighbors of each facet within a duplicated STL (stereolithography) file.  
See the image below, where the red boxes indicate the appended information of three nearest neighbors, i.e., edge-sharing neibhgor facets.
 
![nnbd stl file, data structure](https://github.com/enphysoft/nia-stl/image/cube-stl-nia-boxed.png)
 
 ## File list
- Source code
  - **append-stl-nnbs.f90** - the main f90 program 
- Library files
  - **fufs.f90** - includes frequently used format styles
  - **math.f90** - pure math functions including operator routines
  - **mathstl.f90** - basic math functions for STL format and data structures
- Make Utility 
  - **makefile** - the _makefile_ used to compile the source codes and generate an executable file, e.g., append-stl-nnbs.x
- Input STL files
  - **hinge.stl** - input stl-file used as the first argument, i.e., $ append-stl-nnbs.x hinge.stl, to generate hinge_nnbs.stl.
  - **STL_INPUT.stl** - default input stl-file if there is no argument, i.e., $ append-stl-nnbs.x, which will generate STL_INPUT_NNBS.stl.
- Output STL files: as explained above
  - **hinge_nnbs.stl**
  - **STL_INPUT_NNBS.stl**
- Output Analysis files
  - **NNBcheck.dat** - to compare two vertices, common to two edges of two contacting neighbors
  - **NNBfacet.dat** - to list three contacting facets (i.e., nearest neighbors) to a specific facet.
  - **NNBindex.dat** - to list the full nearest neighbors with vertices paired. 
  - **NNBlists.dat** - to list the three nearest neighbors to a specific fact (with its vertex coordinates).
  
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
$ append-stl-nnbs.x  
```
withtout an argument. In this case, append-stl-nnbs.x will read an input file of a default name, **STL_INPUT.stl**, and generate  **STL_INPUT_NNBS.stl**, which will include the three nearest neighbors per facet. 

###  To run with a specified STL input file
```bash
$ append-stl-nnbs.x hinge.stl 
```
will generate **hinge_nnbs.stl**. The output STL file has a postfix of "\_nnbs" before the dot and extension (.stl). 
