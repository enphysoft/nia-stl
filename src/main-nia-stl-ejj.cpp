/*
Copyright © 2019-present, by Albert S. Kim
Authors: Eun Ji Jun and Albert S. Kim ( http://albertsk.org/ )

Version: 1.0
Package-Version: 20190630.0627
Created: 06/30/2019
Keywords: STL, nearest neighbor
Description:
This program is, for a selected facet, to generate a list of
three nearest facets, sharing an edge, connected by two shared vertices.
For example, if a facet 45 and 72 are nearest neighbors, and
vertices 1 and 2 of facet 45 and vertices 2 and 3 of facet 72 are
paired, respectively, the outcome looks like
  45 72
  1  2
  2  3
  3  0
where the last line indicates that vertex 3 of facet 45 is not
paired with any of vertices of facet 72.

entree = mystl.stl
output = mystlNNB.dat, facetNNB.dat, checkNNB.dat
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include<algorithm>
#include<iterator>
#include<math.h>
#include "main-nia-stl-ejj.h"
#define TOL 1.E-8

int main(int argc, char **argv){
  if (argc == 1) {
    entree_fileSTL = defaultfileSTL;
    cout << "Input  STL  file:" << " " << entree_fileSTL <<"\n";
    cout << "Output STL0 file:" << " " << output_fileSTL0 <<"\n";
    cout << "Output STL1 file:" << " " << output_fileSTL1 <<"\n";
    // cout << "entree_fileSTL  :" << " " << entree_fileSTL <<"\n";
    // cout << "output_fileSTL  :" << " " << output_fileSTL <<"\n";
    // cout << "review_fileNNB:" << " " << review_fileNNB <<"\n";
    // cout << "facets_fileNNB:" << " " << facets_fileNNB <<"\n";
    // cout << "runtst_NNBlist:" << " " << runtst_NNBlist <<"\n";
    // cout << "output_fileNNB:" << " " << output_fileNNB <<"\n";
  }
  else if (argc == 2){
    entree_fileSTL = argv[argc-1];
    arg1_root = entree_fileSTL.substr(0, entree_fileSTL.find(dot));
    arg1_extn = entree_fileSTL.substr(entree_fileSTL.find(dot) + 1);
    output_fileSTL  = arg1_root + postfix_nnbed + dot + arg1_extn;
    output_fileSTL0 = arg1_root + postfix_nnbed0 + dot + arg1_extn;
    output_fileSTL1 = arg1_root + postfix_nnbed1 + dot + arg1_extn;
    cout << "Input  STL file - given  :" << " " << entree_fileSTL <<"\n";
    cout << "Output STL file - C style:" << " " << output_fileSTL0 <<"\n";
    cout << "output STL file - F style:" << " " << output_fileSTL1 <<"\n";
    // cout << "entree_fileSTL:" << " " << entree_fileSTL <<"\n";
    // cout << "review_fileNNB:" << " " << review_fileNNB <<"\n";
    // cout << "facets_fileNNB:" << " " << facets_fileNNB <<"\n";
    // cout << "runtst_NNBlist:" << " " << runtst_NNBlist <<"\n";
    // cout << "output_" <<"\n";
  }
  else{
    cout << "!!!Correct usuage:" << " " << *argv <<"\n";
  }

  ifstream inFile(entree_fileSTL);
  numline	= count(istreambuf_iterator<char>(inFile),istreambuf_iterator<char>(), '\n');
  numtri	= (numline-2)/7;
  N_fct		= numtri;
  N_fct_o10	= N_fct / 10;
  cout << "number of lines in input file = " << " " << numline <<"\n";
  cout << "number of triangles (facets) = " << " " << numtri <<"\n";
  myFacet	= new stl_tri	[numtri];
  loadedFacet	= new stl_tri	[numtri];
  myNNpal	= new stl_nnb	[numtri];
  myNNpal	= new stl_nnb	[numtri];
  loadedNNbor	= new stl_nnb	[numtri];
  numFacets_i	= new int	[numtri];
  numVertices_i = new int	[numtri];
  nnbedFacet	= new stl_tri_nnb[numtri];

  read_ascii_stl (entree_fileSTL,numtri,myFacet,title_stl);
  cout << "numtri = " << numtri << "\n";
  cout << "reading STL file done" << "\n";

  for (i_fct = 0;	i_fct < N_fct;	i_fct++) {
    for (int k = 0;	k < 3;	k++	) {
      myNNpal[i_fct].mvtx[0][k] = myFacet[i_fct].tvtx.vertexA[k];
      myNNpal[i_fct].mvtx[1][k] = myFacet[i_fct].tvtx.vertexB[k];
      myNNpal[i_fct].mvtx[2][k] = myFacet[i_fct].tvtx.vertexC[k];
    }
  }

  for (i_fct = 0; i_fct < N_fct-1 ; i_fct++){
    numFacets_i   [i_fct] = 0;
    numVertices_i [i_fct] = 0;
  }

  for (i_fct = 0; i_fct < N_fct ; i_fct++){
    for (k_fct = i_fct+1 ; k_fct < N_fct; k_fct++){
      numPair_i_l	= 0;

      for (j_vtx = 0; j_vtx < 3; j_vtx++) {
        for (l_vtx = 0; l_vtx < 3; l_vtx++) {
          for (int k = 0; k < 3; k++) {
            vecA[k] = myNNpal[i_fct].mvtx[l_vtx][k];
            vecB[k] = myNNpal[k_fct].mvtx[j_vtx][k];
          }

          if ( abs (vecA[0] - vecB[0]) < TOL && abs (vecA[1] - vecB[1]) < TOL && abs (vecA[2] - vecB[2]) < TOL ){
            numPair_i_l			= numPair_i_l + 1 ;
          }
        }
      }

      if (numPair_i_l == 1){
        myNNpal	[i_fct].vidPair[numVertices_i[i_fct]]	= k_fct;
        numVertices_i [i_fct]					= numVertices_i[i_fct] + 1;
        myNNpal	[i_fct].vidPairMax			= numVertices_i[i_fct];
        myNNpal	[k_fct].vidPair[numVertices_i[k_fct]]	= i_fct;
        numVertices_i [k_fct]					= numVertices_i[k_fct] + 1;
        myNNpal	[k_fct].vidPairMax			= numVertices_i[k_fct];
      }
      else if (numPair_i_l == 2) {
        myNNpal[i_fct].nidPair[numFacets_i[i_fct]]	= k_fct;
        numFacets_i [i_fct]				= numFacets_i[i_fct] + 1;
        myNNpal[k_fct].nidPair[numFacets_i[k_fct]]	= i_fct;
        numFacets_i [k_fct]				= numFacets_i[k_fct] + 1;
      }
    }
  }

  // C-style output: index from 0 to N-1 for N components
  write_ascii_stl_w_nnb_all (output_fileSTL0,myFacet,myNNpal,numtri,title_stl,0);
  // Algebraic and f90-style output: index from 1 to N for N components
  write_ascii_stl_w_nnb_all (output_fileSTL1,myFacet,myNNpal,numtri,title_stl,1);

  read_ascii_stl_w_nnb_all  (output_fileSTL,loadedFacet,loadedNNbor,numtri,title_stl);

  delete myFacet;
  delete myNNpal;
  delete numFacets_i;
  delete numVertices_i;
  delete loadedFacet;
  delete loadedNNbor;
  delete nnbedFacet;

  return 0;
}

/**************************************************************/
/****** read_ascii_stl_file - to read an ascii stl fil ********/
/**************************************************************/
void read_ascii_stl(string inputfile, int numrec, stl_tri *hx, string title){
  cout << "=== The title of the stl file is (from subroutine read_ascii_stl): " <<"\n";
  file.open(entree_fileSTL.c_str());
  // read first line of stl file
  getline (file, line);
  title_stl = line;
  file.close();
  file.open(entree_fileSTL.c_str());
  string sLine;
  getline(file, sLine);
  for (int itri = 0; itri < numrec; itri++){
    myFacet[itri].tid = itri;
    for (int index = 0; index < 21; index ++){
      file >> word;
      if	 (index ==  2) myFacet[itri].tnvec[0] = atof(word.c_str());
      else if (index ==  3) myFacet[itri].tnvec[1] = atof(word.c_str());
      else if (index ==  4) myFacet[itri].tnvec[2] = atof(word.c_str());
      //
      else if (index ==  8) myFacet[itri].tvtx.vertexA[0] = atof(word.c_str());
      else if (index ==  9) myFacet[itri].tvtx.vertexA[1] = atof(word.c_str());
      else if (index == 10) myFacet[itri].tvtx.vertexA[2] = atof(word.c_str());
      //
      else if (index == 12) myFacet[itri].tvtx.vertexB[0] = atof(word.c_str());
      else if (index == 13) myFacet[itri].tvtx.vertexB[1] = atof(word.c_str());
      else if (index == 14) myFacet[itri].tvtx.vertexB[2] = atof(word.c_str());
      //
      else if (index == 16) myFacet[itri].tvtx.vertexC[0] = atof(word.c_str());
      else if (index == 17) myFacet[itri].tvtx.vertexC[1] = atof(word.c_str());
      else if (index == 18) myFacet[itri].tvtx.vertexC[2] = atof(word.c_str());
    }
  }
  file.close();
}

/***********************************************************************/
/******  read_ascii_stl_w_nnb_all - to write an ascii stl file *********/
/***********************************************************************/
void read_ascii_stl_w_nnb_all(string entreefile, stl_tri *hx1, stl_nnb *hxnnb, int numtri, string title){
  file.open(entreefile.c_str());
  // read first line of stl file
  getline (file, line);
  title_stl = line;
  file.close();

  file.open(entreefile.c_str());
  file >> word;
  for (int ifct = 0; ifct < numtri; ifct++){
    loadedFacet[ifct].tid = ifct;
    for (int index = 0; index < 24; index ++){
      file >> word;
      if (index == 2) loadedFacet[ifct].tnvec[0] = atof(word.c_str());
      else if (index == 3) loadedFacet[ifct].tnvec[1] = atof(word.c_str());
      else if (index == 4) loadedFacet[ifct].tnvec[2] = atof(word.c_str());
      else if (index == 8) loadedFacet[ifct].tvtx.vertexA[0] = atof(word.c_str());
      else if (index == 9) loadedFacet[ifct].tvtx.vertexA[1] = atof(word.c_str());
      else if (index == 10) loadedFacet[ifct].tvtx.vertexA[2] = atof(word.c_str());
      else if (index == 12) loadedFacet[ifct].tvtx.vertexB[0] = atof(word.c_str());
      else if (index == 13) loadedFacet[ifct].tvtx.vertexB[1] = atof(word.c_str());
      else if (index == 14) loadedFacet[ifct].tvtx.vertexB[2] = atof(word.c_str());
      else if (index == 16) loadedFacet[ifct].tvtx.vertexC[0] = atof(word.c_str());
      else if (index == 17) loadedFacet[ifct].tvtx.vertexC[1] = atof(word.c_str());
      else if (index == 18) loadedFacet[ifct].tvtx.vertexC[2] = atof(word.c_str());
      else if (index == 21) loadedNNbor[ifct].nidPair[0] = atof(word.c_str());
      else if (index == 22) loadedNNbor[ifct].nidPair[1] = atof(word.c_str());
      else if (index == 23) loadedNNbor[ifct].nidPair[2] = atof(word.c_str());
     }
   }
  file.close();
}

/***********************************************************************/
/******  write_ascii_stl_w_nnb_all - to write an ascii stl file ********/
/***********************************************************************/
void write_ascii_stl_w_nnb_all(string outputfile, stl_tri *hx, stl_nnb *hxnnb, int numtri, string title, int idx_init){
  ofstream outFile(outputfile);
  outFile << "solid" << " " << title <<"\n";
  // cout << title << "\n";
  for (int ifct = 0; ifct < numtri; ifct++){
    outFile << "facet normal" << " " ;
    outFile << myFacet[ifct].tnvec[0] << " ";
    outFile << myFacet[ifct].tnvec[1] << " ";
    outFile << myFacet[ifct].tnvec[2] << "\n";
    outFile << "  outer loop" << "\n";
    outFile << "    vertex" << " " ;
    outFile << myFacet[ifct].tvtx.vertexA[0] << " " ;
    outFile << myFacet[ifct].tvtx.vertexA[1] << " " ;
    outFile << myFacet[ifct].tvtx.vertexA[2] << "\n";
    outFile << "    vertex" << " " ;
    outFile << myFacet[ifct].tvtx.vertexB[0] << " " ;
    outFile << myFacet[ifct].tvtx.vertexB[1] << " " ;
    outFile << myFacet[ifct].tvtx.vertexB[2] << "\n";
    outFile << "    vertex" << " " ;
    outFile << myFacet[ifct].tvtx.vertexC[0] << " " ;
    outFile << myFacet[ifct].tvtx.vertexC[1] << " " ;
    outFile << myFacet[ifct].tvtx.vertexC[2] << "\n";
    outFile << "  endloop" << "\n";
    outFile << "endfacet" << " " ;
    outFile << ifct			+ idx_init << "," ;
    outFile << myNNpal[ifct].nidPair[0] + idx_init << "," ;
    outFile << myNNpal[ifct].nidPair[1] + idx_init << "," ;
    outFile << myNNpal[ifct].nidPair[2] + idx_init << ""  ;

    for (int j = 0 ; j < myNNpal[ifct].vidPairMax ; ++j) {
      outFile  << "," << myNNpal[ifct].vidPair[j] + idx_init ;
    }
    outFile << "\n";
  }
  outFile.close();
}
