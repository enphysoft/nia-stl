//
/* Copyright © 2019-present, by Albert S. Kim
  Author: Albert S. Kim ( http://albertsk.org/ )
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
  output = mystlNNB.dat, facetNNB.dat, checkNNB.dat */



using namespace std;

string entree_fileSTL;
string defaultfileSTL  = "STL_INPUT.stl";
string output_fileSTL  = "STL_INPUT_NNBS.stl";
string output_fileSTL0 = "STL_INPUT_cnia0.stl";
string output_fileSTL1 = "STL_INPUT_cnia1.stl";
string dot = ".";
string arg1_root, arg1_extn;
string postfix_nnbed = "_nia";
string postfix_nnbed0 = "_cnia0";
string postfix_nnbed1 = "_cnia1";
string output_fileNNB = "NNBindex.dat";
string review_fileNNB = "NNBcheck.dat";
string facets_fileNNB = "NNBfacet.dat";
string runtst_NNBlist = "NNBlists.dat";
string title_stl;

int numline, numtri;
int i_fct, l_fct, k_fct;
int i_vtx, l_vtx, j_vtx;
int N_fct, N_fct_o10, numPair_i_l;
int nidPair_tmp[3], vetPair_tmp[3], matPair_tmp[3][3];
int *numFacets_i, *numVertices_i;

double vecA[3], vecB[3];

ifstream file;
string line;
string word;

struct triVertex {
  double vertexA[3];
  double vertexB[3];
  double vertexC[3];
};

struct stl_tri {
  int tid;               //triangle ID, sequential number
  double tnvec[3];         //normal vector
  triVertex tvtx;        //vertex points: vertexA, B, and C.
};

stl_tri *myFacet, *loadedFacet;

struct stl_nnb {
  int nidPair[3];
  int vidPair[30];
  int vidPairMax;
  int matPair[3][3];
  double mvtx[3][3];         //mvtx(i,j) i_th vertex's k_th component
};

stl_nnb *myNNpal, *myNNbor, *loadedNNbor;

struct stl_tri_nnb {
  stl_tri fct;
  stl_nnb nnb;
};

stl_tri_nnb *nnbedFacet;


void read_ascii_stl(string,int,stl_tri *,string);
void write_ascii_stl_w_nnb_all(string, stl_tri *, stl_nnb *, int, string, int );
void read_ascii_stl_w_nnb_all(string, stl_tri *, stl_nnb *, int, string);
