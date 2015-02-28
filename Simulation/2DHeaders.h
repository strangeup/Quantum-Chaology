#ifndef MARIACHI_H
#define MARIACHI_H

//list of headers
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdlib.h> 
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <time.h>
#include <sstream>
#include <utility>
#include <boost/unordered_map.hpp>
#include <boost/scoped_array.hpp>
#include <boost/foreach.hpp>
#include <algorithm>
#include <iterator>




//shorter namespaces for clarity
using namespace std;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
typedef boost::unordered_map <pair<int,int>,double> sparse_matrix;

//----------------------------------------------------------------------------------
//Global Variables
//----------------------------------------------------------------------------------

// Global properties

const string well ="stadium"; //defines potential
extern char symmetry; //b - both or x- reflection x=0 y- reflection y=0 gives n is none
const char symmetryDefault='b'; //x=0 axis
const bool debugging = false;
const int loopStart=0;
const int loopEnd=8;

// Billiard Variables

//Shared Variables
extern int Nx; // defines number of points in x
extern int Ny; // defines number of points in y
const int NxStart=1000;//300 //won't handle 100 million points
const int NyStart=800; //irrational ratio so no degeneracy
const int radiusWell=53; //used for circular wells
const double midY=0;//to set the midpoint of triangle/circle
const double midX=0;

//Sinai Specific Variables
extern int R; //radius of sinai
const int RStart=3;

//Triangle Specific Variables
extern int L;
const int LStart=49;

//Smiloid Specific Variables
const int Reye=3;
const int Rmouth=5;

//Limaçon Specific Variables
const int b=90; //r=b gives cardioid
const int a=90; //r=2b gives limaçon trisectrix

//Other global properties
const double pi=3.141592653589793238;

#endif
