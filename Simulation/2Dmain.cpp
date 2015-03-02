#include "2DHeaders.h"
#include "2Dfunctions.h"
#include "arpack_symm_si.hpp"

//----------------------------------------------------------------------------------
// Global Variables
int Nx;
int Ny;
int R;
int L;
char symmetry;
//----------------------------------------------------------------------------------
// Declarations
void correlationFunction(); //for correlation data
void spectralAnalysis(); //for spectral data
//----------------------------------------------------------------------------------

int main()
{
  // set as start values
  R=RStart;
  Nx=NxStart;
  Ny=NyStart;
  L=LStart;
  int N=numberPoints();
   
  symmetry='n';
  //printWell("well-"+namingConvention()+".csv");
  cout<<"Number of Points: "<<numberPoints()<<endl;
  cout<<"Billiard: "<< well<< endl;
  arpack_solve();
  //spectralAnalysis();
  //correlationFunction();
}//end of main

//----------------------------------------------------------------------------------

void spectralAnalysis()
{
    int NumPoints=numberPoints();
    vector<pair <int,int> > nonZeros;
    nonZeroElem(nonZeros);
    vector<double> evals(NumPoints); //empty vector
    ublas::matrix<double, ublas::column_major> Mat (NumPoints,NumPoints);
    Mat.clear(); //sets all elements == 0

    make_Matrix(Mat,nonZeros); //assigns matrix
    if (debugging==true)	
    	printMatrix(Mat,"matrix");
    lapack::syev('N', 'U', Mat, evals, lapack::optimal_workspace());//solve, N means no eigenvectors
    printEigenvalues(evals,namingConvention()+".csv");
}

void correlationFunction() //the function that calls everything
{
 vector < pair<int,int> > nonZeros; //use a vector
 sparse_matrix map_init; //map to hold "sparse" matrix
 sparse_matrix map_cross; //map to hold "sparse" matrix
 sparse_matrix map_autoCorr; //map to hold "sparse" matrix
  
 nonZeroElem(nonZeros);
 nonzeroMaps(nonZeros,map_init); //fill comparison map
 
 crossCorrelation(map_init,map_init,map_cross); //calculate cross correlation
 const double autoCorr(autoCorrelation(map_init)); //calculate autocorrelation
 printCorrelation(map_cross,autoCorr,'L'); //print the correlationfunction

}//end of correlationFunction

void arpack_solve()
{
vector<vector<double> > evecs;
vector<double> evals;
vector<pair <int,int> > nonZeros;
int N=numberPoints();
ublas::compressed_matrix<double, ublas::column_major, 0,ublas::unbounded_array<int>,ublas::unbounded_array<double> >  Mat(N, N);
nonZeroElem(nonZeros); //assigns non zeros vector
make_Matrix(Mat,nonZeros); //assigns matrix
nonZeros.clear();
text("Beginning Reverse Communication Loop");
begin_rev_communication(Mat,evecs,evals,5,25,3.); //call arnoldi factorisation
printEigenvalues(evals,namingConvention());
printEigenvectors(evecs);
}