/* ARPACK SOLVER FOR SYMMETRIC EIGENPROBLEM IN MODE 3
 * we want to solve A*x = lambda*x in shift-invert mode, where
 * OP = (inv[A - sigma*I]) and  B = I.
 * Use mode 3 of DSAUPD
 *
 * Routines called:
 *     dsaupd  ARPACK reverse communication interface routine.
 *     dseupd  ARPACK routine that returns Ritz values and (optionally) Ritz
 * 			   vectors.   
 * Requires functioning LAPACK, ARPACK, BLAS, AMD and UMFPACK 
 *
 */
#ifndef BINDINGS_ARPACK_SYMMETRIC_SHIFT
#define BINDINGS_ARPACK_SYMMETRIC_SHIFT

#include <iostream>
#include <fstream>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/scoped_array.hpp>
#include <vector>

namespace ublas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;
using namespace std;


extern "C"
{
	// ARPACK functions
	void dsaupd_(int *IDO, char *BMAT, int *N, char *WHICH, int *NEV,
				 double *TOL, double *RESID, int *NCV, double *V, int *LDV,
				 int *IPARAM, int *IPNTR, double *WORKD, double *WORKL,
				 int *LWORKL, int *INFO);
	
	void dseupd_(bool *RVEC, char *HOWMNY, int *SELECT, double *U, double *Z,
				 int *LDZ, double *SIGMA, char *BMAT, int *N, char *WHICH,
				 int *NEV, double *TOL, double *RESID, int *NCV, double *V,
				 int *LDV, int *IPARAM, int *IPNTR, double *WORKD,
				 double *WORKL, int *LWORKL, int *INFO);
}


void begin_rev_communication(ublas::compressed_matrix<double, ublas::column_major, 0,ublas::unbounded_array<int>,ublas::unbounded_array<double> >& mat,
							 //WARNING CHANGES MAT (then changes back)
                             vector<vector<double> >& evecs, //results holder
                             vector<double>& evals, //results holder
                             int nev, //number eigenvalues 
                             int ncv, //length of arnoldi factorisation
                             double sigma, //shift
                             int maxitr=5000, //with default 5000
                             double tol=0 //with default machine precision
                             )

{							 
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * 						  BEGIN ARPACK INTERFACE						   *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	// Translation to C and implementation of symmetric routine courtesy of 
	//Paule Searle
	int N=mat.size1();	// get dimension of matrix
	for (int i=0; i<N;++i)
		mat(i,i)-=sigma; //minus identity*sigma
	// stuff for dsaupd:
	int ido = 0;		// used for reverse communication, initially set to 0
						// ido = -1 or +1: compute Y = OP * X
						// ido = 99: finished reverse communication
	char bmat = 'I';	// standard (non-generalised) eigenvalue problem
	char which[2] = {'L', 'M'};	// find evals of largest magnitude
	int info = 0;		// generate random vector in dsaupd to start Arnoldi
						// iteration
	int ldv = N;		// leading dimension of v
	boost::scoped_array<double> v(new double [ncv*N]);	// the ncv columns of v contain the Lanczos basis
						// vectors
	boost::scoped_array<double> workd (new double [3*N]);	// distributed array to be used in the basic Arnoldi
						// iteration for reverse communication. This is used in
						// the subroutine dseupd
	int lworkl = ncv * (ncv+8);	// dimension of workl
	boost::scoped_array<double> workl(new double [lworkl]);	// work array used in dsaupd as workspace
	boost::scoped_array<double> resid(new double[N]);	// resid contains the final residual vector
	int iparam[11];		// array containing options, we set some below
	int ipntr[11];		// array containing pointers marking the starting
						// locations in the workd and workl arrays for
						// matrices/vectors used by the Lanczos iteration
	
	// iparam settings, only need to set certain ones:
	iparam[1-1] = 1;	// use exact shifts
	iparam[3-1] = maxitr; // max number of iterations
	iparam[7-1] = 3;	// mode 3 is shift-invert mode	
	
	umf::numeric_type<double> Numeric;
	//factorise mat using symbolic analysis
	umf::factor(mat, Numeric);

	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * REVERSE COMMUNICATION LOOP:											   *
	 * Repeatedly call the routine dsaupd and take actions indicated by		   *
	 * parameter ido until either convergence is indicated or maxitr has been  *
	 * exceeded.															   *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	do
	{
		dsaupd_(&ido, &bmat, &N, which, &nev, &tol, resid.get(), &ncv, v.get(), &ldv,
				iparam, ipntr, workd.get(), workl.get(), &lworkl, &info);
		
		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Perform matrix vector multiplication								   *
		 * 		y <--- OP*x													   *
		 * Need to use a matrix vector multiplication routine that takes	   *
		 * workd(ipntr(1)) as the input and returns the result to			   *
		 * workd(ipntr(2))													   *
		 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
		
		// use pointer arithmetic and the fact the fortran arrays start at 1
		double *input = workd.get() + ipntr[1-1] - 1;
		double *output = workd.get() + ipntr[2-1] - 1;
		
		ublas::vector<double> inVec(N);
		copy(input, input+N, inVec.begin());
		
		// solve for outVec
		ublas::vector<double> outVec(N);
		umf::solve(mat, outVec, inVec, Numeric);

		copy(outVec.begin(), outVec.end(), output);

	}
	while(ido == 1 || ido == -1);
	
	// free memory
	Numeric.free();
	
	// stuff for dseupd
	boost::scoped_array<int> select(new int[ncv]);	// used as a workspace for reordering the Ritz values
	boost::scoped_array<double> d(new double [nev]);		// d represents the Ritz values of OP computed by dsaupd
						// transformed to those of the original eigensystem
						// A*z = lambda*B*z
	bool rvec = true;	// compute Ritz vectors
	char howmny = 'A';	// compute all (nev) Ritz vectors
	
	dseupd_(&rvec, &howmny, select.get(), d.get(), v.get(), &ldv, &sigma, &bmat, &N, which, &nev,
			&tol, resid.get(), &ncv, v.get(), &ldv, iparam, ipntr, workd.get(), workl.get(), &lworkl,
			&info);

	 evals=vector<double> (nev,0); //pre-initialise for convenience
	 evecs=vector<vector<double> >(nev,vector<double>(N,0));

	// assign evals
	 for(int i = 0; i < nev; ++i) evals.at(i)=d[i];
	 // assign evecs: have nev of them
	 for(int vecnum = 0; vecnum < nev; ++vecnum){ //loop vectors
	 	for(int i=0; i<N;++i) //loop through elements of vectors
	 		(evecs.at(vecnum)).at(i)=v[vecnum*N+i];
	 }
	 for (int i=0; i<N;++i) //put matrix back to normal
		mat(i,i)+=sigma; //minus identity*sigma
}

#endif 