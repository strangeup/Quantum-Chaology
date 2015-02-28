## Boost Binding for ARPACK solver

### Description
Links and implements the ARPACK symmetric solver for a modest sized sparse matrix, with
shift invert.

### Software Requirements

Needs correct installation of umfpack, AMD, LAPACK, BLAS and ARPACK and up to date BOOST
libraries

* LAPACK can easily be downloaded on a UBUNTU/MINT system using apt-get install liblacpack3 
liblapacke and liblapack-dev 

* BLAS can similarly be set-up using apt get install libblas3 and libblas-dev

* BOOST can be set up following instructions from [here]
(http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html)

* AMD and UMFPACK can easily be set up using apt-get install libsuitesparse-metis 
libsuitesparse-metis-dev

* ARPACK can be set up most easily using ARPACK-ng found [here] 
(http://forge.scilab.org/index.php/p/arpack-ng/) this is a version of ARPACK which has been 
modified slightly with a few bug fixes and an included make file for easy installation

* BOOST-UMFPACK bindings: these allow the usage of sparse routines written in c for the shift 
inverion technique. The bindings used can be found [here] 
(http://mathema.tician.de/software/boost-numeric-bindings/) and for quick use the boost
folder can just be merged into your boost folder (or use the make file) which should get 
you set up. These bindings will also allow you to use the files found elswhere on this git 
which require boost LAPACK bindings.

That is a lot of set-up so good luck!

### Usage

Arguments: 
* a matrix which you wish to solve as a boost::compressed_matrix:
ublas::compressed_matrix<double, ublas::column_major, 0,ublas::unbounded_array<int>,
				ublas::unbounded_array<double> > mat

* a vector for results: 
	vector<double> evals

* a vector for eigenvectors:
	vector<vector<double> > > evecs

* number of eigenvalues/eigenvectors:
	int nev

* length of arnoldi factorisation (look this up in [ARNOLDI documentation](http://www.caam.rice.edu/software/ARPACK/) for more 
information)	
	int ncv 
NB: should be greater than twice the number of eigenvalues but less than size of matrix

* shift (which eigenvalue you want to get values around eg. 3)
	double sigma

* [OPTIONAL] max number of iterations
	int maxit

* [OPTIONAL] tolerance (look this up in [ARNOLDI documentation](http://www.caam.rice.edu/software/ARPACK/) for more information)
	int tol
NB: tol=0 means it uses the machine precision
	

