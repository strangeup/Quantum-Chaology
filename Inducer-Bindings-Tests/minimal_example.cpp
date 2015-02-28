// A test of arpack.hpp
// The input data was gibberish, just a test to check against the compiler the 
// matrix used is probably singular and too small so most likely wouldn't work anyway.
// However I couldnt get Inducers' code to give sensible eigenvalues even for matrices 
// I knew should give me finite results.

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/bindings/arpack.hpp>
#include <complex>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp> 

using namespace std;
const int N=6; //global size variable

int main()
{
	namespace ublas=boost::numeric::ublas;
	namespace arpack=boost::numeric::bindings::arpack; //to make things shorter

	ublas::vector<double> vec(N,0);
	vec(5)=0.81; vec(1)=0.57; vec(2)=(0); vec(0)=0; vec(3)=0;//near the vector

	ublas::mapped_matrix<double> mat(N,N);
	mat(1,5)=3;
	mat(5,1)=3;
	mat(2,2)=-1;
	mat(5,5)=2;
	mat(3,2)=0.5;
	mat(2,3)=0.5;//initalise mapped_vector (sym) (could also use almost any boost vector)

	for(int row=0; row<mat.size1();++row){ //print to terminal
		for(int col=0; col<mat.size2();++col)
			cout<<mat(row,col)<<" ";
		cout<<endl;
	}

	pyublasext::ublas_matrix_operator<ublas::mapped_matrix<double>,ublas::vector<double>,ublas::vector<double>,ublas::mapped_matrix<double> > op(mat);
	pyublasext::ublas_matrix_operator<ublas::mapped_matrix<double>,ublas::vector<double>,ublas::vector<double>,ublas::mapped_matrix<double> >* empt=0;
	//create conatiners for the matrix using matrix_operator header

	arpack::results<ublas::vector<complex<double > > > myresults; //create results container using arpack header
	//needs complex type
  	arpack::arpack_mode mode=arpack::REGULAR_NON_GENERALIZED; //use enu type as mode that is wanted
  	double shift(4.1); //set shift
  	arpack::perform_reverse_communication(op,empt,arpack::SHIFT_AND_INVERT_GENERALIZED,shift,4,6,myresults,vec,arpack::LARGEST_MAGNITUDE,1,0);

  	ublas::vector<complex<double> > evec1=(myresults.m_ritz_vectors).at(0);

  	cout<<"Vectors:";
  	for(ublas::vector<complex<double> >::iterator it= evec1.begin(); it!=evec1.end();++it) //print to terminal
  	 	cout<<real(*it)<<" ";
  	cout<<endl<<"Values: ";
  	for(vector<complex<double> >::iterator it=(myresults.m_ritz_values).begin(); it!=(myresults.m_ritz_values).end();++it)//print to terminal
  		cout<<real(*it)<<" ";
  	cout<<endl;


}
