#ifndef BILLIARD_FUNCTIONS_H
#define BILLIARD_FUNCTIONS_H
#include "2DHeaders.h"

//----------------------------------------------------------------------------------
//Function Declarations
//----------------------------------------------------------------------------------
//Potential Well functions
int matrixElement(const int& i, const int& j, const int& k, const int& l) ; //returns unreduced matrix element
bool pointInWell(const int& y,const int& x); //checks if point in the potential well 
int numberPoints(); //returns number of points in the potential well
void wrap(vector<int> &indexLookup); //wrapping function (passed by ref)
void wrap( boost::unordered_map<int,int>& indiceLookup); //wrapping function (passed by ref)
int unwrap(int gridPoint, char index); //returns index specified by combined index ijkl

//Output Functions
void printEigenvalues(const vector<double>& vec, const string name, const char delimiter); //prints vector
template <typename boost_ublas_matrix>
void printMatrix(const boost_ublas_matrix& mat, const string outputFilename, const char delimiter=',');
void printWell(const string filename); //prints potential well
void printWell(); //shows the potential well used in terminal
void printEigenvectors(const ublas::matrix<double>& mat, const string name, const int eigenNumber, const char delimiter);
void printEigenvectors(const vector<vector <double> >& evecs, const char delimiter); //prints out vectors

//Utility Functions
void text(const string str); //flush text in terminal
string int2string(const int Integer); //converts int to string
string namingConvention(); //naming convention for files
string file_extension(const char& delimiter); //returns the file extension deduced from delimiter

//Matrix Assignment Functions
double reducedElement(const int& i, const int& j, const int& k, const int& l); //reduced matrix element
void nonZeroElem(vector<pair<int,int> >& nonZeros); //non zero i-j combinations
template <typename boost_ublas_matrix>
void make_Matrix(boost_ublas_matrix &Mat, const vector< pair<int,int> >& nonZeros);

//Correlation Specific Functions
int periodicBounds(int index); //maps generic index into periodic matrix index
void nonzeroMaps(const vector<pair <int,int> > &nonZeros, sparse_matrix &map_init); //"sparse" matrix
void crossCorrelation(const sparse_matrix &map_init, const sparse_matrix &map_cmpr, sparse_matrix &map_cross);//calculates cross correlation
double autoCorrelation(const sparse_matrix &map_init); //calculates autocorrelation
void printCorrelation(const sparse_matrix &map_cross, const double& autoCorr, const char& flag);
void periodicSolutions(sparse_matrix &map_cross,int &x,int &y, const double &cross, const int N);

//Arpack Specific Functions
void arpack_solve();

//----------------------------------------------------------------------------------
//POTENTIAL WELL FUNCTIONS

bool pointInWell(const int& y,const int& x) // row, column format
//checks if point is interior (ie. non zero)
{	
	if ( symmetry =='n') // this if excludes boundary points
	  if ( y==0 || x==0 )
		  return false;
	if ( symmetry =='x')
	  if ( y==0 )
	  {
		  return false;
	  }
	if ( symmetry =='y')
	  if ( x==0 )
		  return false;
	  
	if ( well=="rectangle" )
	{
	if ( abs(y)<Ny && abs(x)<Nx )
		return true;
	else 
		return false;
	}
	
	if ( well=="stadium" ) // curved edge and straight edge
	{
	if ( (abs(x)<Nx-Ny && abs(y)<Ny) || pow(x-Nx+Ny,2)+pow(y,2) <pow(Ny,2) ) 
	//Ny is radius (number, not index)
		return true; 
	else 
		return false;
	}
	
	
	else if(well == "sinai") //circle in rectangle
	{
	if(pow(abs(x)-midX,2)+pow(abs(y)-midY,2) <= pow(R,2))
		return false;
	if ( (abs(x)<Nx && abs(y)<Ny))
		return true;
	else 
		return false;
	}
	
	else if (well =="triangle")// equilateral triangle in rectangle
	{
	// note equations defining equilateral triangle are:
	// y=-sqrt(3)x-sqrt(3)*l/3 ; y=sqrt(3)x-sqrt(3)*l/4 ; y=sqrt(3)*l/4 

	if(abs(y) <= sqrt(3)*L/4 && abs(y) >= sqrt(3)*(abs(x))+1+midY-sqrt(3)*L/4)
	  return false;
	if(abs(y) <= midY+sqrt(3)*L/4 && abs(y) >= sqrt(3)*(abs(x))+1+midY-sqrt(3)*L/4)
	  return false;	
	if(abs(x) < Nx && abs(y)< Ny) return true;
	else return false;
	}
	
	else if (well =="interiorTriangle")// equilateral triangle
	{
	// note equations defining equilateral triangle are:
	// y=-sqrt(3)x-sqrt(3)*l/3 ; y=sqrt(3)x-sqrt(3)*l/4 ; y=sqrt(3)*l/4 

	if(abs(y) <= midY-1+R+sqrt(3)*L/4 && abs(y) >= sqrt(3)*(abs(x-midX))+midY+R-1-sqrt(3)*L/4)
	  return true;
	else return false;
	}
	
	else if (well == "lymacon") // interior of lymacon
	{
	if(pow((pow(y-a,2)+pow(x,2)-a*(y-a)),2)-pow(b,2)*(pow(y-a,2)+pow(x,2))<0)
	   return true;
	}
	
	else if (well == "smiloid") // circle and quarter circle
	{
	if( (pow(x-2*R/5,2)+pow(y-4*R/5,2)-pow(Reye,2)<=0))
	  return false;
	else if( (pow(x,2)+pow(y-4*R/3,2)-pow(Rmouth,2)<=0) && y>4*R/3)
	return false;
	else if (pow(y-R,2)+pow(x,2)-pow(R,2)<0)
	   return true;
	}
	
	else if (well == "sinaiTriangle") // equilateral triangle in circle
	{
	if(abs(y) <= midY-1+R+sqrt(3)*L/4 && abs(y) >= sqrt(3)*(abs(x-midX))+midY+R-1-sqrt(3)*L/4)
	  return false;
	else if (pow(y-R,2)+pow(x,2)-pow(R,2)<0)
	   return true;
	}
}


int matrixElement(const int& i, const int& j, const int& k, const int& l) 
// in ijkl for the reflection
// points ij and kl in potential by construction but i-1,j etc. not necessarily
{
	if (i==k && j==l) 
		// if r in original big matrix (with zeros) = c
		return 4; 
	if  (abs(i-k)==1 && j==l) 
		// requires separation of two statements
		return -1;
	if  (i==k && abs(j-l)==1) 
		return -1;	
	return 0;
}	 //end of matrixElement


int numberPoints()
//returns number of interior points
{
int count=0;
for(int x=0; x<Nx; ++x)
{
	for(int y=0; y<Ny; ++y)
	{
		if(pointInWell(y,x)==true)
			++count;
	}
}
return count;
}//end of number of points


void wrap(vector<int> &indiceLookup)
//takes a vector and fills it with 'wrapped' indices of interior points
{
	for(int gridPoint=0; gridPoint<Nx*Ny; ++gridPoint) //initialise indiceLookup
	{
		if (pointInWell(gridPoint / Nx, gridPoint % Nx)==true)
			//only really needs one condition
			indiceLookup.push_back(gridPoint); //adds element to end
			// r*Nx*Ny+c ties everything into one integer to keep things neat
	}
} //end of wrap


void wrap( boost::unordered_map<int,int> &indiceLookup)
//takes a vector and fills it with 'wrapped' indices of interior points
{
	int count=0;
	for(int gridPoint=0; gridPoint<Nx*Ny; ++gridPoint) //initialise indiceLookup
	{
		if (pointInWell(gridPoint / Nx, gridPoint % Nx)==true){
			indiceLookup[gridPoint]=count; //adds element to end
			++count;
		}
			// r*Nx*Ny+c ties everything into one integer to keep things neat
	}
} //end of wrap


int unwrap(int gridPoint, char index) 
//uses combined index ijkl to return single index
{
	switch (index)
	{
		case 'i': //int i= r / Nx
			return gridPoint / Nx;
		case 'j': //int j= r % Nx
			return gridPoint % Nx;
		default:
			cout << "This is not a valid Index" <<endl;
			exit(-1);
			break;	
	}//end of switch
}//end of unwrap

//------------------------------------------------------------------------------
//UTILITY FUNCTIONS

void text(const string str) 
//utility function
{
	cout << "\r" << str;
	cout << flush << "\r" << "                                                                             ";
}

string int2string(const int Integer)
{
    stringstream ss; //define sstream
    string Str; // define string
    ss << Integer;
    ss >> Str;
return Str;
}


string namingConvention()
//contructs name that includes all information for billiard
{
  string name;
  
  switch (symmetry)
  {
    case 'b':
      name="even-even-"+well+"-";
    break;
    case 'n':
      name="odd-odd-"+well+"-";
    break;
    case 'x':
      name="even-odd-"+well+"-";
    break;
    case 'y':
      name="odd-even-"+well+"-";
    break;
  }
  
  if (well=="rectangle" || well=="stadium")
     name+=int2string(Nx)+"-"+int2string(Ny);
  else if (well=="sinai"||well=="smiloid")
    name+=int2string(Nx)+"-"+int2string(Ny)+"-"+int2string(R);
  else if (well=="triangle"|| well=="interiorTriangle")
    name+=int2string(Nx)+"-"+int2string(Ny)+"-"+int2string(L)
		+"-"+int2string(midX)+"-"+int2string(midY);
  else if (well=="lymacon")
    name+=int2string(Nx)+"-"+int2string(Ny)+"-"+int2string(a)
		+"-"+int2string(b);
  else if (well=="sinaiTriangle")
    name+=int2string(Nx)+"-"+int2string(Ny)+"-"+int2string(R)
		+"-"+int2string(L)+"-"+int2string(midX)+"-"+int2string(midY);
                
return name;
}//end of namingConvention


string file_extension(const char& delimiter)
//returns the file extension deduced from delimiter
{	
	string extension;
  	switch (delimiter){ //adds file extension according to delimiter
  		case ',':
  			extension=".csv"; break;
  		case ' ':
  			extension=".dat"; break;
  		default:
			extension=".txt"; break;  	
	}
	return extension;
}

//----------------------------------------------------------------------------------
//OUTPUT FUNCTIONS
void printWell()
//prints well to terminal
{
	for(int y=-Ny; y<Ny+1; ++y) 
	//includes a zero border to illustrate well
	{
		for(int x=-Nx; x<Nx+1; ++x)
			cout << pointInWell(abs(y),abs(x))<<"\t";
	cout<< endl;
	}
}//end of sayWell


void printWell(const string filename)
// output to file
{
ofstream out(filename.c_str());
for(int y=-Ny; y<Ny+1; ++y) 
	//includes a zero border to illustrate well
{
	for(int x=-Nx; x<Nx+1; ++x)
		out << pointInWell(abs(y),abs(x))<<" ";
	out<< endl;
}
	out.close();
}//end of print Well


void printEigenvectors(const ublas::matrix<double>& mat, const string name, const int eigenNumber, const char delimiter=',')
//prints out a particular eigenvector grid
{
  string extension(file_extension(delimiter));
  int count(0);
  ofstream out(name.c_str()+extension);
  
  for(int y=0; y<Ny; ++y)
  {
	for(int x=0; x<Nx; ++x)
	{	
	  if (pointInWell(y, x)==true)
	    {
	    out << pow(mat(count,eigenNumber),2) << (x==Nx-1 ? '\n' : delimiter);
	    ++count;
	    }
	  else 
	    out << 0. << (x==Nx-1 ? '\n' : delimiter);
	}
  }
  out.close();
} //end of print vectors


void printEigenvectors(const vector<vector<double> >& evecs, const char delimiter=',')
//prints out a particular eigenvector grid
{
  string extension(file_extension(delimiter));

  for(auto it=evecs.begin();it<evecs.end();++it)
  {
	  int count(0);
	  ofstream out((int2string(distance(evecs.begin(),it))+extension).c_str());
	  
	  for(int y=0; y<Ny; ++y)
	  {
		for(int x=0; x<Nx; ++x)
		{	
		  if (pointInWell(y, x)==true)
		    {
		    out << pow((*it).at(count),2)	<< (x==Nx-1 ? '\n' : delimiter);
		    ++count;
		    }
		  else 
		    out << 0. << (x==Nx-1 ? '\n' : delimiter); //prints \n if end of line, delim otherwise
		}
	  }
	   out.close();
  }
} //end of print vectors


void printEigenvalues(const vector<double> & vec, const string name, const char delimiter=',')
// print eigenvalues file
{
	// output to file
	ofstream out(name.c_str()+file_extension(delimiter));
	for(auto it=vec.begin(); it<vec.end(); ++it) // print numerical eigenvalues
		out << (*it) << endl;
	out.close();
} //end of print eigenvalues


template <typename boost_ublas_matrix>
void printMatrix(const boost_ublas_matrix& mat, const string outputFilename,const char delimiter=',')
//prints ublas matrix
{
	const int rows = mat.size1();
	const int cols = mat.size2();
	// output to file
	ofstream out(outputFilename.c_str()+file_extension(delimiter));
	for(int r = 0; r < rows; ++r)
		for(int c = 0; c < cols; ++c)
			out << mat(r, c) << (c==cols-1 ? '\n' : delimiter);

	out.close();
} //end of print matrix

//----------------------------------------------------------------------------------
// MATRIX ASSIGNMENT FUNCTIONS


double reducedElement(const int& i, const int& j, const int& k, const int& l)
//element for each symmetry case
{	
	// definitions for use in matrixElement
	double element, norm;
        
	switch (symmetry)
	{
	case 'n': //matrix element for no axis of symmetry
		element = matrixElement( i, j, k, l);
	break;
	case 'x': //matrix element for x=0 axis of symmetry
		norm=0.5;	// normalise
		if (j==0) 
			norm*=sqrt(2)/2;
		if (l==0)
			norm*=sqrt(2)/2;
		element=norm*( 
		+ matrixElement( i , j, k, l) + matrixElement( i ,-j, k, l)
		+ matrixElement( i , j, k,-l) + matrixElement( i ,-j, k,-l) 
		); //end of element assignment
		break;
	case 'y': //matrix element for y=0 axis of symmetry
		norm=0.5; //normalise		
		if (i==0)
			norm*=sqrt(2)/2;
		if (k==0)
			norm*=sqrt(2)/2;
		element=norm*( 
		+ matrixElement( i , j, k, l) + matrixElement(-i , j, k, l)
		+ matrixElement( i , j,-k, l) + matrixElement(-i , j,-k, l)
		); //end of element assignment
		break;
	case 'b': //matrix element for both  symmetry axes
		norm = 0.25; //normalise			
		if (i==0)
			norm*=sqrt(2)/2;
		if (j==0)
			norm*=sqrt(2)/2;
		if (k==0)
			norm*=sqrt(2)/2;
		if (l==0)
			norm*=sqrt(2)/2;
		element=norm*(
		// tilde Hijkl - reflected in x
		+ matrixElement( i , j, k, l) + matrixElement( i ,-j, k, l) 
		+ matrixElement( i , j, k,-l) + matrixElement( i ,-j, k,-l)
		// tilde H-ijkl 
		+ matrixElement(-i , j, k, l) + matrixElement(-i ,-j, k, l) 
		+ matrixElement(-i , j, k,-l) + matrixElement(-i ,-j, k,-l)
		// tilde Hij-kl
		+ matrixElement( i , j,-k, l) + matrixElement( i ,-j,-k, l) 
		+ matrixElement( i , j,-k,-l) + matrixElement( i ,-j,-k,-l)
		// tilde H-ij-kl
		+ matrixElement(-i , j,-k, l) + matrixElement(-i ,-j,-k, l) 
		+ matrixElement(-i , j,-k,-l) + matrixElement(-i ,-j,-k,-l)
		); //end of matrix element input
		break;
		default:
		    cout<<"invalid symmetry option" <<endl;
		    exit(EXIT_FAILURE);
		break;
	}		
	return element;
} //end of reduced element


void nonZeroElem(vector<pair <int,int> > &nonZeros) //non zero i-j combinations
{
    for (int i=0; i<Ny; ++i) //loop through all i
    {
        for(int j=0; j<Nx; ++j) //loop through all j
        {
            if (pointInWell(i,j)) //if point is interior to well
                nonZeros.push_back(pair<int,int>(i*Nx+j,i*Nx+j));
            // self interaction terms
            if (pointInWell(i,j) && pointInWell(i+1,j)) //if both points interior
                nonZeros.push_back(pair<int,int>(i*Nx+j,(i+1)*Nx+j));
             // interaction terms from above
            if (pointInWell(i+1,j) && pointInWell(i,j)) //if both points interior
                nonZeros.push_back(pair<int,int>((i+1)*Nx+j,i*Nx+j));
            // interaction terms from below
            if (pointInWell(i,j) && pointInWell(i,j+1)) //if both points interior
                nonZeros.push_back(pair<int,int>(i*Nx+j,i*Nx+j+1));
            // interaction terms from right
            if (pointInWell(i,j+1) && pointInWell(i,j)) //if both points interior
                nonZeros.push_back(pair<int,int>(i*Nx+j+1,i*Nx+j));
            // interaction terms from left
        }
    }
}


template <typename boost_ublas_matrix>
void make_Matrix(boost_ublas_matrix &Mat, const vector< pair<int,int> > &nonZeros)
// creates then solves the matrix faster than solveMatrix()
{
        boost::unordered_map<int,int> indices; //empty vector
        wrap(indices); //add indexing to empty vector
        text("assign matrix for:" + namingConvention());
        int i,j,k,l,row,col;
        for(auto it = nonZeros.begin() ; it != nonZeros.end(); ++it) //iterate rows
        {
            i=(*it).first / Nx; j=(*it).first % Nx;
            k=(*it).second / Nx; l=(*it).second % Nx;
            row=indices[(*it).first];
            col=indices[(*it).second];
            Mat(row,col)=reducedElement(i,j,k,l); //assign element to matrix
        }
        
        if (debugging==true)
            printMatrix(Mat,"MatrixNew-"+namingConvention()); //for debugging
} //end of construction

//----------------------------------------------------------------------------------
//Correlation Specific Functions

void nonzeroMaps(const vector<pair <int,int> > &nonZeros, sparse_matrix &map_init)
//returns map containing nonzero matrix elements
{
    int NumPoints(numberPoints());
    vector<int> indices;
    wrap(indices); //sets up indexing for use in the reducedElement function
    int count=1;
        for(auto it = nonZeros.begin() ; it != nonZeros.end(); ++it) //iterate rows
        {
            int i=(*it).first / Nx; int j=(*it).first % Nx;
            int k=(*it).second / Nx; int l=(*it).second % Nx;
            int row=distance(indices.begin(),find(indices.begin(),indices.end(),(*it).first));
            int col=distance(indices.begin(),find(indices.begin(),indices.end(),(*it).second));
            int element=reducedElement(i,j,k,l); //define element
            map_init[pair<int,int>(row,col)]=element; //assign element to map_init
        }
}

int periodicBounds(int index)
//returns pair containing mapped row and column to periodic matrix
{
    int NumPoints(numberPoints());
    index = index % numberPoints(); //in 'base' numpoints remainder is new position 
    
    if (index < 0)
        index = numberPoints()+index; //translate negative points
    return index;
}


void periodicSolutions(sparse_matrix &map_cross,int &x,int &y, const double &cross, const int N)
//correct boundary conditions for -N/2 < x,y < N/2
{
   map_cross[pair<int,int>(x,y)]+=cross;       
   map_cross[pair<int,int>(x-N,y)]+=cross;
   map_cross[pair<int,int>(x-N,y-N)]+=cross; //needs work!
   map_cross[pair<int,int>(x-N,y+N)]+=cross;
   map_cross[pair<int,int>(x+N,y)]+=cross;
   map_cross[pair<int,int>(x+N,y-N)]+=cross;
   map_cross[pair<int,int>(x+N,y+N)]+=cross;
   map_cross[pair<int,int>(x,y-N)]+=cross;
   map_cross[pair<int,int>(x,y+N)]+=cross;
 }

void crossCorrelation(const sparse_matrix &map_init, const sparse_matrix &map_cmpr, sparse_matrix &map_cross)
{
    int count(1); //as a progress parameter
    const double MapSize=map_init.size(); //for convenience
    int NumPoints=numberPoints();
     BOOST_FOREACH(sparse_matrix::value_type it_init, map_init) //iterates through map
     {
         text("calculating <g(R)g(r+R)>: "+int2string(count)+"/"+int2string(MapSize));
         BOOST_FOREACH(sparse_matrix::value_type it_cmpr, map_cmpr)//iterates through map
         {
               int x (it_init.first.first-it_cmpr.first.first); //x seperation of two points
               int y (it_init.first.second-it_cmpr.first.second); //y seperation of two points
               double cross (it_init.second*it_cmpr.second/MapSize);
               periodicSolutions(map_cross,x,y,cross,NumPoints);
               //initialise map of cross correlations
          }
          ++count;
      }
cout<<endl;
}


double autoCorrelation(const sparse_matrix &map_init)
{
     double average;
     const double MapSize=map_init.size(); //for convenience
     BOOST_FOREACH(sparse_matrix::value_type it_init, map_init)//iterates through map
        average+=it_init.second/MapSize; //iterates through and averages
     return average;
}


void printCorrelation(const sparse_matrix &map_cross, const double& autoCorr, const char& flag)
{
    string filename="correlation-"+namingConvention()+".csv";
    int NumPoints=numberPoints(); //for convenience
    ofstream out(filename.c_str()); //open file to write to
    
    if (flag=='A') //2D array print
    {
        for(int row=-NumPoints/2; row<NumPoints/2+1;++row)
        {
            for(int col=-NumPoints/2; col<NumPoints/2+1;++col)
            {
            if (map_cross.find(pair<int,int>(row,col)) != map_cross.end())
                out <<map_cross.at(pair<int,int>(row,col))<<" ";
            else
                out<<"0.0"<<" ";
            }
            out<<endl;
        }
    }
    else if (flag=='L') //list print
    {
        for(int row=-NumPoints/2; row<NumPoints/2+1;++row)
        {
            for(int col=-NumPoints/2; col<NumPoints/2+1;++col)
            {
            if (map_cross.find(pair<int,int>(row,col)) != map_cross.end())
                out <<row<<","<< col<<"," 
                    << map_cross.at(pair<int,int>(row,col)) << endl;
            else
                out<<row<<","<<col<<","<<"0.0"<<endl;
            }
        }
     }
     else
         cerr << "Invalid print option"<<endl;
    out.close();   //close file
}

#endif