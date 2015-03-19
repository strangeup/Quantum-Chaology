#ifndef BOX_COUNTER_H 
#define BOX_COUNTER_H

//Headers
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <utility>
#include <boost/unordered_set.hpp>

namespace box_counting{

using namespace std;

class grid_pt{
private:
	double x,y;
public:
//Constructors and Destructor
	grid_pt() {x=0; y=0;}
	grid_pt(double x1, double y1){x=x1; y=y1;}
	~grid_pt() {}
//Operators
	double& operator [] (int comp){
		switch(comp){
		case 1: return x; break;
		case 2: return y; break;
		default: cerr<<"Index out of range\n"; exit(-1);
	}}	
	double operator [] (int comp) const{
		switch(comp){
		case 1: return x; break;
		case 2: return y; break;
		default: cerr<<"Index out of range\n"; exit(-1);
	}}
	grid_pt operator + (const grid_pt pos1)const {return grid_pt(x+pos1.x,y+pos1.y);}
	grid_pt operator - (const grid_pt pos1)const {return grid_pt(x-pos1.x,y-pos1.y);}
	bool   operator == (const grid_pt pos1)const {return (x==pos1.x && y==pos1.y);}
	bool   operator != (const grid_pt pos1)const {return !(x==pos1.x && y==pos1.y);}
//Member Functions
	double pt_length()const{return sqrt(pow(x,2)+pow(y,2));}
	grid_pt grid_posit(const grid_pt origin, const double s)const{
		return grid_pt(floor((x-origin.x)/s),floor((y-origin.y)/s));
	}//return posit rounded down to nearest gridpoint 
	pair<double,double> as_pair()const{return pair<double,double>(x,y);}
//Friend Functions
	friend grid_pt operator - (const grid_pt pos1){return grid_pt(-pos1.x,-pos1.y);}
	friend ostream& operator << (ostream& os, const grid_pt pt){
		os<<"("<<pt.x<<","<<pt.y<<")"; return os;}	
	friend istream& operator >> (istream& is, grid_pt& pt){
		char in_char;
		is>>in_char; if(in_char!='('&& in_char!='{') {cerr<<"Bad Input: '"<<in_char<<"'\n"; exit(-1);}
		is>>pt.x;	 if(is.fail()) {cerr<<"Bad Input (not a double) \n"; exit(-1);}
		is>>in_char; if(in_char!=',') {cerr<<"Bad Input: '"<<in_char<<"'\n"; exit(-1);}
		is>>pt.y;	 if(is.fail()) {cerr<<"Bad Input (not a double) \n"; exit(-1);}
		is>>in_char; if(in_char!=')'&& in_char!='}') {cerr<<"Bad Input: '"<<in_char<<"'\n"; exit(-1);}
		return is;
	}	
};

class contour{
private:
	vector<grid_pt> positions;
	double max_x, max_y, min_x, min_y;
public:
//Constructor
	contour(){max_x=-INFINITY; max_y=-INFINITY; min_x=INFINITY; min_y=INFINITY;}
//Destructor
	~contour(){}
//Mutators
	void add_point(const double x_pos, const double y_pos);
	void add_point(const grid_pt pt);
	void clear(){positions.clear(); max_x=-INFINITY; max_y=-INFINITY; min_x=INFINITY; min_y=INFINITY;}
//Members
	void span_of_contour(boost::unordered_set<pair<double,double> >& box_list,const grid_pt origin, const double s) const;
	double get_max_x()const{return max_x;}
	double get_max_y()const{return max_y;}
	double get_min_x()const{return min_x;}
	double get_min_y()const{return min_y;}
	int size()const{return positions.size();}
//Friends
	friend ostream& operator <<(ostream& out, const contour);
	friend istream& operator >> (istream& is, contour& cont);
};

//CONTOUR CLASS FUNCTIONS
void contour::add_point(const double x_pos, const double y_pos){
	positions.push_back(grid_pt (x_pos,y_pos));
	if (x_pos>max_x) max_x=x_pos; //add max x
	if (x_pos<min_x) min_x=x_pos; //add min x
	if (y_pos>max_y) max_y=y_pos; //add max y
	if (y_pos<min_y) min_y=y_pos; //add min y
}//adds a point to the contour and updates max and min 

void contour::add_point(const grid_pt pt){
	positions.push_back(pt);
	if (pt[1]>max_x) max_x=pt[1]; //add max x
	if (pt[1]<min_x) min_x=pt[1]; //add min x
	if (pt[2]>max_y) max_y=pt[2]; //add max y
	if (pt[2]<min_y) min_y=pt[2]; //add min y
}//adds a point to the contour and updates max and min 

ostream& operator <<(ostream& out, const contour cont){
	for (vector<grid_pt>::const_iterator it = (cont.positions).begin(), end = (cont.positions).end(); it != end;) {
		out << *it;
		if (++it != end) out <<","; //outputs each contour
		else out <<"\n";
		}
		return out;
} 

void contour::span_of_contour(boost::unordered_set<pair<double,double> >& box_list,const grid_pt origin, const double s) const{
	typedef grid_pt pt; //keeps the notation more concise
	for (vector<grid_pt>::const_iterator it=positions.begin(); it<prev(positions.end());++it){ //iterate points terminate one before end
				double x_spaces=((*next(it)).grid_posit(origin,s))[1]-((*it).grid_posit(origin,s))[1];
				double y_spaces=((*next(it)).grid_posit(origin,s))[2]-((*it).grid_posit(origin,s))[2];
				//if(abs(x_spaces)>50||abs(y_spaces)>50)
				//	cout<<(*it)<<*next(it)<<endl;
				if(abs(x_spaces)>0){ //if there is at least 1 x increment
				  	double m=((*next(it))[2]-(*it)[2])/((*next(it))[1]-(*it)[1]); //gradient
				  	double c=-m*((*it)[1])+(*it)[2]; //intercept
				  	//increment x
				  		for(double space=0; abs(space)<abs(x_spaces); (x_spaces > 0)?++space: --space)
				  			box_list.emplace((pt(s*floor((*it)[1]/s)+space,m*(s*floor((*it)[1]/s)+space)+c).grid_posit(origin,s)).as_pair()); //-> map to grid then add point (as pair)
				  	//increment y
				  		for(double space=0; abs(space)<abs(y_spaces); (y_spaces>0) ?++space: --space)
				  			box_list.emplace((pt((s*floor((*it)[2]/s)+space-c)/m,s*floor((*it)[2]/s)+space).grid_posit(origin,s)).as_pair()); //-> map to grid then add point (as pair)
				  	}//end if
				else // straight line upward/downward
					for(double space=0; abs(space)<=abs(y_spaces); (y_spaces > 0)?++space: --space)
				  		box_list.emplace((pt((*it)[1],(*it)[2]+space)).grid_posit(origin,s).as_pair()); //-> map to grid then add point (as pair)dd point
	}//end for loop
}//end span_of_contour

istream& operator >> (istream& is, contour& cont){
	char in_char; grid_pt pt;
	while(is.peek()!='\n'){
		is>>pt; cont.add_point(pt);
		if(is.peek()!='\n'&& is.peek()!=EOF) {
			is>>in_char; if(in_char!=','||is.fail()) {cerr<<"Bad Input: '"<<in_char<<"\n"; exit(-1);}
		}
		else break;
	}
	is.ignore();
	return is;
}	

}//end of namespace box_counting
#endif //BOX_COUNTER.H