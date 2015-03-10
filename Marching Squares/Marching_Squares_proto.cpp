/*
 * marching_proto.cpp
 *
 *  Created on: Oct 31, 2012
 *      Original Author: Niels Walet
 *  Adapted on: Mar 10, 2015
 *		Adapted by: David Robinson
 */


#include "marching_proto.h"

using namespace std;

boost::numeric::ublas::matrix<double> data(Ny,Nx); //global stuff
boost::numeric::ublas::matrix<bool> mask(Ny,Nx);
double contvalue(0);	

class pos //just a container for grid positions
{
public:
	pos(const int ix,const int iy) {x=ix;y=iy;};
	pos() {x=-1;y=-1;};
	void decx() {x--;};	void incx() {x++;};
	void decy() {y--;};	void incy() {y++;};
	friend bool operator <(const pos pos1,const pos pos2)
	{ return ((pos1.x<pos2.x)|((pos1.x==pos2.x)&&(pos1.y<pos2.y)));//x1<x2 bitwiseOR x1==x2 && y1<y2
	};
	friend bool operator !=(const pos pos1,const pos pos2)
	{ return ((pos1.x!=pos2.x)||(pos1.y!=pos2.y)); //x1!=x2 || or y1!=y2
	};
	friend bool operator ==(const pos pos1,const pos pos2)
	{ return ((pos1.x==pos2.x)&&(pos1.y==pos2.y)); //x1==x2 && y1==y2
	};
	bool inside ()
	{ return (0<=x &&x<Ny-1&&0<=y&&y<Nx-1);}; //0<x<Ny-1 && 0<=y<Nx-1
	friend std::ostream& operator<<(std::ostream& os, const pos& posi)
	{ os <<"("<<posi.x<<","<<posi.y<<")";return os;}; //overload << as (x,y)
	int xc() const {return x;}; int yc() const {return y;}; //xcoord ycoord
private:
	int x, y;
};

class rpos //container for interpolated positions between grid points
{
public:
	rpos(const double xi,const double yi) {x=xi;y=yi;}; //constructor
	rpos() {x=-100;y=-100;}; //default
	rpos(const pos& lowercorner, const directions& WhichDir) //constructor
	{//interpolates
		int x1,y1,x2,y2;
		switch (WhichDir) {
		case move_left:
			x1=lowercorner.xc();y1=lowercorner.yc();x2=x1;y2=y1+1;break;
		case move_right:
			x1=lowercorner.xc()+1;y1=lowercorner.yc();x2=x1;y2=y1+1;break;
		case move_down:
			x1=lowercorner.xc();y1=lowercorner.yc();x2=x1+1;y2=y1;break;
		case move_up:
			x1=lowercorner.xc();y1=lowercorner.yc()+1;x2=x1+1;y2=y1;break;
		case no_move:
			cout <<" no_move reached 1";  exit(-1);
		}
		double d1=data(x1,y1)-contvalue,d2=data(x2,y2)-contvalue; //data is input matrix (global const) contvalue I presume is the value of the contour we are interested in
		x=(-d2*x1+d1*x2)/(d1-d2); y=(-d2*y1+d1*y2)/(d1-d2); //interpolations
	};
	rpos& operator +=(const rpos& add){
		x+=add.x;y+=add.y;return *this; //overload +=
	}
	rpos& operator /=(const int& div){
		x=x/div;y=y/div;return *this; //overload /=
	}
	friend std::ostream& operator<<(std::ostream& os, const rpos& posi)
	{ os <<"{"<<posi.x<<","<<posi.y<<"}";return os;}; //overload <<
	friend double dist(const rpos& pos1,const rpos& pos2){ //magnitude of rpos
		return sqrt(pow(pos1.x-pos2.x,2)+pow(pos1.x-pos2.x,2));
	}
	double xc() const {return x;}; double yc() const {return y;};//get x and y

private:
	double x, y;
};


inline void take_step(const directions dir, pos& posi, pos& posn, //increments x and y etc
		list<rpos>& curr_contour, directions& step, const start_dir& startdir)
{
	switch (dir) {
	case move_right:posn.incx();break; //increment x
	case move_left:posn.decx();break; //decrement x
	case move_up:posn.incy();break; //inc y
	case move_down:posn.decy();break;//dec y
	case no_move: exit(-1); //no_move exit with error
	}
	step=dir;
	rpos interpos=rpos(posi,dir);
	if (startdir==fwd) curr_contour.push_back(interpos);
	else curr_contour.push_front(interpos);
};

std::ostream& operator<<(std::ostream& os, const list<rpos>& contour) { //overloaded output
	  os <<"Line[{\n";
	  for (list<rpos>::const_iterator iterator = contour.begin(), end = contour.end();
			  iterator != end;) {
			  os << *iterator;
			  if (++iterator != end) os <<",\n"; //outputs each contour
			  else os <<"}]\n";
		  }
	  return os;
}

void draw_contours() { //this is the "main" function
	map<pos, int>  points;
	ofstream maskout;
	maskout.open("mask.csv");
	if (maskout.fail()){cerr<<"Could not open file\n"; exit(-1);}

	for(int ix=0;ix<Ny;ix++)
		for(int iy=0;iy<Nx;iy++)
		{	
			mask(ix,iy)=data(ix,iy)<contvalue; //mask is a bool. this is the mask as described on wiki
			maskout << mask(ix,iy)<< (iy==Nx-1 ? '\n':',');
		}
	maskout.close();



	for(int ix=0;ix<Ny-1;ix++)
		for(int iy=0;iy<Nx-1;iy++) {
			int label=2*(2*(2*mask(ix,iy+1)+mask(ix+1,iy+1))+
						mask(ix+1,iy))+mask(ix,iy); //convert to a bool between 0 and 16

			if (label == BOOST_BINARY( 1010 ) || label == BOOST_BINARY( 0101 ) ) { //if 10 or 5 => saddle points
				double mid = 0.25*(data(ix,iy)+data(ix+1,iy+1)+data(ix+1,iy)+data(ix,iy+1));//f(ix+0.5),ymap(iy+0.5));
				label+=16*(mid<contvalue); //gives 26 and 21
			}
			if (label != BOOST_BINARY( 0000 ) && label !=BOOST_BINARY( 1111 ) ) { //don't add trivial cases 0 and 15
				points[pos(ix,iy)]=label; //add to points
			//	cout <<ix<<" "<<iy<<" : "<<label<<" ";
			}
	}
	list<list<rpos> > contours;
	map<pos, int>::iterator it;

	while(points.begin()!=points.end()){ //while not on last value
		it=points.begin();
	list<rpos> curr_contour; //current contour
	pos posi, posn; //position initial and end

	start_dir firststep=fwd; //enum: fwd bwd or done

	do {
		posi=it->first;
		directions step=no_move;
		while (posi.inside()&& (posi!=it->first||step==no_move)){
			posn=posi;
//			cout << "LOOKING AT "<<posi<<" : "<<points[posi]<<"\n";
			int current_point_type=points[posi];
			switch(current_point_type)
			{

	/* 11
	   01 : 14 or 1, x-1 or y+1 */
				case BOOST_BINARY(1110): case BOOST_BINARY(0001): //case 1 or 14
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ( (step == no_move && firststep == fwd) |
						(step != no_move && step != move_right))
					take_step(move_left, posi, posn, curr_contour,step, firststep);
				else
					take_step(move_down, posi, posn, curr_contour,step, firststep);
				break;
	/* 11
	   10 : 13 or 2: x+1 or y-1*/
				case BOOST_BINARY(1101):case BOOST_BINARY(0010): //case 2 or 13
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ((step == no_move && firststep == fwd) |
						(step != no_move &&step != move_left))
					take_step(move_right, posi, posn, curr_contour,step, firststep);
				else
					take_step(move_down,posi, posn, curr_contour,step, firststep);

				break;
		/* 10
		   11 : 11 or 4: x+1 or y+1*/
				case BOOST_BINARY(1011):case BOOST_BINARY(0100): //case 4 or 11
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ((step == no_move && firststep == fwd) |
					(step != no_move &&step!=move_down))
					take_step(move_up,posi, posn, curr_contour,step, firststep);
				else
					take_step(move_right,posi, posn, curr_contour,step, firststep);

				break;
		/* 01
		   11 : 7 or 8: x-1 or y+1*/
				case BOOST_BINARY(0111):case BOOST_BINARY(1000): //case 8 or 7
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ((step == no_move && firststep == fwd) |
					(step != no_move &&step != move_right))
					take_step(move_left,posi, posn, curr_contour,step, firststep);
				else
					take_step(move_up,posi, posn, curr_contour,step, firststep);

				break;
		/* 11
		   00 : 12 or 3: x-1 or x+1*/
				case BOOST_BINARY(1100):case BOOST_BINARY(0011): //case 3 or  12
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ((step == no_move && firststep == fwd) |
					(step != no_move &&step != move_right))
					take_step(move_left,posi, posn, curr_contour,step, firststep);
				else
					take_step(move_right,posi, posn, curr_contour,step, firststep);
				break;
		/* 10
		   10 : 9 or 6 y-1 or y+1*/
				case BOOST_BINARY(1001):case BOOST_BINARY(0110): // case 9 or 6
				if ( !(step == no_move && firststep == fwd)) points.erase(posi);
				if ((step == no_move && firststep == fwd) |
					(step != no_move &&step != move_up))
					take_step(move_down,posi, posn, curr_contour,step, firststep);
				else
					take_step(move_up,posi, posn, curr_contour,step, firststep);
				break;
		/* 10
		   01 : 26 or 5 : double, diagonals sloping move_down*/
				case BOOST_BINARY(11010):case BOOST_BINARY(00101): //case 26 or 5
					switch (step) {
					case move_down:
						take_step(move_right,posi, posn, curr_contour,step, firststep);
						points[posi]= BOOST_BINARY(0001); break;
					case move_up:
						take_step(move_left,posi, posn, curr_contour,step, firststep);
						points[posi]= BOOST_BINARY(0100);break;
					case move_left:
						take_step(move_up,posi, posn, curr_contour,step, firststep);
						points[posi]= BOOST_BINARY(0001); break;
					case move_right:
						take_step(move_down,posi, posn, curr_contour,step, firststep);
						points[posi]= BOOST_BINARY(0100); break;
					case no_move:  /* started at degenerate point */
						if (firststep == fwd)
							take_step(move_right,posi, posn, curr_contour,step, firststep);
						/* don't delete half of point yet */
						else
							{take_step(move_up,posi, posn, curr_contour,step, firststep);
							points[posi]= BOOST_BINARY(0100);}
				};
				break;
		/* 01
		   10 : 21 or 10: double, diagonals sloping up*/ //case 21 or 10
		case BOOST_BINARY(10101):case BOOST_BINARY(01010):
			switch (step) {
			case move_down:
				take_step(move_left,posi, posn, curr_contour,step, firststep);
				points[posi]= BOOST_BINARY(0010); break;
			case move_up:
				take_step(move_right,posi, posn, curr_contour,step, firststep);
				points[posi]= BOOST_BINARY(1000); break;
			case move_left:
				take_step(move_down,posi, posn, curr_contour,step, firststep);
				points[posi]= BOOST_BINARY(0010); break;
			case move_right:
				take_step(move_up,posi, posn, curr_contour,step, firststep);
				points[posi]= BOOST_BINARY(1000);  break;
			case no_move:
				if (firststep == fwd)
					take_step(move_left,posi, posn, curr_contour,step, firststep);
				/* don't delete half of point yet */
				else
					{take_step(move_down,posi, posn, curr_contour,step, firststep);
					points[posi]= BOOST_BINARY(0010);}
		};
		break ;
		default:
			cout <<"this should have been caught before? "<<current_point_type<<"\n"<<"position: "<< posn<<"start: "<<(it->first); exit(-1);
	}//switch
//	cout <<"increment\n";
	posi=posn;
	//cout <<posn<<"\n";
	}//while inner
//	cout << "reached end"<<posi<<"\n";
	if (posi==(it->first)) {
//		cout <<"closed loop \n";
		firststep=done;// closed loop, finished
		points.erase(posi);
	}
	else if (firststep==fwd) firststep=bwd;
	else firststep=done;
	} while (firststep!=done); //while do
	contours.push_back(curr_contour);
	} //while loop
	ofstream out;
	out.open("data-test.m");
	if (out.fail()){cerr<<"Could not open file\n"; exit(-1);}
	out<<"Graphics[{\n";
	for (list<list<rpos> >::const_iterator it=contours.begin(),end=contours.end(); it!=end;++it)
			out<<(*it)<<(end==std::next(it) ? "}]":",");

	out.close();
}

int main(){
	fstream in;
	in.open("4.dat");
	if (in.fail()){cerr<<"Could not open file"; exit(-1);}
	for (int i=0; i<Ny; ++i)
		for (int j=0; j<Nx ; ++j)
			in >> data(i,j);
	in.close();
	draw_contours();
	return 0;
}