#include <fstream>
#include <iostream>
#include "box_counter_proto.h"
#include <cstring>
using namespace std;
using namespace box_counting;


int main(int nNumberofArgs, char* pszArgs[]){
	//USER INTERFACE
	if(nNumberofArgs==1){
		cout<<"ARGUMENTS: max_box_size, filename1, filename2..\n";
		exit(-1);
	}
	//print a banner
        cout << "The arguments to " << pszArgs[0] <<  " are:\n";
    //now write out the remaining arguments
    for (int i = 1; i < nNumberofArgs; i++)
    	if (i!=1)
        	cout << "File: "<< i-1 << ": " << pszArgs[i] << "\n";
        else
        	cout<< "Maximum Box Size:"<<": "<<pszArgs[i] <<"\n";
     double max_box_size;
     stringstream ss;
     ss<<pszArgs[1];
     ss>>max_box_size;
     if(ss.fail()){cerr<<"Invalid argument for Max Box Size (double)\n"; exit(-1);}


    //FILE READER		
    for(int filenumber=2; filenumber<nNumberofArgs;++filenumber){
	    string filename=pszArgs[filenumber];
		vector<contour> contours;
		ifstream in; 

		in.open((filename+".dat").c_str()); //Open File
		if(in.fail()){cerr<<"Error Opening File: "<<filename+".dat"<<"\n";exit(-1);}
		contour current_contour;

		int line_count(0);
		 while(!in.eof()){ //read until end of file
		 		cout<<"Adding Contour: "<<line_count+1<<'\r';
		 		cout.flush(); ++line_count;
		 		current_contour.clear();
		 		in>>current_contour;
		 		contours.push_back(current_contour); current_contour.clear();
		 }
		in.close();

		double max_x,max_y,min_x,min_y;
		for(auto it=contours.begin();it<contours.end();++it){ //gets global max and min of contours
			if(max_x<(*it).get_max_x())
				max_x=(*it).get_max_x();
			if(max_y<(*it).get_max_y())
				max_y=(*it).get_max_y();
			if(min_x>(*it).get_min_x())
				min_x=(*it).get_min_x();
			if(min_y>(*it).get_min_y())
				min_y=(*it).get_min_y();
		}
		
		ofstream out;

		out.open((filename+"_boxes.dat").c_str());
		if(out.fail()){cerr<<"Error writing file: "<<filename+"_boxes.dat"<<"\n"; exit(-1);}

		typedef boost::unordered_set<pair<double,double>> my_set_t;
		my_set_t filled_boxes;
		double box_size(10);
		for (double box_size=1;box_size<max_box_size;++box_size){ //Loop through box_size
			filled_boxes=my_set_t();
			for(auto it=contours.begin();it<contours.end();++it){
				(*it).span_of_contour(filled_boxes,grid_pt(0,0),box_size);//add the boxes with contour sections to a set
			}
			out<<filled_boxes.size()<<"\t"<<box_size<<"\n";
		}
	}
	return 0;
}