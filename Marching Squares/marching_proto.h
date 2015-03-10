/*
 * marching_proto.cpp
 *
 *  Created on: Oct 31, 2012
 *      Original Author: Niels Walet
 *  Adapted on: Mar 10, 2015
 *		Adapted by: David Robinson
 */

#ifndef MARCHING_H_
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/binary.hpp>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>  // I/O
#include <fstream>  // I/O
#include <string>
#include <map>
#include <list>

enum directions {move_up, move_down, move_left, move_right, no_move};
enum start_dir {fwd, bwd, done};//forward and backward renamed so compatible with c++11
const int Nx=200;//510;//19;//900; //this is the size of the grid it seems
const int Ny=100;//190;//51;//900;

#define MARCHING_H_


#endif /* MARCHING_H_ */
