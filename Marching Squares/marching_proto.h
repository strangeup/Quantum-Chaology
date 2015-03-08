/*
 * marching.h
 *
 *  Created on: Oct 31, 2012
 *      Author: niels
 */

#ifndef MARCHING_H_
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/binary.hpp>
#include <boost/foreach.hpp>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
//#include <fstream>   // file I/O
#include <iostream>  // I/O
#include <fstream>  // I/O
//#include <iomanip>   // format manipulation
#include <string>
#include <map>
#include <list>

enum directions {move_up, move_down, move_left, move_right, no_move};
enum start_dir {forward, backward, done};
const int nx=510;//19;//900; //this is the size of the grid it seems
const int ny=190;//51;//900;

#define MARCHING_H_


#endif /* MARCHING_H_ */
