### Marching Squares Algorithm (Marching_Squares.h)
Adapted from code by Prof. Niels Walet

Usage: marching::draw_contours(data,contvalue,file_name). 

Arguments:

	* data: `boost::numeric::ublas::matrix<double>` contains 
"image" investigating

	* contvalue: `int` contour value to be investigated

	* file_name: `string` output contours filename

Notes:

	* Requires `-std=c++11` to compile (only line 315)

	* Output to mathematica module by default 
(testing purposes)

	* Still a work in progress 
