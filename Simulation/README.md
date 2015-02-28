#Contents of Folder
This simulates 2D quantum billiards using a finite difference approach. For more information about theory see [my report](https://github.com/strangeup/Quantum-Chaology/blob/master/Theory/Project_1.pdf).

To compile (set-up dependant) g++ 2Dmain -o name -std=c++11 -llapack -larpack -lumfpack -lamd -I /path/to/suiteparse/ 

* 2Dmain.cpp : The function that calls everything three main approachs spectral analysis: for analysis (using Random Matrix Theory) of the eigenvalue spectrum
    * direct approach, uses Fortan77 LAPACK routine _syev()
  
* correlation function: an attempt to characterise complexity of spectrum using a matrix correlation approach. This is ultimately not completely tested and may contain slight errors
  
* arpack function: sets up large sparse matrices and solves using the purpose made arpack bindings. Uses the [symmetric shift invert bindings](https://github.com/strangeup/Quantum-Chaology/blob/master/Bindings/arpack_symm_si.hpp)

* 2Dheaders.h : external headers, global settings and input parameters

* 2Dfunctions.h : all of the functions created for use in model (unfinished in places)
