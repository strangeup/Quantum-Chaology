# Quantum-Chaology
A Finite Difference Approach to Quantum Chaology

#February 24

After failing to find a way of working the ARPACK bindings we initially set out to use (see Inducers bindings [here]( https://github.com/inducer/arpack/blob/master/CPLUSPLUS/arpack.hpp])). We decided to try our hand at our own. These bindings lack the finesse and generality of the ones we orignally attempted to use but using the others we did not manage to get reasonable results (they were either infinite or very close to zero).

Unfortunately a lot of time was wasted trying to get Inducers bindings to work, as the installion of PyUblas and PyUblas-ext did not go smoothly. The code itself was also difficult to work out with no usable documentation and little explanation. It also does not appear to be finish (none of the symmetric routines had the C macros set up for a function call in ARPACK_PROTO), which makes me question whether usable results could have been obtained from it. Hopefully I will get time to write a symmetric macro that will mic overloading of arguments in c and then use that to implement my own version of bindings for this software. At this stage however I will write the bindings fit for purpose and (hopefully) add a few more before the end.


#February 28

A complete version of the bindings (written by my project partner Paul) which I adapted for use as a seperate file was completed today. The previous protoptype worked well for small matrices but for large matrices encountered a segfault. This issue was due to the use of C-Style arrays (which we should have realised would cause a stack overflow for very large arrays). 

This problem was overcome by using boost::scoped_array data type, which functions very similarly but stores that data dynamically, not statically.

The other files for the simulation are currently a work in progress, they are fully functional but unfortunatlely rather verbose and in some places less than intuitive. A more general fashion will be worked on throughout the course of the project.
