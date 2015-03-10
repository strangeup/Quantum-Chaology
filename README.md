# Quantum-Chaology
A Finite Difference Approach to Quantum Chaology. The simulation has been written for a masters project in Physics, but elements of it may be useful in other fields. If you are confused about any of the contents of this git it should at least be partially explained within this document, or in the readme files in each folder.

## Guide to Contents
This is intended to be a partial write-up of my findings and the progression of the simulation my project partner and I have written, with a few extra bits included that could be more generally useful.

Please note that amongst all of this I am still learning as well, any mistakes, contradictions or omissions please leave comments/issues! I don't profess that everything I have done is the best way (in fact I'm almost positive it isn't!) but this is how I have proceeded and I hope that this will be  useful for others starting out out with this manner of project.

The rest of this readme documents the progression of the project. For more detailed information redgarding specific files try my [wiki](https://github.com/strangeup/Quantum-Chaology/wiki) or the readme in the folder containing the file.

## Acknowledgements

Many thanks to Paul for co-writing this simulation with me and to [Inducer](https://github.com/inducer) for making the initial bindings used for LAPACK and UMFPACK so accessible.

## January 29
The correlation code was optimised as previously it was painfully slow. The change that was implemented to speed the system up highlighted a much more pressing problem that any assignment scaled at least as the dimension (D) squared of the matrix. This was obviously unnecessary as the matrix only contains ~5D (in principle known) elements. A more general approach to assignment was found and implemented which resulted in only necessary (non zero) elements being iterated through.  

## February 7

After initial attempts to install [Pyublas](http://mathema.tician.de/software/pyublas/) failed (using both apt-get and the make file) on the version of [xubuntu] (http://cdimage.ubuntu.com/ubuntu-core/releases/14.04/release/) I had installed due to an undetermined bug, a fresh install of the [latest xubuntu] (http://releases.ubuntu.com/14.10/) was performed.

[LAPACK](http://www.netlib.org/lapack/), [BLAS](http://www.netlib.org/blas/), [BOOST](http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html) and [ARPACK-ng](http://forge.scilab.org/index.php/p/arpack-ng/ ) were installed and set-up (with [bindings](https://svn.boost.org/svn/boost/sandbox/numeric_bindings/) for the former) prior to PyUblasExt.

Then PyUBlasExt was installed, using instructions that can be found [here](https://github.com/strangeup/Quantum-Chaology/blob/master/Inducer-Bindings-Tests/procedure.txt).

## February 15

Now that both [Pyublas](http://mathema.tician.de/software/pyublas/) and [PyUblasExt](https://pypi.python.org/pypi/PyUblasExt/0.92.4) were correctly installed investigation of [arpack.hpp] (https://github.com/inducer/arpack/blob/master/CPLUSPLUS/arpack.hpp) continued in more detail. The bindings were thoroughly commented to allow the understanding of the programs structure, the function 
perform\_reverse\_communication's arguments and it's capabilities. 

A minimal example was attempted, which compiled but did not produce any sensible results. It was then realised that this manner of approach may not work for a small matrix so it was implemented in the main body to attempt to try a larger matrix. Again this produced no physically acceptable results. 


## February 24

After failing to find a way of working the ARPACK bindings we initially set out to use (see Inducers bindings [here]( https://github.com/inducer/arpack/blob/master/CPLUSPLUS/arpack.hpp)). We decided to try our hand at our own. These bindings lack the finesse and generality of the ones we orignally attempted to use but using the others we did not manage to get reasonable results (they were either infinite or very close to zero).

Unfortunately a lot of time was wasted trying to get Inducers bindings to work, as the installation of PyUblas and PyUblas-ext did not go smoothly. The code itself was also difficult to work out with no usable documentation and little explanation. It also does not appear to be finished (none of the symmetric routines had the C macros set up for a function call in [arpack_proto](https://github.com/inducer/arpack/blob/master/CPLUSPLUS/arpack_proto.hpp), which makes me question whether usable results could have been obtained from it. Hopefully I will get time to write a symmetric macro that will similarly mimic overloading of arguments in c and then use that to implement my own version of bindings for this software. At this stage however I will write the bindings fit for purpose and (hopefully) add a few more before the end.


## February 28

A complete version of the bindings (written by my project partner Paul) which I adapted for use as a seperate file was completed today. The previous protoptype worked well for small matrices but for large matrices encountered a segfault. This issue was due to the use of C-Style arrays (which we should have realised would cause a stack overflow for very large arrays). 

This problem was overcome by using [boost scoped_array](http://www.boost.org/doc/libs/1_57_0/libs/smart_ptr/scoped_array.htm) data type, which functions very similarly but stores data dynamically, not statically.

The other files for the simulation are currently a work in progress, they are fully functional but unfortunately rather verbose and in some places less than intuitive. A more general fashion will be worked on throughout the course of the project.

## March 8

Several minor changes were implemented in the main code, including a (very) basic interface to allow a multi use static compilation. Investigation of the use of [Intel MKL Pardiso](https://software.intel.com/en-us/node/521677) instead of [umfpack](http://faculty.cse.tamu.edu/davis/suitesparse.html) has been made, attempting to solve a problem with infinite eigenvalues for extremely large matrices. Another project underway is the implementation of a [marching squares algorithm] (http://en.wikipedia.org/wiki/Marching_squares). This algorithm would allow the deduction of the contours of nodal surfaces in the wavefunctions of billiards, which would allow a very sparse storage and extra resolution using averaging.

## March 9

Marching squares algorithm fixed. This was due to a mismatch in the labelling of Nx and Ny.
