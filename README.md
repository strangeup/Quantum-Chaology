# Quantum-Chaology
A Finite Difference Approach to Quantum Chaology

#February 28

A complete version of the bindings (written by my project partner Paul Searle) and which I adapted for use as a seperate file was completed today. The previous protoptype worked well for small matrices but for large matrices encountered a segfault. This issue was due to the use of C-Style arrays (which we should have realised would cause a stack overflow for very large arrays). 

This problem was overcome by using boost::scoped_array data type, which functions very similarly but stores that data dynamically, not statically.

The other files for the simulation are currently a work in progress, they are fully functional but unfortunatlely rather verbose and in some places less than intuitive. A more general fashion will be worked on throughout the course of the project.
