### Marching Squares Algorithm

* Heavily based on code by Prof. Niels Walet

* Output to mathematica module by default (testing purposes)

* Note this is a work in progress 

### Problems with current Marching Algorithm

* Doesn't appear to be anything wrong with my implementation

* It throws the error because it "deletes" points from the map
 object then sometimes tries to access them again (which 
creates a key value pair initialised to zero)

* Doesn't appear to cause a problem for suitably spaced contours but
when contours are closely spaced (e.g. a figure eight) the code 
crashes

* An updated version with a crude fix avoids this error by catching 
the exception and "closing" the loop by placing it back at it's 
start value. This appears to be incorrect as a final incorrect 
point will be added to each contour. 
