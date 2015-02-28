# PROCEDURE OF SETTING UP ARPACK.HPP
* First I set up arpack using the new generation version
* Downloaded and unzip using this:
	zcat arpack-ng_3.1.5.tar.gz | tar -xvf -D

* Then just $~ ./configure
	$~ make 
	$~ make install

* Then I replaced the default siteconf.py file with one I wrote myself
(see dropbox)
	-just find each of the libraries you want (eg libboost_python.a)
	-then for the name put boost_python and put its directory
* then just make then make install
	after this just type easy_install PyUblasExt  
