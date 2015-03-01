## Procedure for setting up [PyUblasExt](http://git.tiker.net/pyublasext.git/snapshot/a10af3278f2ebdf7d396644dafd681ab73000183.tar.gz)

* Download [arpack-ng](http://forge.scilab.org/index.php/p/arpack-ng/downloads/):

* Unzip using this command:

	/:~$ zcat arpack-ng_3.1.5.tar.gz | tar -xvf -D

* Now simply:
       
       /:~$ ./configure
       /:~$ make 
       /:~$ make install

*Download most recent [PyUblasExt](http://git.tiker.net/pyublasext.git/snapshot/a10af3278f2ebdf7d396644dafd681ab73000183.tar.gz) and unzip. Please note that there have been updates to aksetup so it's worthwhile downloading straight from the git.

* Replaced the default siteconf.py file with a custom version: [example](https://github.com/strangeup/Quantum-Chaology/blob/master/Inducer-Bindings-Tests/siteconf.py)
	* just find each of the libraries you want (eg libboost_python.a)
	* then for the name put boost_python and put its directory

* Type `:/~$ make install`
* If there are no errors type `:/~$ easy_install PyUblasExt`  
