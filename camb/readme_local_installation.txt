In order to install camb-ISiTGR and python-ISiTGR, you need to do the following:

* for camb-ISiTGR (fortran source code)
	1) go inside the fortran folder and modify the Makefile to set your compilers
	2) run "make" in your terminal (you might need or not a make clean first)
	3) inside the fortran folder you can execute "./camb ../inifiles/params_MG.ini" to test that the code is working properly
	
* for python-ISiTGR (the python wrapper for the isitgr code)
	1) first you need to install the fortran source code following the previous instructions
	2) in the principal directory (directory where you download the package), run "sudo python setup.py install" for a local installation 
	(do not forget to check that you have the dependencies)
	3) in the principal directory run "echo "export PYTHONPATH=$(pwd):\$PYTHONPATH" >> ~/.bashrc" to add the path of the isitgr module
	to the default python path of your system
	4) to make the changes effectively in your bashrc file, run "source ~/.bashrc"
	5) you can run the ISiTGR_reproducing_previous_results.ipynb notebook to test the wrapper.
