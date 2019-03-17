This is the public repository for the DIAS Airborne Electromagnetic Inversion 
toolbox "aempy"., which was developed starting with the project "Spatially 
constrained Bayesian inversion of frequency- and time-domain electromagnetic 
data from the Tellus projects", funded by the geological survey of Ireland GSI
(2015-sc-004), It is distributed under the GNU GENERAL PUBLIC LICENSE Version 3.
                       
Please keep in mind that this is an experimental software, and may contain 
errors. Use at your own risk! However, we will frequently update the repository 
correcting bugs, and adding additional functionality.                 
 
This repository contains the following subdirectories:

 - 	core
	This directory contains the Fortran 90 source code for the computational
	core run by the Python tollbox. currently it contains wrappers for the two
	systems used in Tellus: GTK4  and CGG Genesis. the numerics is derived from 
	the AMIRA/CSIRO AirBeo software. Wrappers for GEOTEM are work in progress.
	
 -	doc
 	Doumentation for the toolbox, and some useful documentation for python, 
 	including the most important extensions, numpy, scipy, and matplotlib 
 	
 -	modules
 	Contains the modules aemprocs.py and invprocs.py, which are called from the 
 	Python scripts run for different tasks of AEM inversion
 	
 - 	scripts
 	Contains the scripts  for preprocessing, visualization, and one-dimensional inversion of 
 	AEM data
 	
 - 	tutorial
 	Contains several scripts and data demonstratinmg explaining the typical 
 	work flow using the toolbox.      	 

 
Get your working copy via git from the command line:

 git clone https://git.dias.ie/vrath/aempy_public.git

Enjoy!

D. Kiyan & V. Rath
