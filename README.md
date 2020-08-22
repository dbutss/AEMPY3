# Tools for processing, inversion, and interpretation of Airborne ElectroMagnetics

This is the public repository for the  Airborne Electromagnetic Inversion toolbox "aempy" (python3 version), which was originally developed at DIAS starting with the project "Spatially constrained Bayesian inversion of frequency- and time-domain electromagnetic data from the Tellus projects", funded by the geological survey of Ireland GSI (2015-sc-004), It is distributed under the GNU GENERAL PUBLIC LICENSE Version 3.
                       
Please keep in mind that this is an experimental software, and may contain errors. Use at your own risk! However, we will frequently update the repository correcting bugs, and adding additional functionality.                 
 
This repository contains the following subdirectories:

 - 	**core1d**
	This directory contains the Fortran 90 source code for the computational
	core run by the Python toolbox. Currently it contains wrappers for the two
	systems used in Tellus: GTK4  and CGG Genesis. the numerics is derived from 
	the AMIRA/CSIRO AirBeo software. Wrappers for GEOTEM are work in progress.
	
 -	**info**
 	Doumentation for the toolbox, and some useful documentation for python, 
 	including the most important extensions, numpy, scipy, and matplotlib 
 	
 -	**modules**
 	Contains the modules aem.py, inv.py,  and mt.py which are called from the 
 	Python scripts run for different tasks of AEM inversion
 	
 - 	**scripts**
 	Contains the scripts  for preprocessing, visualization, and one-dimensional inversion of 
 	AEM data, explaining the typical work flow using the toolbox.      	 

Get your working copy via git from the command line:

_git clone https://github.com/volkerrath/AEMPY3/_

The scripts and jupyter notebooks are available in the subdirectory AEMPY3/aempy3. 

This version will run under Python up to 3.7 (3.8 not yet tested). To install it in an Linux environment (e.g. Ubuntu, SuSE), you need to the following:

(1) Download the latest Anaconda or Miniconda version (https://www.anaconda.com/distribution/), and install by runing the downloaded bash script. 

(2) Create an appropriate conda environment (including the necessary prerequisites) by:

_conda env create -f AEMpy3.yml_

(3) Activate this environment by:

_conda activate AEMPY3_



Enjoy!

D. Kiyan & V. Rath
