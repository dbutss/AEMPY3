#! /usr/bin/python
"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by dkiyan - Sep 30
edited by vrath  - Aug 13, 2018
"""

#from __future__ import print_function, absolute_import, division

import sys
import os
import random
from datetime import datetime

import numpy as np

import invprocs as inv

import core1d

nan = np.nan # float('NaN')

random.seed(datetime.now())
#for test fix the random seed
#random.seed(2000)`

aem_system = 'genesis'
print ('\n  AEM forward modeling for '+ aem_system+' system\n \n') 


alt=100.
mode = 1
data_err_rel = 0.05
data_err_abs = 3000.


nlyr = 3
true_res = ([100.0, 1000.0, 100.0]) # case A
note='_'+aem_system+\
'_LR-'+str(true_res[0])+'-'+str(true_res[1])+'-'+str(true_res[2])+\
'_E-'+str(data_err_rel)+'-'+str(data_err_abs)


#true_res=np.log10(true_res)
true_mu = ([1.0, 1.0, 1.0])     # relative magnetic permeability 
true_ep = ([1.0, 1. , 1.0])     #  relative dielectric constant 
true_m = ([0.0, 0.0, 0.0])      # chargeability
true_t = ([0.0, 0.0, 0.0])      #time constant
true_c = ([1.0, 1.0, 1.0])      # freq c constant
true_thk = ([50.0, 25.0])
true_depth = np.append( 0., np.cumsum(true_thk)) 


mactive, model_true, model_err, dum, dum = inv.initialize_1dmod(nlyr)


# fill model vector
model_true[0*nlyr:1*nlyr] = true_res
model_true[1*nlyr:2*nlyr] = true_mu
model_true[2*nlyr:3*nlyr] = true_ep
model_true[3*nlyr:4*nlyr] = true_m
model_true[4*nlyr:5*nlyr] = true_t
model_true[5*nlyr:6*nlyr] = true_c
model_true[6*nlyr:7*nlyr-1]= true_thk[0:nlyr-1]


print (' Layer thicknesses: \n', true_thk)
print (' Layer interface depths: \n', true_depth)
print (' Parameter: \n', model_true)
print ('\n ...calculating\n\n')



# calculate forward model
# data_true = core1d.aemfwd1d_genesis(mode,alt,nlyr,model_true)
data_true = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,model_true)')
data_obs= data_true.copy();
print (' True response : \n', data_obs)



sizedat = np.shape(data_obs);
data_err= data_err_rel*data_obs+data_err_abs
data_cov= data_err*data_err
ndata=sizedat[0];

#  perturb data  
data_obs = data_true.copy() + data_err*np.random.randn(1, sizedat[0])
print (' Perturbed  response : \n', data_obs)


Outfile = 'FWDModel'+note
print('\n\n Output to : \n',Outfile)
np.savez_compressed(Outfile,
                    mactive=mactive, 
                    model_true=model_true, 
                    data_true=data_true,
                    data_obs=data_obs,
                    data_err=data_err)
