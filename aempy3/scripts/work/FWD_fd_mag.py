# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by dkiyan - Sep 30

"""
from __future__ import print_function, absolute_import, division

import sys, os
import random
from datetime import datetime
import numpy as np

import aemprocs as aem
import invprocs as inv

import core1d

nan = np.nan #float('NaN')

random.seed(datetime.now())
#for test fix the random seed
#random.seed(2000)`

aem_system = 'aem05'
alt=61.
#sizedat=8
mode = 0
data_err = 30.
#data_err = 0.

identstr='_'+aem_system+'_3L_mag'

nlyr = 3
# resisitivity 
true_res = ([300.0, 10., 1000.0]) 
true_res=np.log10(true_res)
# relative magnetic permeability 
true_mu = ([1.0, 1.02, 1.0]) 
#  relative dielectric constant 
true_ep = ([1.0, 1. , 1.0]) 
# chargeability
true_m = ([0.0, 0.0, 0.0]) 
#time constant
true_t = ([0.0, 0.0, 0.0]) 
# freq c constant
true_c = ([1.0, 1., 1.0]) 

true_thk = ([25.0, 25.0])
true_depth = np.append( 0., np.cumsum(true_thk)) 


mactive, model_true, model_err , dum, dum = inv.initialize_1dmod(nlyr)

# set resistivity 
model_true[0*nlyr:1*nlyr] = true_res
# set relative magnetic permeability 
model_true[1*nlyr:2*nlyr] = true_mu
# set relative dielectric constant 
model_true[2*nlyr:3*nlyr] = true_ep
# set ip cole-cole parameter
model_true[3*nlyr:4*nlyr] = true_m
model_true[4*nlyr:5*nlyr] = true_t
model_true[5*nlyr:6*nlyr] = true_c
# set thickness
model_true[6*nlyr:7*nlyr-1]= true_thk[0:nlyr-1]

print \
(' Layer thicknesses: \n', true_thk)
print \
(' Layer interface depths: \n', true_depth)
print \
(' Parameters: \n', model_true)
print \
('\n ...calculating\n\n')



# calculate forward model
#data_true = core1d.aemfwd1d_aem05(mode,alt,nlyr,model_true)
data_true = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,model_true)')
data_obs= data_true;
print \
(' True response : \n', data_obs)


sizedat = np.shape(data_obs);
data_std= data_err*np.ones(sizedat)
data_cov= data_std*data_std
ndata=sizedat[0];
#  perturb data  
data_obs = data_true + data_err*np.random.randn(1, sizedat[0])
data_true = data_true + 0.*np.random.randn(1, sizedat[0])
print \
(' Perturbed  response : \n', data_obs)

np.savez_compressed('FWDModel'+identstr,mactive=mactive,nlyr=nlyr, model_true=model_true,
                    data_true=data_true,data_obs=data_obs, data_err=data_err,site_alt=alt)


