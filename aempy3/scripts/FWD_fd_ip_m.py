# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 10:57:01 2016

@author: vrath

edited by dkiyan - Sep 30

"""
#from __future__ import print_function, absolute_import, division

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
alt=60.
#sizedat=8
mode = 0

data_err = 0.

identstr='_'+aem_system+'_L2D10_m'
nlyr = 3
true_res = ([300.0,30., 1000.0])  #resisitivity p.shape()
#true_res=np.log10(true_res)   # relative magnetic permeability 
true_mu = ([1.0, 1., 1.0]) # relative magnetic permeability 
true_ep = ([1.0, 1. , 1.0]) #  relative dielectric constant 
true_m = ([0.0, 0.8, 0.0]) # chargeability
true_t = ([0.0, .1, 0.0]) #time constant
true_c = ([1.0, 0.6, 1.0]) # freq c constant
true_thk = ([25.0, 10.0])

true_depth = np.append( 0., np.cumsum(true_thk)) 

mactive, model_true, model_err, dum, dum = inv.initialize_1dmod(nlyr)
# set resistivity 
model_true[0*nlyr:1*nlyr] = true_res # set resistivity 
model_true[1*nlyr:2*nlyr] = true_mu # set relative magnetic permeability 
model_true[2*nlyr:3*nlyr] = true_ep # set relative dielectric constant 
model_true[3*nlyr:4*nlyr] = true_m # set ip cole-cole parameters
model_true[4*nlyr:5*nlyr] = true_t
model_true[5*nlyr:6*nlyr] = true_c
model_true[6*nlyr:7*nlyr-1]= true_thk[0:nlyr-1] # set thickness

print \
(' Layer thicknesses: \n', true_thk)
print \
(' Layer interface depths: \n', true_depth)
print \
(' Parameters: \n', model_true)
print \
('\n ...calculating\n\n')
# calculate forward models
test_m = ([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]) 
size_m = np.shape(test_m)

model = model_true
data_obs = nan*np.ones((size_m[0],8))

for ii in range(size_m[0]):
    model[3*nlyr+1]=test_m[ii]
    # print(model[3*nlyr:4*nlyr])
    data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,model)')
    sizedat = np.shape(data_calc);
    ndata=sizedat[0];
    #  perturb data  
    data_calc = data_calc + data_err*np.random.randn(1, sizedat[0])
    print \
    (' Response : \n', data_calc)

    data_obs[ii,:]=data_calc
    
sizedat = np.shape(data_obs);   
data_std= data_err*np.ones(sizedat)
data_cov= data_std*data_std
print \
('\n Data to : FWDModel'+identstr+'\n')

np.savez_compressed('FWDModel'+identstr,mactive=mactive,nlyr=nlyr, model_true=model_true,
                    data_calc=data_calc,data_obs=data_obs, data_err=data_err,site_alt=alt)

test_t = ([0.1]) 
test_c = ([0.6]) 

