#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 14:39:17 2017

@author: vrath

Edited on Sat March 11 2017 - plots are added

@author: duygu



"""
from __future__ import print_function, absolute_import, division

import sys
import os
import time
import math  as ma
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import random
from datetime import datetime



import aemprocs as aem
import invprocs as inv



FileName='TellusA1'
random.seed(datetime.now())

plots_on = False
#plotformat = '.png'
plotformat = '.png'
#plotformat = '.svg'

Fsize      = 18
Lsize      = 15

SVals = np.empty((0,8))
MVals = np.empty((0,8))

filelist=inv.get_filelist(mystr='A1_FEM2_Sub_L????.dat',mypath='.')
    
nfiles =np.size(filelist)
print(nfiles)
#SVals = np.nan*np.ones((1,8))
#MVals = np.nan*np.ones((1,8))


start = time.time()
nsites = 0
for filename in filelist:
   nfiles = nfiles+1
   print('\n preprocessin0g file '+filename)
   data_obs= np.loadtxt(filename, skiprows=1)
   fline= data_obs[:,0]; 
   #
   D=data_obs[:,:]
   sizework=np.shape(D)
   nvars=sizework[1]
   last=nvars-1
   print(' data block has shape: ',np.shape(D))
   #
   criterion = 'plm'
   threshval=0.5
   columns=[last,last]
   impute='delete'
#   print('\n criterion: '+criterion)
#   print(' columns: ', columns)
#   print(' thresh = ', threshval)
   D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)
   columns=[4,12]
#   print(' impute = ', impute)
#   print(' columns: ', columns)
   D ,nanindex = inv.prep_handle_gaps(D, columns, impute)
   #np.savetxt(filename+'_preprocessed_'+impute+'_plm.dat',D, delimiter=',',header=filename+' preprocessed '+impute+'for power line monitor',fmt=fmtx)
   # print(' data block now has shape: ',np.shape(D))
   #
   criterion = 'greater than'
   threshval=70.
   columns=[4,4]
   nanincol=None
   impute='delete'
#   print('\n criterion: '+criterion)
#   print(' columns: ', columns)
#   print(' thresh = ', threshval)
   D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)
   columns=[4,13]
#   print(' impute = ', impute)
#   print(' columns: ', columns)
   D ,nanindex = inv.prep_handle_gaps(D, columns, impute)
   #np.savetxt(filename+'_preprocessed_'+impute+'_high.dat',D, delimiter=',',header=filename+' preprocessed '+impute+'for power line monitor',fmt=fmtx)
   print(' data block now has shape: ',np.shape(D))
   
   nDfinal=np.shape(D)
   nsites = nsites+nDfinal[0]
   #   
   print(' running pca filter ')
   columns=[5,12]
   ncols= np.size(range(columns[0],columns[1]+1))
   M = np.zeros(8)
   k=0
   while k < ncols:
      k= k+1
      # print(' N pca: ', k)
      Data_k, U, S, V, MSE = inv.prep_pcafilter(D, k, columns, out_full=True)
      #ident = filename+'_Npca'+str(k)
      # np.savez_compressed(filename+ident, k=k, Data=Data_k, S=S)
      M[k-1] = MSE
      # print(Merr)
   
   S = S/S[0]
   SVals = np.append(SVals,[S.T], axis=0)
   MVals = np.append(MVals,[M  ], axis=0)
   
   
elapsed = (time.time() - start)
print (' Used %7.4f sec for %6i lines  - %8i sites\n' % (elapsed, nfiles, nsites))
   #
   
   
x = range(1,8+1)
y = SVals

plt.figure(1,figsize=(10.9,6.5))
plt.plot(x,y.T,'r',linewidth=.1)
plt.title('Singular Values (pcafilter)',fontsize=Fsize)
plt.xlabel(' #SV ',fontsize=Fsize)
plt.ylabel(' SV (-)',fontsize=Fsize)
plt.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
plt.tick_params(labelsize=Lsize)

y = MVals
plt.figure(2,figsize=(10.9,6.5))
plt.plot(x,y.T,'r',linewidth=.1)
plt.title('MSE (pcafilter)',fontsize=Fsize)
plt.xlabel(' #SV ',fontsize=Fsize)
plt.ylabel(' MSE (ppm)',fontsize=Fsize)
plt.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
plt.tick_params(labelsize=Lsize)
F=FileName+'_RMSE'+plotformat
print (' PLot to ',
       F)
fig = plt.savefig(F)