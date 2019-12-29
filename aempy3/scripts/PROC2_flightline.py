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


plt.close('all')



random.seed(datetime.now())

plots_on = True
#plotformat = '.png'
plotformat = '.png'
#plotformat = '.svg'

Fsize      = 18
Lsize      = 15


fmtx='%6i %1.8e %1.6e %5.1f %5.1f   %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f   %7.2f'
impute='med'
#filename  = 'A1_FEM2_Sub_L1379'
filename  = 'bundoran_L7001_aempy'
print(' preprocessing file '+filename)

#filename  = 'A1_FEM2_Sub_L1377_preprocessed_median_nega'
#print(' preprocessing file '+filename)
data_obs= np.loadtxt(filename+'.dat', skiprows=1)
fline= data_obs[:,0]; 

D=data_obs[:,:]
sizework=np.shape(D)
nvars=sizework[1]
last=nvars-1
print(' data block has shape: ',np.shape(D))

######## Plot data along flight line #########

IData                 = np.concatenate((D[:,5], D[:,6], D[:,7], D[:,8]))
QData                 = np.concatenate((D[:,9], D[:,10], D[:,11], D[:,12]))    
IData_min, IData_max  = np.min(IData), np.max(IData)
QData_min, QData_max  = np.min(QData), np.max(QData)

Data_min              = np.min([IData_min, QData_min])
Data_max              = np.max([IData_max, QData_max])
   

site_x                = D[:,1]*0.001  #in km
site_y                = D[:,2]*0.001  #in km


if plots_on: 
   plt.figure(1)
   
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10.9,6.5))    
   ax[0, 0].plot(site_y[:],D[:,5],'r',site_y[:],D[:,9],'g')
   ax[0, 0].set_title('IP912 & Q912',fontsize=Fsize)
   ax[0, 0].set_ylim([Data_min,Data_max])
   ax[0, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   
   ax[0, 1].plot(site_y[:],D[:,6],'r',site_y[:],D[:,10],'g')
   ax[0, 1].set_title('IP3005 & Q3005',fontsize=Fsize)
   ax[0, 1].set_ylim([Data_min,Data_max])
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   
   ax[1, 0].plot(site_y[:],D[:,7],'r',site_y[:],D[:,11],'g')
   ax[1, 0].set_title('IP11962 & Q11962',fontsize=Fsize)
   ax[1, 0].set_ylim([Data_min,Data_max])
   ax[1, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[1, 0].set_xlabel('UTM Northing (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   
   ax[1, 1].plot(site_y[:],D[:,8],'r',site_y[:],D[:,12],'g')
   ax[1, 1].set_title('IP24510 & Q24510',fontsize=Fsize)
   ax[1, 1].set_ylim([Data_min,Data_max])
   ax[1, 1].set_xlabel(' UTM Northing (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   
   plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
   plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
   
   plt.legend(['In-Phase (IP)', 'Quadrature (Q)'],fontsize=Lsize, loc='upper left')
   fig = plt.savefig(filename+'-raw'+plotformat)
   plt.show()


   plt.figure(2)
   plt.plot(D[:,2]*0.001, D[:,4], '-r')
   plt.title('Flight Altitude',fontsize=Fsize)
   plt.xlabel('UTM Northing (km)', fontsize=Fsize)
   plt.ylabel('(m)', fontsize=Fsize)
   plt.grid(color='k',alpha=0.5, linestyle='dotted', linewidth=1.5) 
   plt.tick_params(labelsize=Lsize)
   plt.show

#criterion = 'plm'
#threshval=0.5
#columns=[last,last]
#impute='delete'
#print('\n criterion: '+criterion)
#print(' columns: ', columns)
#print(' thresh = ', threshval)
#D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)
#columns=[5,12]
#print(' impute = ', impute)
#print(' columns: ', columns)
#D ,nanindex = inv.prep_handle_gaps(D, columns, impute)
#np.savetxt(filename+'_preprocessed_'+impute+'_plm.dat',D, delimiter=',',header=filename+' preprocessed '+impute+'for power line monitor',fmt=fmtx)
#print(' data block now has shape: ',np.shape(D))
   
   

criterion = 'greater than'
threshval=100.
columns=[4,4]
nanincol=None
impute='delete'
print('\n criterion: '+criterion)
print(' columns: ', columns)
print(' thresh = ', threshval)
D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)
columns=[4,13]
print(' impute = ', impute)
print(' columns: ', columns)
D ,nanindex = inv.prep_handle_gaps(D, columns, impute)
np.savetxt(filename+'_preprocessed_'+impute+'_high.dat',D, delimiter=',',header=filename+' preprocessed '+impute+'for altitude',fmt=fmtx)
print(' data block now has shape: ',np.shape(D))




if plots_on: 
   plt.figure(1)
   
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(14.0,7.5))    
   
   ax[0, 0].plot(D[:,2]*0.001,D[:,5],'r',D[:,2]*0.001,D[:,9],'g')
   ax[0, 0].set_title('IP912 & Q912',fontsize=Fsize)
   ax[0, 0].set_ylim([Data_min,Data_max])
   ax[0, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   
   ax[0, 1].plot(D[:,2]*0.001,D[:,6],'r',D[:,2]*0.001,D[:,10],'g')
   ax[0, 1].set_title('IP3005 & Q3005',fontsize=Fsize)
   ax[0, 1].set_ylim([Data_min,Data_max])
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   
   ax[1, 0].plot(D[:,2]*0.001,D[:,7],'r',D[:,2]*0.001,D[:,11],'g')
   ax[1, 0].set_title('IP11962 & Q11962',fontsize=Fsize)
   ax[1, 0].set_ylim([Data_min,Data_max])
   ax[1, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[1, 0].set_xlabel('UTM Northing (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   
   ax[1, 1].plot(D[:,2]*0.001,D[:,8],'r',D[:,2]*0.001,D[:,12],'g')
   ax[1, 1].set_title('IP24510 & Q24510',fontsize=Fsize)
   ax[1, 1].set_ylim([Data_min,Data_max])
   ax[1, 1].set_xlabel(' UTM Northing (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   
   plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
   plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
   
   plt.legend(['In-Phase (IP)', 'Quadrature (Q)'],fontsize=Lsize, loc='upper left')
   fig = plt.savefig(filename+'-preprocessed'+plotformat)
   plt.show()
   
  
   

k = 3                    
columns=[5,12]
print('\n PCA filter: ')
print(' N pca: ', k)
Data_k = inv.prep_pcafilter(D, k, columns)
np.savetxt(filename+'_Npca'+str(k)+'.dat',Data_k, delimiter=',',header=filename+' Npca'+str(k),fmt=fmtx)


if plots_on: 
   plt.figure(1)
   
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(14.0,7.5))    
   ax[0, 0].plot(Data_k[:,2]*0.001,Data_k[:,5],'r',Data_k[:,2]*0.001,Data_k[:,9],'g')
   ax[0, 0].set_title('IP912 & Q912 Npca'+str(k),fontsize=Fsize)
   ax[0, 0].set_ylim([Data_min,Data_max])
   ax[0, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   
   ax[0, 1].plot(Data_k[:,2]*0.001,Data_k[:,6],'r',Data_k[:,2]*0.001,Data_k[:,10],'g')
   ax[0, 1].set_title('IP3005 & Q3005 Npca'+str(k),fontsize=Fsize)
   ax[0, 1].set_ylim([Data_min,Data_max])
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   
   ax[1, 0].plot(Data_k[:,2]*0.001,Data_k[:,7],'r',Data_k[:,2]*0.001,Data_k[:,11],'g')
   ax[1, 0].set_title('IP11962 & Q11962 Npca'+str(k),fontsize=Fsize)
   ax[1, 0].set_ylim([Data_min,Data_max])
   ax[1, 0].set_ylabel('(ppm)',fontsize=Fsize)
   ax[1, 0].set_xlabel('UTM Northing (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   
   ax[1, 1].plot(Data_k[:,2]*0.001,Data_k[:,8],'r',Data_k[:,2]*0.001,Data_k[:,12],'g')
   ax[1, 1].set_title('IP24510 & Q24510 Npca'+str(k),fontsize=Fsize)
   ax[1, 1].set_ylim([Data_min,Data_max])
   ax[1, 1].set_xlabel(' UTM Northing (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   
   plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
   plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
   
   plt.legend(['In-Phase (IP)', 'Quadrature (Q)'],fontsize=Lsize, loc='upper left')
   fig = plt.savefig(filename+'-PCAfilterNpca'+str(k)+plotformat)
   plt.show()

#plt.close('all')
