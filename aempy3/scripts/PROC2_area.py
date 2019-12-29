#!/usr/bin/env python2

"""
Created on Tue Jan  3 14:04:37 2017

@author: Duygu

Updated on Sun March 12 2017

PLOT_scatter and PLOT_interpolate are included...

"""
from __future__ import print_function, absolute_import, division
#import os
#import sys
#import time
#import math  as ma
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
#import random

import aemprocs as aem
import invprocs as inv


plt.close('all')


fmtx='%6i %1.8e %1.6e %5.1f %5.1f   %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f   %7.2f'
filename = 'StGormans_Rectangle'

plots_on = True
mycmap='magma'
#plotformat = '.jpg'data_obs= np.loadtxt(filename+'.dat', skiprows=1)
#
#plotformat = '.eps'
#plotformat = '.png'
#plotformat = '.pdf'
plotformat = '.eps'


Fsize = 16   # Axis label size 
Lsize = 14   # Axis tick label size
s0    = 10
m0 = 's'

#
#aem.process_aem05_data('GSI___15IRL_DLV1563_FEM2.xyz', 'A1', 'A1_FEM2')
#aem.choose_data_rect('A1_FEM2_All.dat','638968.67 641519.17 5922331.93 5924940.46','StGormans')

print (' Reading observed data from ' + filename )
data_obs= np.loadtxt(filename+'.dat', skiprows=1)

D=data_obs[:,:]
sizework=np.shape(D)
nvars=sizework[1]
last=nvars-1
print(' data block has shape: ',np.shape(D))
impute='delete'
#impute='interpolate'
# Flag bad Values
criterion = 'plm'
threshval=0.5
columns=[last,last]
print('\n criterion: '+criterion)
print(' columns: ', columns)
print(' thresh = ', threshval)
D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)
criterion = 'greater than'
threshval=70.
columns=[4,4]
nanincol=None
print('\n criterion: '+criterion)
print(' columns: ', columns)
print(' thresh = ', threshval)
D, nanindex = inv.prep_insert_flag(D, criterion,threshval,columns,incol=None)

columns=[4,12]
print('\n impute = ', impute)
print(' columns: ', columns)
D ,nanindex = inv.prep_handle_gaps(D, columns, impute)
np.savetxt(filename+'_preprocessed_'+impute+'_high.dat',D, delimiter=',',header=filename,fmt=fmtx)
print(' data block now has shape: ',np.shape(D))
#
#
#
I912                  = D[:,5]
I3005                 = D[:,6]
I11962                = D[:,7]
I24510                = D[:,8]
Q912                  = D[:,9]
Q3005                 = D[:,10]
Q11962                = D[:,11]
Q24510                = D[:,12]
    
IData                 = np.concatenate((I912, I3005, I11962, I24510))
QData                 = np.concatenate((Q912, Q3005, Q11962, Q24510))    
IData_min, IData_max  = np.min(IData), np.max(IData)
QData_min, QData_max  = np.min(QData), np.max(QData)
  

site_x                = D[:,1]  
site_y                = D[:,2]  
site_x                = site_x*0.001       #in km
site_y                = site_y*0.001       #in km
sitex_min, sitex_max  = np.min(site_x), np.max(site_x)
sitey_min, sitey_max  = np.min(site_y), np.max(site_y)

site_alt              = D[:,4]  
min_alt               = np.min(site_alt)
max_alt               = np.max(site_alt)
 


r=1



## Interpolated Data Map Plot ##
#  multiquadric': sqrt((r/self.epsilon)**2 + 1)
# 'inverse': 1.0/sqrt((r/self.epsilon)**2 + 1)
# 'gaussian': exp(-(r/self.epsilon)**2)
# 'linear': r
# 'cubic': r**3
# 'quintic': r**5
# 'thin_plate': r**2 * log(r)
# RBFfunction='thin_plate'
RBFfunction='inverse'
numIndexes            = 300
smooth = 0.


xi= np.linspace(sitex_min,sitex_max,numIndexes)    
yi= np.linspace(sitey_min,sitey_max,numIndexes)
XI, YI = np.meshgrid(xi, yi)
plt.close('all')
fI912                = interpolate.Rbf(site_x, site_y, I912, function=RBFfunction)
interpI912           = fI912(XI, YI)

fI3005               = interpolate.Rbf(site_x, site_y, I3005, function=RBFfunction)
interpI3005          = fI3005(XI, YI)

fI11962              = interpolate.Rbf(site_x, site_y, I11962, function=RBFfunction)
interpI11962         = fI11962(XI, YI)

fI24510              = interpolate.Rbf(site_x, site_y, I24510, function=RBFfunction)
interpI24510         = fI24510(XI, YI)



fQ912                = interpolate.Rbf(site_x, site_y, Q912, function=RBFfunction)
interpQ912           = fQ912(XI, YI)

fQ3005               = interpolate.Rbf(site_x, site_y, Q3005, function=RBFfunction)
interpQ3005          = fQ3005(XI, YI)

fQ11962              = interpolate.Rbf(site_x, site_y, Q11962, function=RBFfunction)
interpQ11962         = fQ11962(XI, YI)

fQ24510              = interpolate.Rbf(site_x, site_y, Q24510, function=RBFfunction)
interpQ24510         = fQ24510(XI, YI)



r=1 

if plots_on:
#==============================================================================
   #plt.figure(1)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10.9,6.5))
   cmap = plt.get_cmap('RdBu')
   im1= ax[0, 0].pcolor(XI, YI, interpI912, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('In-Phase - 912 Hz / '+RBFfunction,fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
#   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpI3005, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('In-Phase - 3005 Hz / '+RBFfunction,fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
#   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpI11962, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('In-Phase - 11962 Hz / '+RBFfunction,fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)', fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
#   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpI24510, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('In-Phase - 24510 Hz / '+RBFfunction)
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
#   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Iinterpolated_'+RBFfunction+plotformat) #dpi=1200)
   plt.show()
   plt.close()
   #==============================================================================
   
   #plt.figure(2)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10.9,6.5))
   cmap = mpl.cm.RdBu
   im1= ax[0, 0].pcolor(XI, YI, interpQ912, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('Quadrature - 912 Hz / '+RBFfunction,fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
#   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpQ3005, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('Quadrature - 3005 Hz / '+RBFfunction,fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
#   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpQ11962, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('Quadrature - 11962 Hz / '+RBFfunction,fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
#   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpQ24510, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('Quadrature - 24510 Hz / '+RBFfunction,fontsize=Fsize)
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
#   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Qinterpolated_'+RBFfunction+plotformat) #dpi=1200)
   plt.show()
   plt.close()
   1
   
k =2
print('\n PCA filter: ')
print(' N pca: ', k)
columns=[5,12]
print(' columns: ', columns)
D = inv.prep_pcafilter(D, k, columns)
np.savetxt(filename+'_Npca'+str(k)+'.dat',D, delimiter=',',header=filename+'_rewrite',fmt=fmtx)

#
I912                  = D[:,5]
I3005                 = D[:,6]
I11962                = D[:,7]
I24510                = D[:,8]
Q912                  = D[:,9]
Q3005                 = D[:,10]
Q11962                = D[:,11]
Q24510                = D[:,12]
    
IData                 = np.concatenate((I912, I3005, I11962, I24510))
QData                 = np.concatenate((Q912, Q3005, Q11962, Q24510))    
#IData_min, IData_max  = np.min(IData), np.max(IData)
#QData_min, QData_max  = np.min(QData), np.max(QData)
  

site_x                = D[:,1]  
site_y                = D[:,2]  
site_x                = site_x*0.001       #in km
site_y                = site_y*0.001       #in km
sitex_min, sitex_max  = np.min(site_x), np.max(site_x)
sitey_min, sitey_max  = np.min(site_y), np.max(site_y)

site_alt              = D[:,4]  
min_alt               = np.min(site_alt)
max_alt               = np.max(site_alt)


xi= np.linspace(sitex_min,sitex_max,numIndexes)    
yi= np.linspace(sitey_min,sitey_max,numIndexes)
XI, YI = np.meshgrid(xi, yi)
 
fI912                = interpolate.Rbf(site_x, site_y, I912, function=RBFfunction)
interpI912           = fI912(XI, YI)

fI3005               = interpolate.Rbf(site_x, site_y, I3005, function=RBFfunction)
interpI3005          = fI3005(XI, YI)

fI11962              = interpolate.Rbf(site_x, site_y, I11962, function=RBFfunction)
interpI11962         = fI11962(XI, YI)

fI24510              = interpolate.Rbf(site_x, site_y, I24510, function=RBFfunction)
interpI24510         = fI24510(XI, YI)



fQ912                = interpolate.Rbf(site_x, site_y, Q912, function=RBFfunction)
interpQ912           = fQ912(XI, YI)

fQ3005               = interpolate.Rbf(site_x, site_y, Q3005, function=RBFfunction)
interpQ3005          = fQ3005(XI, YI)

fQ11962              = interpolate.Rbf(site_x, site_y, Q11962, function=RBFfunction)
interpQ11962         = fQ11962(XI, YI)

fQ24510              = interpolate.Rbf(site_x, site_y, Q24510, function=RBFfunction)
interpQ24510         = fQ24510(XI, YI)



r=1 

if plots_on:
#==============================================================================
   #plt.figure(1)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10.9,6.5))
   cmap = plt.get_cmap('RdBu')
   im1= ax[0, 0].pcolor(XI, YI, interpI912, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('In-Phase - 912 Hz Pc'+str(k),fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
#   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpI3005, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('In-Phase - 3005 Hz Pc'+str(k),fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
#   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpI11962, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('In-Phase - 11962 Hz Pc'+str(k),fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)', fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
#   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpI24510, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('In-Phase - 24510 Hz Pc'+str(k),fontsize=Fsize)
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
#   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Iinterpolated_Npc'+str(k)+plotformat) #dpi=1200)
   plt.show()
   plt.close()
   #==============================================================================
   
   #plt.figure(2)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(10.9,6.5))
   cmap = mpl.cm.RdBu
   im1= ax[0, 0].pcolor(XI, YI, interpQ912, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('Quadrature - 912 Hz Pc'+str(k),fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
#   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpQ3005, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('Quadrature - 3005 Hz Pc'+str(k),fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
#   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpQ11962, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('Quadrature - 11962 Hz Pc'+str(k),fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
#   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpQ24510, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('Quadrature - 24510 Hz Pc'+str(k),fontsize=Fsize)
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
#   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Qinterpolated_Npc'+str(k)+plotformat) #dpi=1200)
   plt.show()
   plt.close()



# plt.close('all')
