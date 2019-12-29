#!/usr/bin/env python2

"""
Created on Tue Jan  3 14:04:37 2017

@author: Duygu

Updated on Sun March 12 2017

PLOT_scatter and PLOT_interpolate are included...

"""

import os
import sys
import time
import math  as ma
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


filename = 'StGormans_allLines'
blocksize = 10


import aemprocs as aem
import invprocs as inv

import core1d



"""
First step is to transform AEM05 Data to AEMPY-Compatible Format
Takes 3 arguments: 
    
     @1 is the input file name (with extendion)
     @2 is the subformat (A1 = A1 survey, TF = Overlap lines, 
        TB = Tellus Border, TE =Tellus)
     @3 is the output file name 
"""
#First import FEM2 Data from Tellus Website, then run the command line below

#aem.process_aem05_data('GSI___15IRL_DLV1563_FEM2.xyz', 'A1', 'A1_FEM2')


"""
#The resulting A1_FEM2_All.dat file contains all the data. Each flight line data are separately stored in A1_FEM2_Sub_L*.dat files. 
#Second step is to specify an area of interest. Takes three arguments:
   
     @1 is the input file name 
     @2 are the left lower and right uper corners in m as [minX maxX minY maxY]
     @3 is the output file name 
"""
#aem.choose_data_rect('A1_FEM2_All.dat','638968.67 641519.17 5922331.93 5924940.46','Sub_StGormans')


""" 
Scatter Data Map Plot of the ALL Flight Lines in the Rectangle
Interpolated Data Map Plot of the ALL Flight Lines in the Rectangle

"""

## Scatter Data Map plot ## 

filename = 'Sub_StGormans_Rectangle'

plots_on = True
mycmap='magma'
#plotformat = '.jpg'
#plotformat = '.eps'
plotformat = '.png'
#plotformat = '.pdf'
#plotformat = '.svg'


Fsize = 16   # Axis label size 
Lsize = 14   # Axis tick label size
s0    = 10   # Marker size
m0 = 's'     # Marker symbol


blocksize = 1     # Blocksize can be set to >1 if you would like to average over certain number of stations using inv.prep_blockdata_aem05
writedata = False



print (' Reading observed data from ' + filename )
data_obs, site_gps, site_alt, site_x, site_y, comment = inv.prep_blockdata_aem05(filename,blocksize,writedata)


I912                  = np.array([data_obs[:,0]])
I3005                 = np.array([data_obs[:,1]])
I11962                = np.array([data_obs[:,2]])
I24510                = np.array([data_obs[:,3]])
Q912                  = np.array([data_obs[:,4]])
Q3005                 = np.array([data_obs[:,5]])
Q11962                = np.array([data_obs[:,6]])
Q24510                = np.array([data_obs[:,7]])

    
IData                 = np.concatenate((I912, I3005, I11962, I24510))
QData                 = np.concatenate((Q912, Q3005, Q11962, Q24510))    
IData_min, IData_max  = np.min(IData), np.max(IData)
QData_min, QData_max  = np.min(QData), np.max(QData)
  

site_x                = np.array([site_x])
site_y                = np.array([site_y])
site_x                = site_x*0.001       #in km
site_y                = site_y*0.001       #in km
#site_x                = site_x.T
#site_y                = site_y.T

site_alt              = site_alt
min_alt               = np.min(site_alt)
max_alt               = np.max(site_alt)
 



r=10 #every 10th value is considered for plotting

if plots_on:
   plt.figure(1)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(12.0,10.0))  #sharex=True and sharey=True are not working anymore - issue date: November 27, 2018 
   cmap = mpl.cm.RdBu
   im1= ax[0, 0].scatter(site_x[::r], site_y[::r], c= I912[::r], s = s0, marker = m0, edgecolors = "None", vmin = IData_min, vmax =IData_max, cmap=cmap)
   ax[0, 0].set_title('In-Phase - 912 Hz',fontsize=Fsize)
   ax[0, 0].set_ylabel('UTM Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   ax[0, 0].axis('equal')
   im2= ax[0, 1].scatter(site_x[::r], site_y[::r], c= I3005[::r], s = s0, marker = m0, edgecolors = "None", vmin = IData_min, vmax =IData_max, cmap=cmap)
   ax[0, 1].set_title('In-Phase - 3005 Hz',fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   ax[0, 1].axis('equal')
   im3= ax[1, 0].scatter(site_x[::r], site_y[::r], c= I11962[::r], s = s0, marker = m0, edgecolors = "None", vmin = IData_min, vmax =IData_max, cmap=cmap)
   ax[1, 0].set_title('In-Phase - 11962 Hz',fontsize=Fsize)
   ax[1, 0].set_ylabel('UTM Northing (km)',fontsize=Fsize)
   ax[1, 0].set_xlabel('UTM Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   ax[1, 0].axis('equal')
   im4=ax[1, 1].scatter(site_x[::r], site_y[::r], c= I24510[::r], s = s0, marker = m0, edgecolors = "None", vmin = IData_min, vmax =IData_max, cmap=cmap)
   ax[1, 1].set_title('In-Phase - 24510 Hz',fontsize=Fsize)
   ax[1, 1].set_xlabel('UTM Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   ax[1, 1].axis('equal')
   
   #fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
  # cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Iscatter'+plotformat)
   plt.show()
   
   plt.figure(2)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(12.0,10.0))
   im1= ax[0, 0].scatter(site_x[::r], site_y[::r], c= Q912[::r], s = s0, marker = m0, edgecolors = "None", vmin = QData_min, vmax =QData_max, cmap=cmap)
   ax[0, 0].set_title('Quadrature - 912 Hz',fontsize=Fsize)
   ax[0, 0].set_ylabel('UTM Northing (km)', fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   ax[0, 0].axis('equal')
   
   im2= ax[0, 1].scatter(site_x[::r], site_y[::r], c = Q3005[::r], s = s0, marker = m0, edgecolors = "None", vmin = QData_min, vmax =QData_max, cmap=cmap)
   ax[0, 1].set_title('Quadrature -  3005 Hz',fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   ax[0, 1].axis('equal')
   im3= ax[1, 0].scatter(site_x[::r], site_y[::r], c = Q11962[::r], s = s0, marker = m0, edgecolors = "None", vmin = QData_min, vmax =QData_max, cmap=cmap)
   ax[1, 0].set_title('Quadrature - 11962 Hz',fontsize=Fsize)
   ax[1, 0].set_ylabel('UTM Northing (km)',fontsize=Fsize)
   ax[1, 0].set_xlabel('UTM Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   ax[1, 0].axis('equal')
   im4= ax[1, 1].scatter(site_x[::r], site_y[::r], c = Q24510[::r],s = s0, marker = m0, edgecolors = "None", vmin = QData_min, vmax = QData_max, cmap=cmap)
   ax[1, 1].set_title('Quadrature - 24510 Hz',fontsize=Fsize)
   ax[1, 1].set_xlabel('UTM Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   ax[1, 1].axis('equal')
   
   fig.subplots_adjust(right=0.91)    
   cb = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cb)
   cb.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Qscatter'+plotformat)
   plt.show()


## Interpolated Data Map Plot ##

RBFfunction='inverse'
numIndexes            = 300

site_x                = np.array([site_x])  #in km
site_y                = np.array([site_y])  #in km
site_x                = site_x.T
site_y                = site_y.T
sitex_min, sitex_max  = np.min(site_x), np.max(site_x)
sitey_min, sitey_max  = np.min(site_y), np.max(site_y)


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



r=1 #since r is set to 1, interpolation takes some time 

if plots_on:
#==============================================================================
   plt.figure(1)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(12.0,10.0))
   cmap = mpl.cm.RdBu
   im1= ax[0, 0].pcolor(XI, YI, interpI912, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('In-Phase - 912 Hz',fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpI3005, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('In-Phase - 3005 Hz',fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpI11962, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('In-Phase - 11962 Hz',fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)', fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpI24510, vmin=IData_min, vmax=IData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('In-Phase - 24510 Hz')
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Iinterpolated'+plotformat) #dpi=1200)
   plt.show()
   #==============================================================================
   
   plt.figure(2)
   fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(12.0,10.0))
   cmap = mpl.cm.RdBu
   im1= ax[0, 0].pcolor(XI, YI, interpQ912, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 0].set_title('Quadrature - 912 Hz',fontsize=Fsize)
   ax[0, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 0].tick_params(labelsize=Lsize)
   ax[0, 0].axis('equal')
   im2= ax[0, 1].pcolor(XI, YI, interpQ3005, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[0, 1].set_title('Quadrature - 3005 Hz',fontsize=Fsize)
   ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[0, 1].tick_params(labelsize=Lsize)
   ax[0, 1].axis('equal')
   im3= ax[1, 0].pcolor(XI, YI, interpQ11962, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 0].set_title('Quadrature - 11962 Hz',fontsize=Fsize)
   ax[1, 0].set_ylabel('Northing (km)',fontsize=Fsize)
   ax[1, 0].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 0].tick_params(labelsize=Lsize)
   ax[1, 0].axis('equal')
   im4= ax[1, 1].pcolor(XI, YI, interpQ24510, vmin=QData_min, vmax=QData_max, cmap=mpl.cm.rainbow)
   ax[1, 1].set_title('Quadrature - 24510 Hz',fontsize=Fsize)
   ax[1, 1].set_xlabel('Easting (km)',fontsize=Fsize)
   ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
   ax[1, 1].tick_params(labelsize=Lsize)
   ax[1, 1].axis('equal')
   fig.subplots_adjust(right=0.91)    
   cbar_ax = fig.add_axes([0.92, 0.14, 0.03, 0.7])
   fig.colorbar(im4, cax=cbar_ax)
   cbar_ax.set_label('(ppm)') # does not work!
   fig = plt.savefig(filename+'-Qinterpolated'+plotformat) #dpi=1200)
   plt.show()

# plt.close('all')