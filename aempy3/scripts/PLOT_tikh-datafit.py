#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 2018

@author: duygu

"""



import sys
import os
import math  as ma
import numpy as np
import scipy as sc

#import aemplots as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col


mycmap=mpl.cm.gist_rainbow


plotformat = '.png'
#plotformat = '.eps'
Fsize = 20
Lsize = 16



filename='A1_3098D0060m_ERTsection_Rectangle_Npca3_TikhGCV_BlkSiz1_nlyr35_res100.0_tau0_-2.0_tau1_0.698970004336_DErr_30.0_Results.npz'
print (' Data read from '+filename)

FileName,filext0 = os.path.splitext(filename)

#title=FileName
title = 'Bundoran Test Line 3098 D0 60m'

tmp = np.load(filename)
site_model = tmp['site_model']
site_error = tmp['site_error']
site_sens  = tmp['site_sens']
site_num   = tmp['site_num']
site_num   = abs(site_num)
site_data  = tmp['site_data']

site_x     = tmp['site_x']
site_y     = tmp['site_y']
                 
             
site_rms         = site_data[:,0]
Icalc_9kHz       = site_data[:,1]
Icalc_3kHz       = site_data[:,2]
Icalc_12kHz      = site_data[:,3]
Icalc_25kHz      = site_data[:,4]
Qcalc_9kHz       = site_data[:,5]
Qcalc_3kHz       = site_data[:,6]
Qcalc_12kHz      = site_data[:,7]
Qcalc_25kHz      = site_data[:,8]


Iobs_9kHz        = site_data[:,9]
Iobs_3kHz        = site_data[:,10]
Iobs_12kHz       = site_data[:,11]
Iobs_25kHz       = site_data[:,12]
Qobs_9kHz        = site_data[:,13]
Qobs_3kHz        = site_data[:,14]
Qobs_12kHz       = site_data[:,15]
Qobs_25kHz       = site_data[:,16]


Data_min = -500.0
Data_max = 2500.0  

#plt.figure(1)
#fig, ax =  plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(12,10))  
#
#ax[0, 0].set_title('In-Phase - 912 Hz',fontsize=Fsize)
##ax1.set_ylim([Data_min,Data_max])
#ax[0, 0].set_ylabel('(ppm)',fontsize=Fsize)
#ax[0, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
#ax[0, 0].tick_params(labelsize=Lsize)
#
#   
#ax[0, 1].plot(site_num,Icalc_3kHz,'-r',site_num,Iobs_3kHz,'.k')
#ax[0, 1].set_title('In-Phase - 3005 Hz',fontsize=Fsize)
##ax[0, 1].set_ylim([Data_min,Data_max])
#ax[0, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
#ax[0, 1].tick_params(labelsize=Lsize)
#
#ax[1, 0].plot(site_num,Icalc_12kHz,'-r',site_num,Iobs_12kHz,'.k')
#ax[1, 0].set_title('In-Phase - 11962 Hz',fontsize=Fsize)
#ax[1, 0].set_xlabel('site number',fontsize=Fsize)
#ax[1, 0].set_ylabel('(ppm)',fontsize=Fsize)
##ax[0, 1].set_ylim([Data_min,Data_max])
#ax[1, 0].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
#ax[1, 0].tick_params(labelsize=Lsize)
#
#
#ax[1, 1].plot(site_num,Icalc_25kHz,'-r',site_num, Iobs_25kHz,'.k')
#ax[1, 1].set_title('In-Phase - 24510 Hz',fontsize=Fsize)
#ax[1, 1].set_xlabel('site number',fontsize=Fsize)
##ax[0, 1].set_ylim([Data_min,Data_max])
#ax[1, 1].grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
#ax[1, 1].tick_params(labelsize=Lsize)
#
#plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
#
#plt.legend(['Calculated', 'Observed'], fontsize=Fsize, loc=4)
#fig = plt.savefig(filename+'-1'+plotformat)
#plt.show()
#


plt.figure(1)
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

ax1.plot(site_num,Icalc_9kHz,'-r',site_num,Iobs_9kHz,'.r',linewidth=2)
ax1.plot(site_num,Icalc_3kHz,'-g',site_num,Iobs_3kHz,'.g',linewidth=2)
ax1.plot(site_num,Icalc_12kHz,'-k',site_num,Iobs_12kHz,'.k', linewidth=2)
ax1.plot(site_num,Icalc_25kHz,'-b',site_num,Iobs_25kHz,'.b', linewidth=2)

ax1.set_xlabel('site number',fontsize=Lsize)
ax1.set_ylabel('(ppm)',fontsize=Lsize)
ax1.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
ax1.tick_params(labelsize=Lsize)
ax1.set_ylim([Data_min,Data_max])
ax1.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
ax1.tick_params(labelsize=Lsize)



ax2.plot(site_num,Qcalc_9kHz,'-r',site_num,Qobs_9kHz,'.r',linewidth=2)
ax2.plot(site_num,Qcalc_3kHz,'-g',site_num,Qobs_3kHz,'.g',linewidth=2)
ax2.plot(site_num,Qcalc_12kHz,'-k',site_num, Qobs_12kHz,'.k',linewidth=2)
ax2.plot(site_num,Qcalc_25kHz,'-b',site_num,Qobs_25kHz,'.b',linewidth=2)

ax2.set_xlabel('site number',fontsize=Lsize)

ax2.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
ax2.tick_params(labelsize=Lsize)
ax2.set_ylim([Data_min,Data_max])
ax2.set_xlabel('site number',fontsize=Lsize)
ax2.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
ax2.tick_params(labelsize=Lsize)


#plt.legend(['912 Hz', '3005 Hz', '11962 Hz', '24510 Hz'], fontsize=Fsize, loc=4)
fig = plt.savefig(filename+'-4'+plotformat)
plt.show()



f = plt.figure(figsize=(15,2))
ax = f.add_subplot(111)
ax.plot(site_num, site_rms,'.k')
ax.set_ylabel('nRMS', fontsize=Lsize)
ax.set_xlabel('site number', fontsize=Lsize)
ax.tick_params(labelsize=Lsize)
ax.grid(color='k', alpha=0.5, linestyle='dotted', linewidth=1.5)
fig = plt.savefig(filename+'-3'+plotformat)






