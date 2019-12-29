#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Show several 1d block models as (stitched) section."""
# for  running under python 2.X andf 3.X
from __future__ import print_function, absolute_import, division 

import sys
import os
import math  as ma
import numpy as np
import scipy as sc

#import aemplots as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from invprocs import centers

mycmap=mpl.cm.gist_rainbow


plotformat = '.png'
#plotformat = '.eps'
FSize = 16
LSize = 12
min_lrho=1.
max_lrho=3.5 
sl = [min_lrho, max_lrho]

topo_use_average = False # True

low_sens = True
if low_sens:
   lowsens =1e-3
   
max_depth = True
if max_depth:
   maxdepth=50
   
high_rms = True
if high_rms:
   highrms = 4
   
high_err = False
if high_err:
   higherr = 0.5 
   
plotdepth = 60.
plot_depth_adapt= False# True

plotheight = 140.
plot_height_adapt= False# True
blank = 10



#filename='TIKH_results_A1_FEM2_Sub_L1379_Npca3_nlyr-40_tau0-0.01_tau1-100_block5.npz'
filename='A1_FEM2_Sub_L1379_Npca3_TikhFIX_BlkSiz5_nlyr36_900Hz_Results.npz'
print (' Data read from '+filename)

FileName,filext0 = os.path.splitext(filename)

#title=FileName
title = 'A1_FEM2_Sub_L1379_Pca3'

tmp = np.load(filename)
site_model = tmp['site_model']
site_error = tmp['site_error']
site_sens  = tmp['site_sens']
site_num   = tmp['site_num']
site_data  = tmp['site_data']
site_rms=site_data[:,0]
site_conv=site_rms/np.abs(site_rms)

m_active=tmp['m_active']
nlyr=tmp['nlyr']

site_x=tmp['site_x']
site_y=tmp['site_y']
site_alt=tmp['site_alt']
site_gps=0.1*tmp['site_gps']
site_topo=site_gps-site_alt -56#-(site_gps[0]-site_alt[0])
if topo_use_average:
   site_tref=np.mean(site_topo)
else:
   site_tref=0
site_topo=site_topo-site_tref



models=np.shape(site_model)
sites=models[0]
param=models[1]
site_val=site_model[:,m_active==1]
site_thk=site_model[:,6*nlyr:7*nlyr-1]
site_r=np.sqrt(np.power(site_x,2.) +np.power(site_y,2.))
site_r=site_r-site_r[0]
#site_r=np.absolute(site_r-site_r[0])
#site_r=site_r-site_r[sitelen[0]-1]
scale=np.max(np.abs(site_sens.flat))
site_sens=site_sens/scale
site_sens=site_sens[:,m_active==1]

nlay = nlyr
dxmed2 = np.median(np.diff(site_r)) / 2.


# Plots
plt.style.use('ggplot')

#mpl.rcParams['text.usetex'] = True


fig, ax = plt.subplots(figsize=(20, 10))



max_topo = max(site_topo)

if low_sens:
#    site_val[np.abs(site_sens) < lowsens]=np.nan
   nmod=0
   plotmin = max_topo
   while nmod < sites:
       invalid=[]
       thk = site_thk[nmod,:]
       thk_halfspace =  thk[-1]
       thk = np.hstack((thk.T, thk_halfspace)) 
       z0 = np.hstack((0., np.cumsum(thk)))
       zm= 0.5*(z0[0:nlay]+z0[1:nlay+1])
       invalid = np.abs(site_sens[nmod,:])<lowsens
       site_val[nmod,invalid]=np.nan       
       zm = site_topo[nmod] -zm
       if any(invalid) :
          plotmin=np.min([plotmin,np.min(zm[invalid])])
       nmod=nmod+1
   plotmin1 = plotmin
# print(plotmin1)



if max_depth:
   nmod=0
   plotmin =max_topo 
   while nmod < sites:
       invalid=[]
       thk = site_thk[nmod,:]       
       thk_halfspace =  thk[-1]
       thk = np.hstack((thk.T, thk_halfspace)) 
       z0 = np.hstack((0., np.cumsum(thk)))
       zm = 0.5*(z0[0:nlay]+z0[1:nlay+1]) 
       zm = np.hstack((zm,zm[nlay-1]))  
       invalid = zm>maxdepth
       site_val[nmod,invalid[0:nlay]]=np.nan
       zm = site_topo[nmod] -zm
       if any(invalid) :
          plotmin=np.min([plotmin,np.min(zm[invalid])])
       nmod=nmod+1
   plotmin2 = plotmin
# print(plotmin2)
if high_err:
#    site_val[np.abs(site_sens) < lowsens]=np.nan
   nmod=0
   plotmin = max_topo
   while nmod < sites:
       invalid=[]
       thk = site_thk[nmod,:]
       thk_halfspace =  thk[-1]
       thk = np.hstack((thk.T, thk_halfspace)) 
       z0 = np.hstack((0., np.cumsum(thk)))
       zm= 0.5*(z0[0:nlay]+z0[1:nlay+1])
       invalid = np.abs(site_sens[nmod,:])<lowsens
       site_val[nmod,invalid]=np.nan       
       zm = site_topo[nmod] -zm
       if any(invalid) :
          plotmin=np.min([plotmin,np.min(zm[invalid])])
       nmod=nmod+1
       
   plotmin3 = plotmin
    
site_val = np.ma.masked_invalid(site_val)
val=np.zeros(np.shape(site_val))

patches = []
nmod = 0
while nmod < sites:
    if high_rms and site_rms[nmod]<highrms:
       val[nmod, :] = site_val[nmod,:]
       thk = site_thk[nmod,:]
       thk_halfspace =  thk[-1]
       thk = np.hstack((thk.T, thk_halfspace)) 
       z0 = np.hstack((0., np.cumsum(thk))) 
       z = site_topo[nmod]-z0
       for j in range(nlay):
           rect = Rectangle((site_r[nmod] - dxmed2, z[j]), dxmed2 * 2, thk[j])
           patches.append(rect)
    nmod=nmod+1
       
p = PatchCollection(patches, cmap=plt.get_cmap(mycmap) , linewidths=0,edgecolor=None)
p.set_clim(sl)
p.set_array(site_val.ravel())
ax.add_collection(p)

if plot_depth_adapt:
   plotmin = np.max([plotmin1,plotmin2,plotmin3])
else:
   plotmin = plotdepth
   
if plot_height_adapt:
   plotmax = max_topo
else:
   plotmax = plotheight
   
   
ax.set_ylim((plotmax+blank,plotmin-blank))
ax.set_xlim((min(site_r) - dxmed2, max(site_r) + dxmed2))
ax.invert_yaxis()
#ax.axis('equal')
#    if title is not None:
ax.set_title(title, fontsize=FSize)
ax.set_ylabel('depth (m)', fontsize=FSize)
ax.yaxis.set_label_position("right")
ax.set_xlabel(' profile distance (m)', fontsize=FSize)
ax.tick_params(labelsize=FSize)
ax.grid(b='on')
#
#    pg.mplviewer.createColorbar(p, cMin=cmin, cMax=cmax, nLevs=5)
#
cb = plt.colorbar(p, orientation='horizontal',aspect=40,pad=0.1)
xt = [-1, -1.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
cb.set_ticks( xt, [str(xti) for xti in xt] )
cb.set_label(r'log10($\Omega$ m)', size=FSize)
#cb.set_label('log10(Omega m)', size=FSize)
cb.ax.tick_params(labelsize=FSize) 
#
#    plt.draw()
#]#    return fig, ax
print (' PLot to '+FileName+plotformat)
fig = plt.savefig(FileName+plotformat)