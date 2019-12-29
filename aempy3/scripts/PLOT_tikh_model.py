#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Show several 1d block models as (stitched) section.

@authors: duygu, vrath

"""
# for  running under python 2.X andf 3.X
from __future__ import print_function, absolute_import, division 

import sys
import os
import math  as ma
import numpy as np
import scipy as sc


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from invprocs import centers

#mycmap=mpl.cm.gist_rainbow
mycmap = mpl.cm.jet                  # Define colourmap


#plotformat = '.eps'
plotformat = '.png'
FSize = 20
LSize = 17
min_lrho=1.5
max_lrho=3.5
sl = [min_lrho, max_lrho]

topo_use_average = False # True

low_sens = True
if low_sens:
   lowsens =1e-3 #1e-4 ????
   
   
max_depth = False
if max_depth:
   maxdepth=80
   
high_rms = True
if high_rms:
   highrms = 10.
   
plotdepth = -10.
plot_depth_adapt= False

plotheight = 60.
plot_height_adapt= True 
blank = 5


filename='A1_3098D0060m_ERTsection_Rectangle_Npca3_TikhGCV_BlkSiz1_nlyr35_res100.0_tau0_-2.0_tau1_0.698970004336_DErr_30.0_Results.npz'
print (' Data read from '+filename)

FileName,filext0 = os.path.splitext(filename)


title = 'Bundoran Test Line: A1 3098 D0 60m'    




tmp          = np.load(filename)
site_model   = tmp['site_model']
site_error   = tmp['site_error']
site_sens    = tmp['site_sens']
site_num     = tmp['site_num']
site_num     = abs(site_num)
site_data    = tmp['site_data']
site_rms     = site_data[:,0]
site_conv    = site_num/np.abs(site_num)
site_conv[0] = 1.

m_active     = tmp['m_active']
nlyr         = tmp['nlyr']

site_x       = tmp['site_x']
site_y       = tmp['site_y']
site_alt     = tmp['site_alt']
site_gps     = tmp['site_gps']


site_ref     = site_alt + 57.0            # Geoid Height is approximately 57.0 for Ireland
site_topo    = site_gps - site_ref        # Assuming that Digital Terrain Elevation is not distributed by GSI


models=np.shape(site_model)
sites=models[0]
param=models[1]
site_val=site_model[:,m_active==1]
site_thk=site_model[:,6*nlyr:7*nlyr-1]
site_r=np.sqrt(np.power(site_x,2.) +np.power(site_y,2.))
site_r=site_r-site_r[0]
site_sens=site_sens[:,m_active==1]
scale=np.max(np.abs(site_sens.flat))
site_sens=site_sens/scale
nlay = nlyr
dxmed2 = np.median(np.diff(site_r)) / 2.


#**************************Plots******************************
#plt.style.use('ggplot')


fig, ax = plt.subplots(figsize=(15, 8))   

max_topo = max(site_topo)

if low_sens:
##   site_val[np.abs(site_sens) < lowsens]=np.nan
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
## print(plotmin1)

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
   plotmin = np.max([plotmin1,plotmin2])
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
ax.set_ylabel('elevation (m)', fontsize=FSize)
ax.yaxis.set_label_position("right")
ax.set_xlabel(' profile distance (m)', fontsize=FSize)
ax.tick_params(labelsize=FSize)
#ax.grid(b='on')
#
#    pg.mplviewer.createColorbar(p, cMin=cmin, cMax=cmax, nLevs=5)
#
cb = plt.colorbar(p, orientation='horizontal',aspect=35,pad=0.15)
xt = [-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5]
cb.set_ticks( xt, [str(xti) for xti in xt] )
cb.set_label('$\log_{10}$ ($\Omega$ m)', size=FSize)
cb.ax.tick_params(labelsize=FSize) 

#    plt.draw()
#    return fig, ax
print (' PLot to '+FileName+plotformat)
fig = plt.savefig(FileName+plotformat,dpi=300)