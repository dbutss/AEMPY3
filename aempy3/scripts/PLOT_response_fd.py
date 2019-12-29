# /usr/bin/python

"""
Plots MH Results of TDEM Data

"""

import numpy as np
import matplotlib.pyplot as plt


Fsize = 18
Lsize = 15

# Load observed data
#filename       = 'FWDModel_genesis_3LRES5_E01_thkfirst-15m.npz' 
filename       = 'FWDModel_aem05_LR-0.3-0.3-0.3_E-60.0.npz'
tmp            = np.load(filename)
data_obs       = tmp['data_obs']  
data_true       = tmp['data_true'] 
data_err        = tmp['data_err'] 

I_obs          = np.concatenate(([data_obs[:,0]], [data_obs[:,1]], [data_obs[:,2]], [data_obs[:,3]]))
Q_obs          = np.concatenate(([data_obs[:,4]], [data_obs[:,5]], [data_obs[:,6]], [data_obs[:,7]]))
I_true          = np.concatenate(([data_true[0]], [data_true[1]], [data_true[2]], [data_true[3]]))
Q_true          = np.concatenate(([data_true[4]], [data_true[5]], [data_true[6]], [data_true[7]]))
freq= np.array([912.0, 3005.0, 11962.0, 24510])
# 
fig = plt.figure()
ax1 = fig.add_subplot(111) 
ax1.plot(freq, I_obs, 'o', color= 'red', markersize=7, label='In-Phase')
ax1.plot(freq, Q_obs, 'o',  color= 'black', markersize=7, label='Quadrature')
ax1.plot(freq, I_true, '+', color= 'red', markersize=9, label='In-Phase')
ax1.plot(freq, Q_true, '+',  color= 'black', markersize=9, label='Quadrature')

#upperlimits = X_obs+data_err
#lowerlimits = X_obs-data_err
ax1.errorbar(freq, I_obs, yerr=data_err,color= 'red')
ax1.errorbar(freq, Q_obs, yerr=data_err,color= 'black')
#
ax1.set_xlabel('freq (Hz)',fontsize=Fsize)
ax1.set_ylabel('ppm',fontsize=Fsize)
ax1.xaxis.set_label_position('top')
ax1.xaxis.set_ticks_position('both')
plt.xscale('log')
##plt.yscale('log')
ax1.tick_params(labelsize=Lsize)
#
ax1.grid(True)
#plt.legend(loc='best', numpoints=1, fontsize=12)
