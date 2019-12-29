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
filename       = 'FWDModel_genesis_3L_RES_derr-0.01-500.npz'
tmp            = np.load(filename)
data_obs       = tmp['data_obs']  
data_err       = tmp['data_err'] 
X_obs          = np.concatenate(([data_obs[:,0]],[data_obs[:,1]],[data_obs[:,2]],[data_obs[:,3]],[data_obs[:,4]],
                                 [data_obs[:,5]],[data_obs[:,6]],[data_obs[:,7]],[data_obs[:,8]],[data_obs[:,9]], [data_obs[:,10]]))
Z_obs          = np.concatenate(([data_obs[:,11]],[data_obs[:,12]],[data_obs[:,13]],[data_obs[:,14]],[data_obs[:,15]],
                                 [data_obs[:,16]],[data_obs[:,17]],[data_obs[:,18]],[data_obs[:,19]],[data_obs[:,20]], [data_obs[:,21]]))
 
T_obs= np.array([0.009, 0.026, 0.052, 0.095, 0.156, 0.243, 0.365, 0.547, 0.833, 1.259, 1.858])

 
fig = plt.figure()
ax1 = fig.add_subplot(111) 
#

upperlimits = X_obs+data_err
lowerlimits = X_obs-data_err
ax1.errorbar(T_obs, X_obs, yerr=data_err[0:11],color= 'red')
ax1.errorbar(T_obs, Z_obs, yerr=data_err[11:22],color= 'black')

#
ax1.plot(T_obs, X_obs, 'o', color= 'red', markersize=7, label='In-line Component')
ax1.plot(T_obs, Z_obs, 'o', color= 'black', markersize=7, label='Vertical Component')


ax1.set_xlabel('time (ms)',fontsize=Fsize)
ax1.set_ylabel('data (fT)',fontsize=Fsize)
ax1.xaxis.set_label_position('top')
ax1.xaxis.set_ticks_position('both')
plt.xscale('log')
#plt.yscale('log')
ax1.tick_params(labelsize=Lsize)

ax1.grid(True)
plt.legend(loc='best', numpoints=1, fontsize=12)
