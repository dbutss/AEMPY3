# /usr/bin/python

"""
Plots results from 'MH_fd_resdz.py' 
Single-site results: Running Metropolis-Hastings algorithm for one site
 
"""

import numpy as np
import pylab
import matplotlib.pyplot as plt

nan= float('NaN')

# Load "observed" data
filename       = 'FWDModel_gtk4_3LRES5_thkfirst-15m.npz' 
tmp            = np.load(filename)
data_obs       = tmp['data_obs']  
I_obs          = np.concatenate(([data_obs[:,0]], [data_obs[:,1]], [data_obs[:,2]], [data_obs[:,3]]))
Q_obs          = np.concatenate(([data_obs[:,4]], [data_obs[:,5]], [data_obs[:,6]], [data_obs[:,7]]))


filename ='FWDModel_gtk4_3LRES5_thkfirst-15m_dz3l.npz'
tmp = np.load(filename)
nlyr= tmp['nlyr']
avg_res=tmp['avg_res']
med_res=tmp['med_res']
std_res=tmp['std_res']
prc_res=tmp['prc_res']
med_dz =tmp['med_dz'] 
mchain=tmp['mchain']

res_val = mchain[:,0:nlyr]
res_val = np.power(10,res_val)
sizemchain = np.shape(mchain)
ncol       = sizemchain[1] 
data_calc  = mchain[:,ncol-9:-1]


I912   = data_calc[:,0]
I3005  = data_calc[:,1]
I11962 = data_calc[:,2]
I24510 = data_calc[:,3]
Q912   = data_calc[:,4]
Q3005  = data_calc[:,5]
Q11962 = data_calc[:,6]
Q24510 = data_calc[:,7]

#Plot every 50th results
step= 50;
sample_accptd= np.arange(1,len(I912),step);
freq= np.array([912.0, 3005.0, 11962.0, 24510])


fig = plt.figure()
ax1 = fig.add_subplot(111) 

for n in sample_accptd:
    data_calc_I = np.concatenate(([I912[n]], [I3005[n]], [I11962[n]], [I24510[n]]))
    data_calc_Q = np.concatenate(([Q912[n]], [Q3005[n]], [Q11962[n]], [Q24510[n]]))
   
    ax1.plot(freq, data_calc_I[:], linewidth= .7, color= [.75, .75, .75]);
    ax1.plot(freq, data_calc_Q[:], linewidth= .7, color= [.55, .55, .55])    
    plt.xscale('log')
    plt.xlim([100,100000])
    ax1.set_xlabel('frequency (Hz)',fontsize=18)
    ax1.set_ylabel('data (ppm)',fontsize=18)
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(labelsize=15)

ax1.plot(freq, I_obs, 'o', color= 'red', markersize=7, label='In-Phase')
ax1.plot(freq, Q_obs, 'o',  color= 'black', markersize=7, label='Quadrature')


ax1.grid(True)
plt.legend(loc='best', numpoints=1, fontsize=12)
