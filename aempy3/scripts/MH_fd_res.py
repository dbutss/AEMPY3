#! /usr/bin/python
# for  running under python 2.X andf 3.X
from __future__ import print_function, absolute_import, division

import sys, os
import numpy as np
import scipy as sc
import random
from datetime import datetime

import invprocs as inv
import core1d
#


print ('\n  AEM Metropolis-Hastings \n \n') 


random.seed(datetime.now())
#for test fix the random seed
#random.seed(2000)`

aem_system = 'aem05'
call_fwd = 'core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,m_initial)'


nan = np.nan 
alt=60.

#sizedat=8
mode = 0
data_err = 30.
#data_err = 0.


nsample = 1000000
nburnin = 10000
accp_check = 1000
accp_max = 60. 
accp_min = 30. 
step0 =.1
stepfac=0.6666
mode=0
onlyactive = False
Datafile='FWDModel_aem05_LR-100.0-10.0-100.0_E-30.0.npz'
print (' Reading observed data from ' + Datafile)
true_res = ([100.0, 5.0, 100.0]) # case A

tmp= np.load(Datafile)
data_obs=tmp['data_obs']
sizedat = np.shape(data_obs);
data_std=tmp['data_err']*np.ones(sizedat)
data_cov= data_std*data_std;
ndata=sizedat[1];

data_old = data_obs
data_obs = data_old #    + data_err*np.random.randn(1, sizedat[0])


sizedat = np.shape(data_obs);
nsite=sizedat[0]
ndata=sizedat[1]
dactive = np.ones(ndata,dtype='int8')
#dactive[0]=0
#dactive[4]=0


nlyr = 12
znlyr = np.zeros(nlyr);iznlyr=znlyr.astype(int)
onlyr = np.ones(nlyr) ;ionlyr=onlyr.astype(int)


mactive, prior_avg, prior_std, m_upper, m_lower = inv.initialize_1dmod(nlyr)


mactive[0*nlyr:1*nlyr]   = ionlyr[0:nlyr] 
sizepar  = np.shape(mactive); 
sizeact  = np.shape(np.flatnonzero(mactive))
mpara=sizeact[0];

dzstart = .7
dzend   = .7
dz = np.logspace(dzstart,dzend,nlyr-1)+0.25*np.random.randn(nlyr-1)
z = np.append( 0., np.cumsum(dz)) 
x=np.zeros(np.shape(z))
y=np.zeros(np.shape(z))
#

#dz = ([25.0, 25.0])
#z = np.append( 0., np.cumsum(dz)) 


guess_r = 100. #initial guess for resistivity in prior_avg
guess_mu = 1.0 # relative magnetic permeability 
guess_ep = 1.0 #  relative dielectric constant 
guess_m = 0.0 # chargeability
guess_t = 0.0 # time constant
guess_c = 1.0 # freq c constant


prior_avg[0*nlyr:1*nlyr] = np.log10(guess_r*onlyr[0:nlyr])
prior_avg[1*nlyr:2*nlyr] = guess_mu*onlyr[0:nlyr]
prior_avg[2*nlyr:3*nlyr] = guess_ep*onlyr[0:nlyr]
prior_avg[3*nlyr:4*nlyr] = guess_m*onlyr[0:nlyr]
prior_avg[4*nlyr:5*nlyr] = guess_t*onlyr[0:nlyr]
prior_avg[5*nlyr:6*nlyr] = guess_c*onlyr[0:nlyr]
prior_avg[6*nlyr:7*nlyr-1] = dz[0:nlyr-1]

#prior_std defines standard deviation of prior_avg
guess_sr=0.2
res_std=guess_sr*onlyr[0:nlyr];


prior_std[0*nlyr:1*nlyr] = res_std;
#prior_atd[1*nlyr:2*nlyr] = guess_mu*onlyr[0:nlyr]
#prior_std[2*nlyr:3*nlyr] = guess_ep*onlyr[0:nlyr]
#prior_std[3*nlyr:4*nlyr] = guess_m*onlyr[0:nlyr]
#prior_std[4*nlyr:5*nlyr] = guess_t*onlyr[0:nlyr]
#prior_std[5*nlyr:6*nlyr] = guess_c*onlyr[0:nlyr]
#prior_std[6*nlyr:7*nlyr-1] = 0.1*dz[0:nlyr-1]

Lx=40.
Ly=40.
Lz=40.
L = np.array([Lx,Ly,Lz])

Cm_prior=np.identity(np.size(prior_avg))
prior_var= res_std*res_std
Cm_prior_res = inv.covp3d(x,y,z,L,var=prior_var)
Cm_prior[0*nlyr:1*nlyr,0*nlyr:1*nlyr] = Cm_prior_res


# sc.linalg.block_diag()
LC = sc.linalg.cholesky(Cm_prior)



m_upper[0*nlyr:1*nlyr]  =  6. #prior_avg[mactive==1] + 3*prior_std[mactive==1]
m_lower[0*nlyr:1*nlyr]  = -1. #prior_avg[mactive==1] - 3*prior_std[mactive==1]

print \
(' Parameter set for inverting: \n', mactive)
print \
(' Layer thicknesses: \n', dz)
print \
(' Layer interface depths: \n', z)
print \
(' Initial homogeneous halfspace resistivity of %6.2f Ohm*m' % (guess_r))
print \
(' Log10 Standard error of %6.2f ' % (guess_sr))

print \
(res_std)

print \
(' Upper limits: \n', m_upper)
print \
(' Lower limits: \n', m_lower)

print ('\n  starting simulation \n') 

#Metropolis-Hastings algorithm with nsample iteration

m_gauss   = np.random.randn(sizepar[0])
m_initial = prior_avg # + step0*prior_std*mactive*m_gauss #*r

fp =np.sum(np.power((m_initial-prior_avg)/prior_std,2)); 
fp_old=fp;

data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,m_initial)')
dc_old = data_calc
ss_old = np.sum(np.power(abs(data_obs - data_calc)/data_std,2))

m_old   = m_initial;


scount = 0 
step = step0

if onlyactive==False:
    mpara = 7*nlyr
    mchain = np.zeros((nsample+1,mpara+ndata+1))
    mchain[scount,0:mpara]= m_old
else:
    mchain = np.zeros((nsample+1,mpara+ndata+1))
    mchain[scount,0:mpara]= m_old[mactive==1]

mchain[scount,mpara:mpara+ndata]= dc_old
mchain[scount,mpara+ndata]= ss_old


reject = 0
accept = 0
for n in range(nsample):
    
#    print('\n',n)
    #draw random numbers from the prior model (prior_avg)
    m_g = np.random.randn(sizepar[0])
#    print(' m_g   ',m_g[0*nlyr:1*nlyr])
    m_c = np.matmul(LC,m_g)
    m_s = step*mactive*m_c
#    print(' m_s   ', m_s[0*nlyr:1*nlyr])
#    print(' m_g   ', m_g[0*nlyr:1*nlyr])
#    print(' m_c   ', m_c[0*nlyr:1*nlyr])
    #m_sample = m_old + tep*prior_std*mactive*np.random.randn(sizepar[0])
    m_sample = m_old + m_s
#    print(' m_smp   ', m_sample[1*nlyr:2*nlyr])
    if np.all(m_sample < m_upper) and np.all(m_sample > m_lower):
        scount = scount+1
        m_new = m_sample; 
        fp_new=np.sum(np.power((m_new-m_old)/prior_std,2)); 
        #data_calc = core1d.aemfwd1d_gtk4(mode,alt,nlyr,m_new)
        data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,alt,nlyr,m_new)')
        dc_new = data_calc
        ss_new = np.sum(np.power(abs(data_obs - data_calc)/data_std,2))
    #fl_new=np.exp(-0.5*ss_new)
        alpha = min(1,np.exp(-0.5*(ss_new-ss_old) -0.5*(fp_new-fp_old)));
    #print (alpha) 
        p = np.random.rand(1,1)      
        if p < alpha:
        # accept new candidate

            accept = accept+1

            if onlyactive:
                  mchain[scount,0:mpara]= m_new[mactive==1]
            else:
                  mchain[scount,0:mpara]= m_new

            mchain[scount,mpara:mpara+ndata]= dc_new
            mchain[scount,mpara+ndata]= ss_new

            m_old  = m_new
            ss_old = ss_new;
            fp_old = fp_new;
            dc_old = dc_new
        else:      
            if onlyactive:
                  mchain[scount,0:mpara]= m_old[mactive==1]
            else:
                  mchain[scount,0:mpara]= m_old

            mchain[scount,mpara:mpara+ndata]= dc_old
            mchain[scount,mpara+ndata]= ss_old
    else:
            reject=reject+1
    
    if np.mod(n,accp_check)==0 and not n==0:        
        accpp=100*accept/scount
        if accpp > accp_max: step = step/stepfac
        if accpp < accp_min: step = step*stepfac
        print (' percentage of samples accepted = %6.2f of %6i, with  %6i out of bounds, step set to %6.4f'\
        % (accpp, scount, reject, step))

           
print (' final percentage of samples accepted = %6.2f ' % (100*accept/scount))
#np.savez('Mh_chain',mchain)
print(np.shape(mchain[nburnin:scount,1]))
# np.savez('Mh_chain1',mchain=mchain)
#np.savez_compressed('Mh_chain2',mchain=mchain[nburnin:scount,:])


avg_res = np.zeros(nlyr)
std_res = np.zeros(nlyr)
med_res = np.zeros(nlyr)
prc_res =  np.zeros((nlyr,2))
for ilyr in range(nlyr):
    avg_res[ilyr] = np.mean(mchain[nburnin:scount,ilyr])
    med_res[ilyr] = np.median(mchain[nburnin:scount,ilyr])
    std_res[ilyr] = np.std(mchain[nburnin:scount,ilyr])
    prc_res[ilyr,0] = np.percentile(mchain[nburnin:scount,ilyr],2.3) # 15.9)
    prc_res[ilyr,1] = np.percentile(mchain[nburnin:scount,ilyr],97.7)# 84.1)
    print (' Layer %3i with center at %4.1f m Resistivity - Mean: %6.2f Median: %6.2f' % (ilyr, z_center[ilyr],avg_res[ilyr], med_res[ilyr]) )

filen,filext = os.path.splitext(Datafile) 
filout=filen+'_MHchain'+'_'+str(nsample)+'_Lz'+str(Lz)+filext
np.savez_compressed(filout,mchain=mchain[nburnin:scount,:],nlyr=nlyr,\
                     avg_res = avg_res, med_res = med_res, std_res=std_res,\
                     prc_res = prc_res)
   

