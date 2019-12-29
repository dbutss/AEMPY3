#! /usr/bin/python
# for  running under python 2.X andf 3.X
from __future__ import print_function, absolute_import, division

import os
import sys

import time
import math  as ma
import numpy as np

import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import random
from datetime import datetime


import aemprocs as aem
import invprocs as inv

import core1d

print(' AEM MAP Estimation ' )

nan = np.nan # float('NaN')
random.seed(datetime.now())


#for test fix the random seed
#random.seed(2000)`



mode=0
invmeth= 'par'
maxiter = 50
threshrms1 =0.1
threshrms2 =1.e-4
linsmin=5
linsfac=0.1

meth=invmeth.lower()[0:3]

if   meth == 'dat':
    print (' Running data space iteration ')
elif meth == 'par':
    print (' Running parameter space iteration ')
elif meth == 'pqn':
    print (' Running qn-parameter space iteration ')
else:
    sys.exit(' inversion method '+invmeth+' not implemented!')


filename = 'BTML_L1455_Radar-less100m'
blocksize = 20
writedata = False
aem_system = 'gtk4'
# filename = 'StGormans_L1382'
# blocksize = 10
print(' Reading observed data from ' + filename) 
data_obs, site_gps, site_alt, site_x, site_y, comment =eval('inv.prep_blockdata_'+aem_system+'(filename,blocksize,writedata)')

# print 'data_obs',data_obs[0]
sizedat = np.shape(data_obs);
nsite=sizedat[0]
ndata=sizedat[1]
dactive = np.ones(sizedat,dtype='int8')
dactive[data_obs==np.nan]=0

filename = 'StGormans_L1382'
blocksize = 10
writedata = False
aem_system = 'gtk4'
ErrDat = 30.

#data_std= 30*np.ones(sizedat[1])
#data_std= data_std[dactive==1]
#Nerr= la.norm(data_std)
#Wd = np.diagflat(1./data_std,0)
##Wd = sp.csr_matrix(Wd)
#
#ndata=sizedat[1];
#Cd_prior = np.diagflat(data_std*data_std,0)
#if invmeth.lower()[1:3] == 'par' or invmeth.lower()[1:2] == 'qn':
#    Cd_inv   =  la.inv(Cd_prior)
#    Wd = la.sqrtm(Cd_inv)
#

# setup model 
nlyr=16
znlyr = np.zeros(nlyr);iznlyr=znlyr.astype(int)
onlyr = np.ones(nlyr) ;ionlyr=onlyr.astype(int)


mactive, prior_avg, prior_std, m_upper, m_lower = inv.initialize_1dmod(nlyr)

mactive[0*nlyr:1*nlyr]   = ionlyr[0:nlyr]
sizepar  = np.shape(mactive); #sizepar = sizepar[:]
sizeact  = np.shape(np.flatnonzero(mactive))
mpara=sizeact[0];


dzstart = 0.5
dzend   = 1.
dz = np.logspace(dzstart,dzend,nlyr)
# dz = np.linspace(dzstart,dzend,nlyr)
z = np.append( 0., np.cumsum(dz)) 
prior_avg[6*nlyr:7*nlyr-1] = dz[0:nlyr-1]

#initial guess for resistivity in prior_avg
guess_r = 100.
res=guess_r*onlyr[0:nlyr];
res = np.log10(guess_r);
prior_avg[0*nlyr:1*nlyr]   = res

#prior_std defines standard deviation of prior_avg
guess_s=0.2
res_std=guess_s*onlyr[0:nlyr];
prior_std[0*nlyr:1*nlyr] = res_std
prior_cov= prior_std*prior_std


print \
(' Parameter set for inverting: \n', mactive)
#print \
#(' Layer thicknesses: \n', dz)
#print \
#(' Layer interface depths: \n', z)
#print \
#(' Initial homogeneous halfspace resistivity of %6.2f Ohm*m' % (guess_r))
#print \
#(' Log10 Standard error of %6.2f ' % (guess_s))
#print \
#(' Upper limits: \n', m_upper)
#print \
#(' Lower limits: \n', m_lower)



Lz=66.;Ly=30;
print \
(' Assuming exponential z covariance with correlation length %6.2f m' \
% (Lz))
#print \
#(' Assuming exponential y (inline) covariance with correlation length %6.2f m' \
#% (Ly))
std=np.mat(res_std*res_std).T
Cm_prior  = inv.covp_expnt(dz,Lz,std)
#print(np.shape(Cm_prior))
# Cm_prior = la.block_diag(Cr_prior,np.identity(6*nlyr-1))
Cm_prior = inv.collapse_covar(Cm_prior,mactive)
#print(Cm_prior)
        
if meth == 'par' or meth == 'pqn':
    Cm_inv   = la.inv(Cm_prior)
    Wm = la.sqrtm(Cm_inv)
    
# print(Cm_inv)
# Bayesian inversion with deterministic method 
mode= 0
# ser initial values
model_gauss   = np.random.randn(sizepar[0])
model_ini     = prior_avg + prior_std*mactive*model_gauss #*r


m_prior       = inv.collapse_model(prior_avg,mactive);
m_prior=np.mat(m_prior).T
site_results = np.zeros((nsite,2+7*nlyr-1+2*ndata))
site_sens = np.zeros((nsite,mpara))


start = time.time()

for ii in range(sizedat[0]):

    data_std= ErrDat*np.ones(sizedat[1])
    data_std= data_std[dactive[ii,:]==1]
    Nerr= la.norm(data_std)
    Wd = np.diagflat(1./data_std,0)
    Wd = sp.csr_matrix(Wd)
    Cd_prior = np.diagflat(data_std*data_std,0)
    
    if invmeth.lower()[1:3] == 'par' or invmeth.lower()[1:2] == 'qn':
       Cd_inv   =  la.inv(Cd_prior)
       Wd = la.sqrtm(Cd_inv)

    sizact = np.shape(data_std);
    nd=sizact[0]


    print ('\n'+filename+':   site #%6i' % (ii) )           
    

    
    dobs = np.mat(data_obs[ii,dactive[ii,:]==1])
    dstd = np.mat(data_std)
    
    tau = 1.
    bet = 1.
    
    Cm_post = []
    Gi      = []
    Jac     = []
    J       = []
    
    niter = 0
    nrms_iter = 1.e30
    nrms_old  = nrms_iter
    model_new     = model_ini
    model_old     = model_new 

    while (niter < maxiter): 
    
        niter =niter +1               
   
        m_old = model_new[mactive==1]
        model_old = model_new
 

        data_calc, Jac = eval('core1d.aemjac1d_'+aem_system+'(mode,site_alt[ii],nlyr,model_new,mpara,mactive)')
        
        dcalc = np.mat(data_calc[ii,dactive[ii,:]==1])
        resid = (dobs - dcalc)/dstd

        nrms_iter, srms_iter = inv.rms(dcalc,dobs,Wd) 
        print (' iteration %6i NRMS =  %7.3f  SRMS = %4.1f percent' % (niter, nrms_iter,srms_iter) )
        
        if nrms_iter > 
        
        
        Jac=Jac[dactive[ii,:]==1,:]
        Jac=sp.csr_matrix(Jac)
#        sens=np.sum(Wd*Jac,axis=0)
        

        m_iter = inv.collapse_model(model_new,mactive)
        m_iter=np.mat(m_iter).T
        
    
        
        if   invmeth.lower()[0:3] == 'dat':
             d =(dobs-dcalc).T + Jac*np.mat(m_iter-m_prior)
             J = sp.csr_matrix(Cm_prior*Jac.T)
             A = Jac*J + tau*Cd_prior
             r = d
             m_update = m_prior+ J*spla.spsolve(A,r)
             m_delta  = m_update - m_old
             
        elif invmeth.lower()[0:3] == 'par':
             d =(dobs-dcalc).T + (Jac*np.mat(m_iter-m_prior))
             J = sp.csr_matrix(Jac.T*Cd_inv)
             A = J*Jac+ tau*Cm_inv
             r = (Jac.T*Cd_inv)*d
             print (np.shape(A))   
             print (np.shape(r))     
             m_update = m_prior+spla.spsolve(A,r)
             m_delta  = m_update - m_old

        elif invmeth.lower()[0:3] == 'pqn':
             J = sp.csr_matrix(Jac.T*Cd_inv)
             A = (Jac.T*Cd_inv*Jac+Cm_prior)
             r = (Jac.T*Cd_inv)*(dobs-dcalc)-Cm_inv*np.mat(m_iter-m_prior)
             m_update = m_old+spla.spsolve(A,r)
             m_delta  = m_update - m_old
    #         print mm
                     
    
    #             #sys.exit('parameter space inversion not yet implemented!')
        model_new.flat[mactive==1] = m_update
        
#                dcalc = np.mat(data_calc[ii,dactive[ii,:]==1])
#        resid = (dobs - dcalc)/dstd
#
#        nrms_iter, srms_iter = inv.rms(dcalc,dobs,Wd) 
#        print (' iteration %6i NRMS =  %7.3f  SRMS = %4.1f percent' % (niter, nrms_iter,srms_iter) )
#        
        
        print (nrms_iter,nrms_old)
        
        
        if nrms_iter > nrms_old:
           print (' iteration %6i, starting line search...' %(niter))
           bet = 1.
           nrms=0.
           while bet > linsmin and nrms > nrms_iter:
               bet=bet*linsfac
               print (' iteration %6i, beta = %5.3f' %(niter,bet))
               model_new.flat[mactive==1]=bet*m_delta + m_old
               data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,site_alt[ii],nlyr,model_new)')
               nrms, dummy = inv.rms(dcalc[dactive==1],dobs,Wd)
    #    print  model_new
        
    elapsed = (time.time() - start)
    print (' Used %7.4f sec for %6i sites' % (elapsed, ii))
    np.savez_compressed('MAP_'+meth+'_'+comment+'_site'+str(ii),mactive=mactive, model_new=model_new,
                        data_new=data_calc,data_obs=data_obs,data_std=data_std)

    line = np.concatenate((np.array([ii]),np.array([nrms_iter]),\
                           np.array(model_new),np.array(data_calc),\
                           np.array(data_obs[ii])))
    site_results[ii] =line 
    site_sens[ii] = sens
             
comment = ''             
np.savez_compressed('MAP_results'+filename+'_'+meth,\
                    mactive=mactive, site_results=site_results,\
                    site_sens=site_sens,site_conv=conv,\
                    nlyr=nlyr, site_y=site_y,site_x=site_x,site_gps=site_gps,site_alt=site_alt)

    
