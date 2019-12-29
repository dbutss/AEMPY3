#! /usr/bin/python

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

print(' AEM TIKH Estimation ')

#folder = os.path.relpath('../core/')
#if folder not in sys.path:
#    sys.path.append(folder)
#    
#folder = os.path.relpath('../modules/')
#if folder not in sys.path:
#    sys.path.append(folder)    
    
    
nan = np.nan # float('NaN')
random.seed(datetime.now())
#for test fix the random seed
#random.seed(2000)`

out = True
mode=0


# system 
aem_system = 'aem05'
pars_system= ''

# data
filename = 'A1_FEM2_Sub_L1379_Npca3'
blocksize = 10
writedata = False
DataErr=30.
# window_width=3
# window_shift=2

read_prior = False
filename='A1_FEM2_Sub_L1379_Npca3_TikhFIX_BlkSiz5_nlyr36_900Hz_Results.npz'



if read_prior:
   print (' Prior model read from '+filename)
   FileName,filext0 = os.path.splitext(filename)
   tmp = np.load(filename)
   prior_model = tmp['site_model']

# model definition
# mesh
nlyr=36
dzstart = 4.
dzend   = 12.
#initial guess for resistivity in prior_avg
guess_r = 300.
#prior_std defines standard deviation of prior_avg
guess_s=0.5



mode= 0
maxiter = 30
# dont overfit 
threshrms1 =0.5
# small changee
threshrms2 =1.e-2
# converged?
threshrms3 =3.e-0
# 'FIX', 'GCV', 'UPRE'
Gfunction = 'FIX'
# Gfunction = 'GCV'

n0=1
tau0min=np.log10(0.001) #0.001
tau0max=np.log10(0.001)   #10
n1=1
tau1min=np.log10(60.)
tau1max=np.log10(60.)
tau0 = np.logspace(tau0min,tau0max,n0)
tau1 = np.logspace(tau1min,tau1max,n1)


identstr ='_Tikh'+Gfunction+'_BlkSiz'+str(blocksize)+'_nlyr'+str(nlyr)+'_900Hz'


print(' Reading observed data from ' + filename) 
data_obs, site_gps, site_alt, site_x, site_y, identstr_in = inv.prep_blockdata_aem05(filename,blocksize,writedata)
# print 'data_obs',data_obs[0]
sizedat = np.shape(data_obs);
nsite=sizedat[0]
ndata=sizedat[1]
d_active = np.ones(sizedat,dtype='int8')
#d_active[:,0] = 0    # 912 Hz
#d_active[:,4] = 0    # 912 Hz
#

znlyr = np.zeros(nlyr);iznlyr=znlyr.astype(int)
onlyr = np.ones(nlyr) ;ionlyr=onlyr.astype(int)

m_active, prior_avg, prior_std, m_upper, m_lower = inv.initialize_1dmod(nlyr)
m_active[0*nlyr:1*nlyr]   = ionlyr[0:nlyr]
sizepar  = np.shape(m_active); 
sizeact  = np.shape(np.flatnonzero(m_active))
mpara=sizeact[0];



dz = np.logspace(np.log10(dzstart),np.log10(dzend),nlyr)
# dz = np.linspace(dzstart,dzend,nlyr)
z = np.append( 0., np.cumsum(dz)) 
prior_avg[6*nlyr:7*nlyr-1] = dz[0:nlyr-1]


res=guess_r*onlyr[0:nlyr];
res = np.log10(guess_r);
prior_avg[0*nlyr:1*nlyr]   = res


res_std=guess_s*onlyr[0:nlyr];
prior_std[0*nlyr:1*nlyr] = res_std
prior_cov= prior_std*prior_std 

max_val= 6
min_val=-1
m_upper[m_active==1] =  max_val #prior_avg[m_active==1] + 3*prior_std[m_active==1]
m_lower[m_active==1] =  min_val #prior_avg[m_active==1] - 3*prior_std[m_active==1]


if out:
#   print \
#   (' Parameter set for inverting: \n', m_active)
   print \
   (' Layer thicknesses: \n', dz)
   print \
   (' Layer interface depths: \n', z)
   print \
   (' Initial homogeneous halfspace resistivity of %6.2f Ohm*m' % (guess_r))
   print \
   (' Log10 Standard error of %6.2f ' % (guess_s))
#   print \
#   (' Upper limits: \n', m_upper)ndata
#   print \
#   (' Lower limits: \n', m_lower)
   print \
   (' Assuming conventional difference operators L0 and L1' )
   print \
   (' RegparFunction is ', Gfunction)
   print \
   (' tau 0 = \n', tau0)
   print \
   (' tau 1 = \n', tau1)




# tikhonov quasi-inversion with deterministic method 

L0 = inv.diffops(dz,der=False,otype='L0') 
L1 = inv.diffops(dz,der=False,otype='L1') 


#model_new_gauss   = np.random.randn(sizepar[0])
model_ini     = prior_avg # + model_std*m_active*model_gauss #*r
error_ini     = prior_std
m_apr         = np.mat(inv.collapse_model(prior_avg,m_active))


site_model = np.zeros((nsite,sizepar[0]))
site_error = np.zeros((nsite,sizepar[0]))
site_sens  = np.zeros((nsite,sizepar[0]))
site_data  = np.zeros((nsite,1+2*sizedat[1]))
site_num   = np.zeros((nsite))

m_test = np.zeros((n0*n1,mpara))
m_err  = np.zeros((n0*n1,mpara))
gcv    = np.zeros((n0*n1,1))
upr    = np.zeros((n0*n1,1))
dnorm  = np.zeros((n0*n1,1))
mnorm  = np.zeros((n0*n1,1))
tau    = np.zeros((n0*n1,2))




start = time.time()

for ii in range(sizedat[0]):
    if out:
       print ('site: %i of %i, blksiz=%i' % (ii,sizedat[0],blocksize) )
    
    d_active_site = d_active[ii,:]
    d_active_site = np.ravel(d_active_site,1)
    d_obs_site = data_obs[ii,:]
    d_obs_site = d_obs_site[d_active_site==1]
    d_obs_site = np.ravel(d_obs_site,1)
    d_std_site=  DataErr*np.ones(np.shape(d_obs_site))
    dvar= d_std_site*d_std_site
    Nerr= la.norm(d_std_site)
    Wd = np.diagflat(1./d_std_site,0)
    Wd = sp.csr_matrix(Wd)
    Cdi=Wd.T* Wd

    sizact = np.shape(d_std_site);
    nd=sizact[0]

    
    niter = 0
    nrms_iter = 1.e30
    nrms_old  = nrms_iter
    model     = model_ini
    error     = error_ini
    calc      = np.nan*np.ones(sizedat[1])
    sens      = np.nan*np.ones(sizepar[0])
    model_old = model     #        print \
    error_old = error
    calc_old  = calc
       
    
    while (niter < maxiter) and (nrms_iter > threshrms1): 
    
        niter =niter +1
        
        data_calc = \
        eval('core1d.aemfwd1d_'+aem_system+'(mode,site_alt[ii],nlyr,model'+pars_system+')')
        d_calc_site = np.mat(data_calc[d_active_site==1])
        calc = data_calc
        
        resid     = Wd*(d_obs_site-d_calc_site).T
#        print ' observed data: ', d_obs
#        print ' calculated data: ', d_calc
#        print ' r_iter: ', resid
        #nrms_iter     = np.sqrt(np.sum(np.power(abs(resid),2))/(nd-1))
        nrms_iter, srms_iter = inv.rms(d_calc_site,d_obs_site,Wd)
        if out:
           print (' iteration %6i NRMS =  %7.3f  SRMS = %4.1f percent' % (niter, nrms_iter,srms_iter) )
#       print (' iteration %6i NRMS =  %7.3f' % (niter, nrms_iter) )
        
        if nrms_iter < nrms_old*(1.-threshrms2):
            site_num[ii] = ii
            nrms_old  = nrms_iter
            model_old = model
            error_old = error
            calc_old  = calc

            
        else:
            site_num[ii] = ii
            if nrms_iter > threshrms3: site_num[ii] = -ii        

            if out:
               print (' iteration %6i NRMS(new) =  %7.3f > NRMS(old) =  %7.3f'\
               % (niter, nrms_iter,nrms_old))
            model = model_old
            error     = error_old
            nrms_iter = nrms_old
            calc      = calc_old
            break
            
    
        
        
#        data_calc, Jac = core1d.aemjac1d_aem05(mode,site_alt[ii],nlyr,model,mpara,m_active)
        data_calc, Jac\
        = eval('core1d.aemjac1d_'+aem_system+'(mode,site_alt[ii],nlyr,model,mpara,m_active'+pars_system+')')
        m_iter = np.mat(inv.collapse_model(model,m_active))
#        print ' m_iter: ', m_iter
        diff_m     = np.mat(m_iter-m_apr)
        Jac=Jac[d_active[ii,:]==1,:]
        Jac=sp.csr_matrix(Jac)
        #sens=np.sum(Wd*Jac,axis=0)
        
        sensmod=np.sum(np.abs(Wd*Jac),axis=0)
        
        m_test = np.zeros((n0*n1,mpara))
        m_err  = np.zeros((n0*n1,mpara))
        
        gcv    = np.zeros((n0*n1,1))
        upr    = np.zeros((n0*n1,1))
        ufc    = np.zeros((n0*n1,1))
        dnorm  = np.zeros((n0*n1,1))
        mnorm  = np.zeros((n0*n1,1))
        tau    = np.zeros((n0*n1,2))
           
        
        
        it = -1
        for t0 in tau0:
            S0=ma.sqrt(t0)*L0
            Cmi0 = S0.T*S0
            for t1 in tau1:
                S1=ma.sqrt(t1)*L1
                Cmi1 = S1.T*S1
                
                model_test= model
                it = it + 1 
                tau[it,0]=t0
                tau[it,1]=t1
                A = (Jac.T*Cdi*Jac+Cmi0+Cmi1)
                r = (Jac.T*Cdi)*(d_calc_site-d_obs_site).T+Cmi0*diff_m.T+Cmi1*diff_m.T
                m_delta = spla.spsolve(A,r)
                
                
                m_test[it] = m_iter - m_delta
                model_test.flat[m_active==1] = m_test[it]
#                if ii==28:
#                   print(model_test.flat[m_active==1])
                # d_calc = core1d.aemfwd1d_aem05(mode,site_alt[ii],nlyr,model_test)
                data_calc = eval('core1d.aemfwd1d_'+aem_system+'(mode,site_alt[ii],nlyr,model_test)')
                d_calc_site = np.mat(data_calc[d_active_site==1])
                resit = Wd*(d_obs_site-d_calc_site).T
                
                dnorm[it] = la.norm(Wd*(d_obs_site-d_calc_site).T)
                mnorm[it] = la.norm(S0.T*S0*diff_m.T+S1.T*S1*diff_m.T)
                                #print (np.shape(A))
#               # Cov a posteriori
                C = la.inv(A)  
                m_err[it]=np.sqrt(np.diag(C))
#                C =sp.linalg.inv(A)
#               # Generalized Inverse
                G=  C*Jac.T*Wd.T
                M = Wd*Jac*G
                
                gcv[it] = nd*ma.pow(dnorm[it],2)/ma.pow(np.trace(np.eye(nd)-M),2)
                upr[it] = ma.pow(dnorm[it],2)/nd - 2.*ma.pow(Nerr,2)*(1.-np.trace(M)/nd)
                #ufc[it] = 1./ dnorm[it] + 1./ mnorm[it] 
 
            if n0*n1 == 1: 
                g_index = 0
            
            elif Gfunction.lower()[0:3]=='upr':
                g_index = np.argmin(gcv)
                g0=tau[g_index,0]
                g1=tau[g_index,1]
                if out: 
                   print \
                   ( ' Min GCV at tau0 = %7.3g, tau1 - %7.3g ' % (g0,g1))

            elif Gfunction.lower()[0:3]=='upr':
                g_index = np.argmin(upr)
                g0=tau[g_index,0]
                g1=tau[g_index,1]
                if out:
                   print \
                   ( ' Min UPR at tau0 = %7.3g, tau1 - %7.3g ' % (g0,g1))

            else:
                sys.exit(Gfunction.lower()[0:3]+' not iplemented!')
#        g_model = m_test[g_index,:]
#        g_err   = m_err[g_index,:] 
        model.flat[m_active==1] = m_test[g_index,:]
        error.flat[m_active==1] = nrms_iter*nrms_iter*m_err[g_index,:]
        sens.flat[m_active==1]  = sensmod

        
    elapsed = (time.time() - start)
    print (' Used %7.4f sec for %6i sites \n' % (elapsed, ii+1))
    
    site_data[ii]  = np.concatenate((np.array([nrms_iter]),calc.flat, data_obs[ii,:].flat))
    site_model[ii] = model.flat 
    site_sens[ii]  = sens.flat
    site_error[ii] = error.flat
    #site_regpar   =  
np.savez_compressed(filename+identstr+'_Results',\
                    m_active=m_active,\
                    d_active =d_active,\
                    nlyr=nlyr,\
                    blocksize=blocksize,\
                    site_num=site_num,\
                    site_model=site_model,\
                    site_sens=site_sens,\
                    site_error=site_error,\
                    site_data=site_data,\
                    site_y=site_y,site_x=site_x,site_gps=site_gps,site_alt=site_alt)
