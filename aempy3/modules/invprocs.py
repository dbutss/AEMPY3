# -*- coding: utf-8 -*-
"""

Tools for Airborne EM

authors: VR/DK/RD

"""
#==============================================================================
import sys
import os
from sys import exit as error
from datetime import datetime
import random
import math  as ma
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.signal import medfilt
import core1d
    
def get_filelist(mystr='*',mypath='.'):
    """
    geberates filelist from path and wildcard
    author: vrath
    """
    import sys
    import os
    import fnmatch
    
    filelist = fnmatch.filter(os.listdir(mypath), mystr)
    return filelist
    
    
def set_relpath(relpath="../modules"):

     sys.path.insert(1,os.path.join(os.path.dirname(__file__), relpath))  

#def l_corner(r_norm,m_norm,regpar):
#    '''
#    L-curve method for determining the 'optimim' regularization parameter 
#    '''
#    alpha = 0.
#    npoints=np.shape(r_norm)
#    r_norm=np.log10(r_norm)
#    m_norm=np.log10(m_norm)
#    
#    return alpha
 
def generate_random_cov(m0, mactive, covar, seed=None):
   '''
   generates random models with prescribed covariance
   '''
   sizepar=np.shape(m0)
   
   if seed==None:
      random.seed(datetime.now())
   else:
      random.seed(seed)
      
   m_gauss   = np.random.randn(sizepar[0])
    
   L = la.cholesky(covar)
   
   m = m0+mactive*m_gauss*L
   
   return m


#def aniso_euclidean_norm(x1, x2, scale=None):
#    '''
#    anisotropic scaling for radial basic function interpolation
#    '''   
#    if scale==None:
#        fn=np.sqrt( ((x1 - x2)**2).sum(axis=0) )
#    else:
#        scale=scale/np.sum(scale)
#        f = scale*(x1 - x2)**2
#        fn=np.sqrt( (f).sum(axis=0) )
#    return fn
    
def rnormp(data_obs,data_calc,data_std=None,p=2):
    '''
     Calculates the p-norm of the residuals 
    '''  
    if data_std is None:
        data_std=np.ones(np.shape(data_obs))
        
    resid=(data_obs - data_calc)/data_std

    rnormp=np.power(resid,p)
    rnorm=np.sum(rnormp)
#    return {'rnorm':rnorm, 'resid':resid }
    return (rnorm,resid)
    
def rms(dcalc,dobs,Wd=1.):
    '''
     Calculates the NRMS ans SRMS
    '''  
    sizedat=np.shape(dcalc)
    nd=sizedat[1]
    rscal     = Wd*(dobs-dcalc).T
    #print(sizedat,nd)
    # normalized root mean square error
    nrms      = np.sqrt(np.sum(np.power(abs(rscal),2))/(nd-1))  
    #sum squared scaled symmetric error
    serr      = 2.*nd*np.abs(rscal)/(abs(dobs.T)+abs(dcalc.T)) 
    ssq       = np.sum(np.power(serr,2))
    srms      = 100. * np.sqrt(ssq/nd)
    
    return nrms, srms

def prep_tiles(DataIn=None, TileOut=None, TileSize=250.,Header=None):       
    """
        blocks data into tiles for fast interpolation 
        takes 2 arguments: 
    
            @1 is the input data (with extension)
            @2 is the output tile file 
               
    """   
    X = DataIn[1,:]
    Y = DataIn[2,:]

    X_min, X_max  = np.min(X), np.max(X)
    Y_min, Y_max  = np.min(Y), np.max(Y)
    
    X_N = int(ma.floor((X_max-X_min)/TileSize)) 
    Y_N = int(ma.floor((Y_max-Y_min)/TileSize)) 
    
    X_T = X_min +range(X_N)*TileSize
    Y_T = X_min +range(Y_N)*TileSize
    
    
    for ii in range(range(X_N)
         for jj in range(range(Y_N) 
            
            indices = np.where(x > 10)
             
def initialize_obsdata(D,dvalue = None):
    '''
    Initializes observational data structure inversion as all active, 
    all incative, or NaN.
    '''
    onobs = np.ones(np.shape(D))
    if dvalue.lower() =='active':
       dactive =onobs.astype(int)
    elif dvalue.lower() =='inactive':
       dactive = 0*onobs.astype(int)
    elif dvalue.lower() =='nan':
       dactive =np.nan*onobs.astype(int)
    else:
       error('d_active initial  '+dvalue+' not implemeted !')

    return dactive
 
    
def initialize_1dmod(nlyr):
    '''
    Initializes data structures for 1D inversion 
    parameters as all inactive (mactive)
    priors (here: gaussian) to reasonable values for res/dz-only inversion
    bounds to practical +-infinite
    
    '''
    # import numpy as np
    znlyr = np.zeros(nlyr)
    iznlyr=znlyr.astype(int)
    onlyr = np.ones(nlyr) 
    ionlyr=onlyr.astype(int)

#mactive determines the active paramteters, currently 7 per layer
    mactive = np.zeros((7*nlyr))
    mactive[0*nlyr:1*nlyr]   = iznlyr[0:nlyr]
    mactive[1*nlyr:2*nlyr]   = iznlyr[0:nlyr]
    mactive[2*nlyr:3*nlyr]   = iznlyr[0:nlyr]
    mactive[3*nlyr:4*nlyr]   = iznlyr[0:nlyr]
    mactive[4*nlyr:5*nlyr]   = iznlyr[0:nlyr]
    mactive[5*nlyr:6*nlyr]   = iznlyr[0:nlyr]
    mactive[6*nlyr:7*nlyr-1] = iznlyr[0:nlyr-1]
    
    # sizepar  = np.shape(mactive); 
    mprior_p1 = np.zeros(np.shape(mactive))
    mprior_p1[0*nlyr:1*nlyr]   = onlyr[0:nlyr]
    mprior_p1[1*nlyr:2*nlyr]   = onlyr[0:nlyr]
    mprior_p1[2*nlyr:3*nlyr]   = onlyr[0:nlyr]
    mprior_p1[3*nlyr:4*nlyr]   = znlyr[0:nlyr]
    mprior_p1[4*nlyr:5*nlyr]   = znlyr[0:nlyr]
    mprior_p1[5*nlyr:6*nlyr]   = znlyr[0:nlyr]
    mprior_p1[6*nlyr:7*nlyr-1] = onlyr[0:nlyr-1]
    
    mprior_p2= np.ones(np.shape(mactive))
    
    u_max=1e30
    u_min=-1e30
    mupper=u_max*np.ones(np.shape(mactive))
    mlower=u_min*np.ones(np.shape(mactive))
  
#    m_upper[mactive==1] =  6. #prior_avg[mactive==1] + 3*prior_std[mactive==1]
#    m_lower[mactive==1] = -1. #prior_avg[mactive==1] - 3*prior_std[mactive==1]

    
    return mactive, mprior_p1, mprior_p2, mupper, mlower
    
    
def collapse_covar(M, isactive):
    
    tmp = M[isactive==1,:]
    A = tmp[:,isactive==1]
    return A
    
    
def collapse_model(M, isactive):
    A = M[isactive==1]
    return A
    
    
#def do_linesearch(delta_m, m0, alt, dobs, mactive, dactive, nrms_0, rfac=0.666, maxreduce=6):
#    
#    niter = 0
#    fact  = 1.  
#    nrms  = 0.
#
#    mode =  0
#    while (niter < maxreduce) and (nrms >= nrms_0):
#        m_new = m0 + fact *delta_m# and (nrms_iter > nrms_0): 
#        d_calc = core1d.aemfwd1d_aem05(mode,alt,nlyr,m_new)
#        dcalc = np.mat(d_calc[dactive==1])
#        resid     = Wd*(dobs-dcalc).T
#        nrms     = np.sqrt(np.sum(np.power(abs(resid),2))/(nd-1))
#        niter-niter+1
#        fact=fact*rfac
#    return m_iter 
#    
def shift(l, n):
    return l[n:] + l[:n]
 
def diffops(dz,der=True,otype='L0',mtype='dense',mform='csr'):
    '''
    Differential operators L0-L2 based on dz
    '''
    layers=np.shape(dz)   
    nlyr=layers[0]
    
    if otype == 'L0':
        d = np.ones((1, nlyr))
        L = sp.spdiags(d, [0], nlyr, nlyr,format=mform)
    elif otype == 'L1' :
        z = np.append( 0., np.cumsum(dz))
        zc = 0.5*(z[0:nlyr-1]+z[1:nlyr])
        zc = np.append(zc,zc[nlyr-2]+dz[nlyr-2])
        h  = 1./np.diff(zc)
        if der != True:
            h=np.ones(np.shape(h))
        d = np.zeros((2, nlyr-1))
        d[0,:] =-h        
        d[1,:] = np.roll(h,1)
        #L = sp.spdiags(d, [0,1], nlyr, nlyr,format='csr')
        #L[nlyr-1,nlyr-1] = 0.
        L = sp.spdiags(d, [0,-1], nlyr, nlyr-1,format=mform)
        L=L.transpose()
        
#        print(np.shape(d))
#        print(np.shape(L))
        
    else:
        error('DiffOperator '+otype+' not implemeted !')

    if mtype=='dense':
        L =L.todense()

    return L#  
    
def diffops3d(dz,der=True,otype='L0',mtype='dense',mform='csr'):
    '''
    Differential operators L0-L2 based on dz
    '''
    layers=np.shape(dz)   
    nlyr=layers[0]
    
    if otype == 'L0':
        d = np.ones((1, nlyr))
        L = sp.spdiags(d, [0], nlyr, nlyr,format=mform)
    elif otype == 'L1' :
        z = np.append( 0., np.cumsum(dz))
        zc = 0.5*(z[0:nlyr-1]+z[1:nlyr])
        zc = np.append(zc,zc[nlyr-2]+dz[nlyr-2])
        h  = 1./np.diff(zc)
        if der != True:
            h=np.ones(np.shape(h))
        d = np.zeros((2, nlyr-1))
        d[0,:] =-h
        d[0,0] = 1.
        d[0,1] = 0.
        d[1,:] = h
        L = sp.spdiags(d, [0,-1], nlyr, nlyr-1,format=mform)
        
        
        L=L.transpose()
    else:
        error('DiffOperator '+otype+' not implemeted !')

    if mtype=='dense':
        L =L.todense()

    return L#  

    
    
def suppops(dz, m, seps=1.e-4,otype='MS',mtype='sparse',mform='csr',der=True):
    '''
    calculate minimum support operators MS or MGS 
    '''
    
    if otype == 'MS':
        mabs=np.abs(m)
        M = mabs/(mabs+seps)

    elif otype == 'MGS' :
        layers=np.shape(dz)   
        nlyr=layers[0]+1
        z = np.append( 0., np.cumsum(dz))
        zc = 0.5*(z[0:nlyr-1]+z[1:nlyr])
        zc = np.append(zc,zc[nlyr-2]+dz[nlyr-2])
        h  = 1./np.diff(zc)        
        d = np.zeros((2, nlyr))
        d[0,:] =-h
        d[1,:] = h
        L=sp.spdiags(d, [0,-1], nlyr, nlyr-1,format=mform)
        L=L.transpose()
        
        dabs =np.abs(L*m)
        M= dabs/(dabs+seps)
    
    W=la.inv(M)
    
    if mtype=='dense':
        M =M.todense()
        W =W.todense()
    
    return M, W 
    
    
def barrops(m,mlower,mupper,otype='log',mtype='sparse',mform='csr'):
    '''
    Barrier operator used for enforcing upper/lower limits
    for parameters
    
    '''
    sizem = np.shape(m);
    nm=sizem[0]
    if otype=='log':
        d=np.log(m-mlower)+np.log(mupper-m)
    elif otype=='inv':
        d=1./(m-mlower)+1./(mupper-m)
    else:
        error('BarrOperator '+otype+' not inplemeted !')

    B=sp.spdiags(d, [0], nm, nm,format=mform)

    return B
    
def error_model(data_obs,data_add,data_mult):
    '''
    Error model including multiplicative and additive noise
    following Brodie (2015) GALEISBSTDEM
    
    '''
    nd = np.size(data_obs)
    na = np.size(data_add)
    nm = np.size(data_mult)
    
    if na[0]==1: 
       data_a=data_add*np.ones(nd[0], nd[1])
    else:
       data_a=data_add
    if nm[0]==1: 
       data_m=data_mult*np.ones(nd[0], nd[1])
    else:
       data_m=data_mult

    data_err=np.sqrt(np.power(data_m*data_obs,2)+np.power(data_a,2))

    return data_err

def centers(dz):
    '''
    defines cell centers 
    '''
    layers=np.shape(dz)
    nlyr= layers[0]
    nlyr1=nlyr+1
    
#    print(dz[nlyr-1])
    dz =np.append(dz,dz[nlyr-1])
    z = np.append( 0., np.cumsum(dz))
    zc = 0.5*(z[0:nlyr1-1]+z[1:nlyr1])
    
    return zc

def covd3d(x,y,z,rot=([0., 0., 0.]),ctype='exp',L=([1., 1., 1.]),var=[]):
    '''
    Exponential 3D covariance with correlation correlation  lengths Lx,Ly,Lz
    '''
#    npoints=np.shape(z)
#    
    shapez=np.shape(z)
#    shapex=np.shape(x)
#    shapey i=np.shape(y)
    npoints=shapez[0] 
         
    Cd=np.zeros((npoints,npoints))
    
    for il in range(npoints-1):
        for jl in range(il,npoints-1):
            Cd[jl,il]= np.exp(\
            -np.abs(x[il]-x[jl])/L[1] \
            -np.abs(y[il]-x[jl])/L[2] \
            -np.abs(z[il]-z[jl])/L[3] \
            )
            Cd[il,jl]= Cd[jl,il]
            
          
    if len(var)==0:
       return Cd
    else:
       Cd = Cd*np.diag(var)
    
    return Cd

def covp3d(x,y,z,rot=([0., 0., 0.]),ctype='exp',L=([1., 1., 1.]),var=[]):
    '''
    Exponential 3D covariance with correlation correlation  lengths Lx,Ly,Lz
    '''
    # npoints=np.shape(z)
    shapez=np.shape(z)
#    shapex=np.shape(x)
#    shapey=np.shape(y)
    npoints=shapez[0] 
#    print(shapex,shapey,shapez)
#    print(np.shape(L))
    Cp_expnt=np.zeros((npoints,npoints))
#    print(np.shape(Cp_expnt))
    for il in range(npoints):
        for jl in range(il,npoints):
#            print(np.abs(z[il]))
#            print(np.abs(z[jl]))
#            print(np.abs(L[3]))
            Cp_expnt[jl,il]= np.exp(\
            -np.abs(x[il]-x[jl])/L[0] \
            -np.abs(y[il]-x[jl])/L[1] \
            -np.abs(z[il]-z[jl])/L[2] \
            )
            Cp_expnt[il,jl]= Cp_expnt[jl,il]
      
    if len(var)==0:
       return Cp_expnt
    else:
       Cp_expnt = Cp_expnt*np.diag(var)
                             
    return Cp_expnt
     

def covd_expnt3d(x,y,z,L=([1., 1., 1.]),var=[]):
    '''
    Exponential 3D covariance with correlation correlation  lengths Lx,Ly,Lz
    '''
#    npoints=np.shape(z)
#    
    shapez=np.shape(z)
#    shapex=np.shape(x)
#    shapey i=np.shape(y)
    npoints=shapez[0] 
         
    Cd_expnt=np.zeros((npoints,npoints))
    
    for il in range(npoints-1):
        for jl in range(il,npoints-1):
            Cd_expnt[jl,il]= np.exp(\
            -np.abs(x[il]-x[jl])/L[1] \
            -np.abs(y[il]-x[jl])/L[2] \
            -np.abs(z[il]-z[jl])/L[3] \
            )
            Cd_expnt[il,jl]= Cd_expnt[jl,il]
            
          
    if len(var)==0:
       return Cd_expnt
    else:
       Cd_expnt = Cd_expnt*np.diag(var)
    
    return Cd_expnt
    
def covp_expnt3d(x,y,z,L=([1., 1., 1.]),var=[]):
    '''
    Exponential 3D covariance with correlation correlation  lengths Lx,Ly,Lz
    '''
    # npoints=np.shape(z)
    shapez=np.shape(z)
#    shapex=np.shape(x)
#    shapey=np.shape(y)
    npoints=shapez[0] 
#    print(shapex,shapey,shapez)
#    print(np.shape(L))
    Cp_expnt=np.zeros((npoints,npoints))
#    print(np.shape(Cp_expnt))
    for il in range(npoints):
        for jl in range(il,npoints):
#            print(np.abs(z[il]))
#            print(np.abs(z[jl]))
#            print(np.abs(L[3]))
            Cp_expnt[jl,il]= np.exp(\
            -np.abs(x[il]-x[jl])/L[0] \
            -np.abs(y[il]-x[jl])/L[1] \
            -np.abs(z[il]-z[jl])/L[2] \
            )
            Cp_expnt[il,jl]= Cp_expnt[jl,il]
      
    if len(var)==0:
       return Cp_expnt
    else:
       Cp_expnt = Cp_expnt*np.diag(var)
                             
    return Cp_expnt
 
   
def covd_gauss3d(x,y,z,L=([1., 1., 1.]),var=[]):
    '''
    Gaussian 3D covariance with correlation  lengths Lx,Ly,Lz
    '''
    shapez=np.shape(z)
#    shapex=np.shape(x)
#    shapey=np.shape(y)
    npoints=shapez[0] 
    
    Cd_gauss=np.zeros((npoints,npoints))
    for il in range(npoints):
        for jl in range(il,npoints):
            Cd_gauss[jl,il]= np.exp(\
            -0.5*np.power(np.abs(x[il-1]-x[jl-1])/L[0],2)\
            -0.5*np.power(np.abs(y[il-1]-y[jl-1])/L[1],2)\
            -0.5*np.power(np.abs(z[il-1]-z[jl-1])/L[2],2)\
            )
            Cd_gauss[il,jl]= Cd_gauss[jl,il]
        #enforce positive definiteness
    thresh=1.e2*np.finfo(float).eps
    eigval, eigvec = np.linalg.eig(Cd_gauss)            
    Q = np.matrix(eigvec)
    xdiag = np.matrix(np.diag(np.maximum(eigval, thresh)))
    Cd_gauss = Q*xdiag*Q.T
    

    if len(var)==0:
       return Cd_gauss
    else:
       Cd_gauss = Cd_gauss*np.diag(var)
    
    return Cd_gauss


#def covp_expnt(dz,Lz=1.,var=[]):
#    '''
#    Exponential 1D covariance with correlation length L
#    '''
#    layers=np.shape(dz)
#    nlyr= layers[0]
#    
##    print(dz[nlyr-1])
##    dz =np.append(dz,dz[nlyr-1])
##    z = np.append( 0., np.cumsum(dz))
##    zc = 0.5*(z[0:nlyr1-1]+z[1:nlyr1])
#    zc= centers(dz)
#    Cp_expnt=np.zeros((nlyr,nlyr))
#    for il in range(nlyr):
#        for jl in range(il,nlyr):
#            Cp_expnt[jl,il]= np.exp(-np.abs(zc[il]-zc[jl])/Lz)
#            Cp_expnt[il,jl]=Cp_expnt[jl,il]
#       
#    print(np.shape(Cp_expnt))
#    if len(var)==0:
#       return Cp_expnt
#    else:
#       Cp_expnt = Cp_expnt*np.diag(var)
#       return Cp_expnt
#
#def covp_gauss(dz,Lz=1.,var=[]):
#    '''
#    Exponential 1D covariance with correlation  length L
#    '''
#    layers=np.shape(dz)
#    nlyr= layers[0]
#    
##    print(dz[nlyr-1])
##    dz =np.append(dz,dz[nlyr-1])
##    z = np.append( 0., np.cumsum(dz))
##    zc = 0.5*(z[0:nlyr1-1]+z[1:nlyr1])
#
#    zc= centers(dz)
#    
#    
#    Cp_gauss=np.zeros((nlyr,nlyr))
#    for il in range(nlyr-1):
#       for jl in range(il,nlyr-1):
#           Cp_gauss[jl,il]= np.exp(-0.5*np.power(np.abs(zc[il]-zc[jl])/Lz,2))
#           Cp_gauss[il,jl]=Cp_gauss[jl,il]
#    #enforce positive definiteness
#    thresh=1.e2*np.finfo(float).eps
#    eigval, eigvec = np.linalg.eig(Cp_gauss)
#    Q = np.matrix(eigvec)
#    xdiag = np.matrix(np.diag(np.maximum(eigval, thresh)))
#    Cp_gauss = Q*xdiag*Q.T
#    
#    if len(var)==0:
#        return Cp_gauss
#    else:
#       Cp_gauss = Cp_gauss*np.diag(var)
#        
#        
#    return Cp_gauss
##    
#    
    
    
def prep_blockdata_aem05(filename,blocksize,writedata=False,istart=1,siteshift=1):
    '''
     Reads data from file filename and averages into blocks of size blocksize. 
    '''
    
    numlines= blocksize
    
#skips the header, reads each column in out    

    fname=filename+'.dat'
    data_obs= np.loadtxt(fname, skiprows=5)
    
    
    fline= data_obs[:,0]; 
    lonUTM= data_obs[:,1]  
    latUTM= data_obs[:,2]
    gps= data_obs[:,3]
    alt= data_obs[:,4]
    ip912= data_obs[:,5] 
    ip3005= data_obs[:,6]
    ip11962= data_obs[:,7]
    ip24510= data_obs[:,8]
    q912= data_obs[:,9]
    q3005= data_obs[:,10]
    q11962= data_obs[:,11]
    q24510= data_obs[:,12]
    
    site_x= []
    for i,d in enumerate(lonUTM[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(lonUTM[i - numlines:i])
           site_x.append(avg_numlines)
    
    site_y= []
    for i,d in enumerate(latUTM[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(latUTM[i - numlines:i])
           site_y.append(avg_numlines)
           
    gps_avg= []
    for i,d in enumerate(gps[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(gps[i - numlines:i])
           gps_avg.append(avg_numlines)
    
    alt_avg= []
    for i,d in enumerate(alt[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(alt[i - numlines:i])
           alt_avg.append(avg_numlines)

    q912_avg= []
    for i,d in enumerate(q912[:]):   
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(q912[i - numlines:i])
           q912_avg.append(avg_numlines)

    q3005_avg= []
    for i,d in enumerate(q3005[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(q3005[i - numlines:i])
           q3005_avg.append(avg_numlines)
        
    q11962_avg= []
    for i,d in enumerate(q11962[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(q11962[i - numlines:i])
           q11962_avg.append(avg_numlines)
        
    q24510_avg= []
    for i,d in enumerate(q24510[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(q24510[i - numlines:i])
           q24510_avg.append(avg_numlines)
        
    ip912_avg= []
    for i,d in enumerate(ip912[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(ip912[i - numlines:i])
           ip912_avg.append(avg_numlines)

    ip3005_avg= []
    for i,d in enumerate(ip3005[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(ip3005[i - numlines:i])
           ip3005_avg.append(avg_numlines)
        
    ip11962_avg= []
    for i,d in enumerate(ip11962[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(ip11962[i - numlines:i])
           ip11962_avg.append(avg_numlines)
        
    ip24510_avg= []
    for i,d in enumerate(ip24510[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(ip24510[i - numlines:i])
           ip24510_avg.append(avg_numlines)
        
    
    datablock=  np.column_stack((ip912_avg, ip3005_avg, ip11962_avg, ip24510_avg, q912_avg, q3005_avg, q11962_avg, q24510_avg))
    #datablock=  np.column_stack((q912_avg, q3005_avg, q11962_avg, q24510_avg,ip912_avg, ip3005_avg, ip11962_avg, ip24510_avg))
  
    #np.savetxt('data_obs_MH_block.dat', DataOut, fmt= ('%6.2f'))
    if writedata:
        np.savez_compressed('Line'+str(fline)+'_block'+str(blocksize),datablock=datablock, gps=gps_avg, alt=alt_avg, lon=site_x, lat=site_y)
     
    comment=filename+'_block'+str(blocksize)
    return datablock, gps_avg, alt_avg, site_x, site_y, comment
    
def prep_blockdata_genesis(filename,blocksize,writedata=False, istart=1,siteshift=1):
    '''
     Reads time-domain data from file filename and averages into blocks of size blocksize. 
    '''
    
    numlines= blocksize
    
    fname=filename+'.dat'
    data_obs= np.loadtxt(fname, skiprows=1 )
    
    
    fline= data_obs[:,0]; 
    lonUTM= data_obs[:,1]  
    latUTM= data_obs[:,2]
    gps= data_obs[:,3]
    alt= data_obs[:,4]
    X01= data_obs[:,5] 
    X02= data_obs[:,6]
    X03= data_obs[:,7]
    X04= data_obs[:,8]
    X05= data_obs[:,9]
    X06= data_obs[:,10]
    X07= data_obs[:,11]
    X08= data_obs[:,12]
    X09= data_obs[:,13]
    X10= data_obs[:,14]
    X11= data_obs[:,15]
    Z01= data_obs[:,16] 
    Z02= data_obs[:,17]
    Z03= data_obs[:,18]
    Z04= data_obs[:,19]
    Z05= data_obs[:,20]
    Z06= data_obs[:,21]
    Z07= data_obs[:,22]
    Z08= data_obs[:,23]
    Z09= data_obs[:,24]
    Z10= data_obs[:,25]
    Z11= data_obs[:,26]
    
    
    
    site_x= []
    for i,d in enumerate(lonUTM[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(lonUTM[i - numlines:i])
           site_x.append(avg_numlines)
    
    site_y= []
    for i,d in enumerate(latUTM[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(latUTM[i - numlines:i])
           site_y.append(avg_numlines)
           
    gps_avg= []
    for i,d in enumerate(gps[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(gps[i - numlines:i])
           gps_avg.append(avg_numlines)
    
    alt_avg= []
    for i,d in enumerate(alt[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(alt[i - numlines:i])
           alt_avg.append(avg_numlines)

    X01_avg= []
    for i,d in enumerate(X01[:]):   
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X01[i - numlines:i])
           X01_avg.append(avg_numlines)

    X02_avg= []
    for i,d in enumerate(X02[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X02[i - numlines:i])
           X02_avg.append(avg_numlines)
        
    X03_avg= []
    for i,d in enumerate(X03[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X03[i - numlines:i])
           X03_avg.append(avg_numlines)
        
    X04_avg= []
    for i,d in enumerate(X04[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X04[i - numlines:i])
           X04_avg.append(avg_numlines)
        
    X05_avg= []
    for i,d in enumerate(X05[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X05[i - numlines:i])
           X05_avg.append(avg_numlines)

    X06_avg= []
    for i,d in enumerate(X06[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(X06[i - numlines:i])
          X06_avg.append(avg_numlines)
        
    X07_avg= []
    for i,d in enumerate(X07[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(X07[i - numlines:i])
           X07_avg.append(avg_numlines)

    X08_avg= []
    for i,d in enumerate(X08[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(X08[i - numlines:i])
          X08_avg.append(avg_numlines)
          
    X09_avg= []
    for i,d in enumerate(X09[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(X09[i - numlines:i])
          X09_avg.append(avg_numlines)
          
    X10_avg= []
    for i,d in enumerate(X10[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(X10[i - numlines:i])
          X10_avg.append(avg_numlines)
          
    X11_avg= []
    for i,d in enumerate(X11[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(X11[i - numlines:i])
          X11_avg.append(avg_numlines)
   

    Z01_avg= []
    for i,d in enumerate(Z01[:]):   
       if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z01[i - numlines:i])
          Z01_avg.append(avg_numlines)

    Z02_avg= []
    for i,d in enumerate(Z02[:]):
       if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z02[i - numlines:i])
          Z02_avg.append(avg_numlines)
        
    Z03_avg= []
    for i,d in enumerate(Z03[:]):
       if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z03[i - numlines:i])
          Z03_avg.append(avg_numlines)
        
    Z04_avg= []
    for i,d in enumerate(Z04[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(Z04[i - numlines:i])
           Z04_avg.append(avg_numlines)
        
    Z05_avg= []
    for i,d in enumerate(Z05[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(Z05[i - numlines:i])
           Z05_avg.append(avg_numlines)

    Z06_avg= []
    for i,d in enumerate(Z06[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z06[i - numlines:i])
          Z06_avg.append(avg_numlines)
        
    Z07_avg= []
    for i,d in enumerate(Z07[:]):
       if (i % numlines) == 0 and i != 0:
           avg_numlines= np.mean(Z07[i - numlines:i])
           Z07_avg.append(avg_numlines)

    Z08_avg= []
    for i,d in enumerate(Z08[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z08[i - numlines:i])
          Z08_avg.append(avg_numlines)
          
    Z09_avg= []
    for i,d in enumerate(Z09[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z09[i - numlines:i])
          Z09_avg.append(avg_numlines)
          
    Z10_avg= []
    for i,d in enumerate(Z10[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z10[i - numlines:i])
          Z10_avg.append(avg_numlines)
          
    Z11_avg= []
    for i,d in enumerate(Z11[:]):
      if (i % numlines) == 0 and i != 0:
          avg_numlines= np.mean(Z11[i - numlines:i])
          Z11_avg.append(avg_numlines)


    datablock=  np.column_stack((X01_avg, X02_avg, X03_avg, X04_avg, X05_avg, X06_avg, X07_avg, X08_avg, X09_avg, X10_avg, X11_avg,\
    Z01_avg, Z02_avg, Z03_avg, Z04_avg, Z05_avg, Z06_avg, Z07_avg, Z08_avg, Z09_avg, Z10_avg, Z11_avg))
  
    #np.savetxt('data_obs_MH_block.dat', DataOut, fmt= ('%6.2f'))
    if writedata:
        np.savez_compressed('Line'+str(fline)+'_block'+str(blocksize),datablock=datablock, gps=gps_avg, alt=alt_avg, lon=site_x, lat=site_y)
     
    comment=filename+'_block'+str(blocksize)
    return datablock, gps_avg, alt_avg, site_x, site_y, comment    
    
def prep_insert_flag(M, criterion='neg',threshval=0., columns=[0,0],flag=None,incol=None):
    '''
     replaces bad values in data by nans, and returns a boolean array
     for bad values
    '''
    if flag==None: flag= np.nan
    
    cols= range(columns[0],columns[1]+1)
    T=M[:,cols]

    oldflags = np.where(np.isnan(T))

    if criterion.lower()[0:3]   ==  'neg':
        T[oldflags]= -9999.
        T[T<=0] = flag
    elif criterion.lower()[0:3] == 'les':
        T[oldflags]= threshval-9999.
        T[T<=threshval] = flag
    elif criterion.lower()[0:3] == 'gre':
        T[oldflags]= threshval+9999.
        T[T>=threshval] = flag
    elif criterion.lower()[0:3] == 'plm': 
        T[oldflags]= threshval+9999.
        T[np.abs(T) > threshval] = flag 
    elif criterion.lower()[0:3] == 'nan':  
        if incol == None:
           error(' no column  given!')
        newindex = np.isnan(M[:,incol])
        T[newindex,:] = flag
        
    T[oldflags] =  flag
    M[:,cols]=T
    NaNindex=np.where(np.isnan(M))    
    
    return M, NaNindex
      
def prep_handle_gaps(M,columns=[0,0],impute_val='delete',stdval=None):
    '''
     replaces bad values in data by nans, and returns a boolean array
     for bad values
    '''
    cols= range(columns[0],columns[1]+1)
    T=M[:,cols] 
 
    sizeT=np.shape(T)
   
    nanindex = np.where(np.isnan(T))
    nrows = np.ravel(nanindex[0])
    ncols = np.ravel(nanindex[1])
    if impute_val.lower()[0:3]   == 'del':
         T = T[~np.isnan(M).any(axis=1)]
         M = M[~np.isnan(M).any(axis=1)]
    if impute_val.lower()[0:3]   == 'ave':
        val=np.nanmean(T,0)
        for irow in nrows:
          for icol in ncols:
             T[irow,icol]=val[icol] 
    elif impute_val.lower()[0:3] == 'med':
        val=np.nanmedian(T,0)
        for irow in nrows:
          for icol in ncols:
             T[irow,icol]=val[icol] 
    elif impute_val.lower()[0:3] == 'noi':
        v_avg=np.nanmean(T,0)
        if stdval == None: 
           v_std=np.nanstd(T,0)
        else:
           v_std=stdval
        val=v_avg+v_std*np.random.randn(sizeT[0],sizeT[1])
        for irow in nrows:
          for icol in ncols:
              T[irow,icol]=val[irow,icol] 
              
    M[:,cols]=T
     
    NaNindex=np.where(np.isnan(M))    
    # NaNindex=np.isnan(M)

    return M, NaNindex
    
def prep_interpolate_gaps(M, method='cubic',param=0.):
    '''
         interpolates bad values or nans in M, which is a block of data in the 
     internal aempy format the first 4 vars are the a,y,alt, gps  
     method: (str or int) specifies the kind of interpolation as a string 
     (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’) where 
     ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline interpolation of 
     first, second or third order) or as an integer specifying the order 
     of the spline interpolator to use. Default is ‘linear’.
    '''
    #from  numpy import nan 
    import scipy.interpolate as si
    sizeM=np.shape(M)
    nvars=sizeM[1]

    yd= M[:,1]  
    xd= M[:,2]
    rd = np.sqrt(xd*xd+yd*yd)
    for var in range [5,nvars]:
          var = M[:,var]

          if not np.all(np.isnan(var)):
            r0 = rd[not np.isnan(var)]  
            v0 = var[not np.isnan(var)]   
            vf = si.interpolate(r0,v0,kind=method)
            M[:,var]=vf(rd)

    return M 


def prep_pcafilter(D, k, columns = None, out_full =False, out = False, thresh_MSE=0.):
    '''
     PCA filter as used by Minsley et al.(2012)
     IN:
        M block of data (usually flightline)
        k number of PCs < ndata
    '''
    
    if columns == None:
        error('pcafilter: no columns definded! ') 
           
    cols= range(columns[0],columns[1]+1)

    if k > np.size(cols): k = np.size(cols)

    P = D.copy()
    T = D[:,cols] 
   
    T_avg=np.nanmean(T,0)
       # T_std=np.nanvar(T,0)
        
    X = (T-T_avg ) #/M_std
    #print(np.mean(X,1))    
    #U,s,Vt = np.linalg.svd(X, full_matrices=False)
    #print(T)
    U,S,Vt = la.svd(X, full_matrices=False)
    #print(s)
    SS = np.diag(S)
    V = Vt.T

    #Tk = np.dot(U[:, :k], np.dot(S[:k, :k], V[:,:k].T))
    Tk = np.dot(U[:, :k], np.dot(SS[:k, :k], V[:,:k].T))
    
    MSE = np.sqrt(np.mean((X - Tk)**2))
    #MSE = (np.square(X - Tk)).mean(axis=None)
    if out: print (" using "+str(k)+" PCs, MSE = %.6G" %(MSE))
   
    P[:,cols]=Tk+T_avg
    
        
    if out_full:
        return P, U, S, V, MSE
    else:    
        return P
    


def prep_blockdata_aem05_v1(filename,blocksize,writedata=False, istart=1,siteshift=1):
    
    window_width = blocksize   
    
    fname    = filename+'.dat'
    #data_obs = np.loadtxt(fname, skiprows=5)
    data_obs = np.loadtxt(fname, skiprows=1)

    
    fline    = data_obs[:,0]; 
    lonUTM   = data_obs[:,1]  
    latUTM   = data_obs[:,2]
    gps      = data_obs[:,3]
    alt      = data_obs[:,4]
    ip912    = data_obs[:,5] 
    ip3005   = data_obs[:,6]
    ip11962  = data_obs[:,7]
    ip24510  = data_obs[:,8]
    q912     = data_obs[:,9]
    q3005    = data_obs[:,10]
    q11962   = data_obs[:,11]
    q24510   = data_obs[:,12]



    lonUTM_cumsumvec        = np.cumsum(np.insert(lonUTM,0,0))
    site_x                  = (lonUTM_cumsumvec[window_width:] - lonUTM_cumsumvec[:-window_width]) / window_width

    latUTM_cumsumvec        = np.cumsum(np.insert(latUTM,0,0))
    site_y                  = (latUTM_cumsumvec[window_width:] - latUTM_cumsumvec[:-window_width]) / window_width

    gps_cumsumvec           = np.cumsum(np.insert(gps,0,0))
    gps_avg                 = (gps_cumsumvec[window_width:] - gps_cumsumvec[:-window_width]) / window_width
    
    alt_cumsumvec           = np.cumsum(np.insert(alt,0,0))
    alt_avg                 = (alt_cumsumvec[window_width:] - alt_cumsumvec[:-window_width]) / window_width

    ip912_cumsumvec         = np.cumsum(np.insert(ip912,0,0)) 
    ip912_avg               = (ip912_cumsumvec[window_width:] - ip912_cumsumvec[:-window_width]) / window_width
     
    ip3005_cumsumvec        = np.cumsum(np.insert(ip3005,0,0)) 
    ip3005_avg              = (ip3005_cumsumvec[window_width:] - ip3005_cumsumvec[:-window_width]) / window_width

    ip11962_cumsumvec       = np.cumsum(np.insert(ip11962,0,0)) 
    ip11962_avg             = (ip11962_cumsumvec[window_width:] - ip11962_cumsumvec[:-window_width]) / window_width
     
    ip24510_cumsumvec       = np.cumsum(np.insert(ip24510,0,0)) 
    ip24510_avg             = (ip24510_cumsumvec[window_width:] - ip24510_cumsumvec[:-window_width]) / window_width
     
    q912_cumsumvec          = np.cumsum(np.insert(q912,0,0)) 
    q912_avg                = (q912_cumsumvec[window_width:] - q912_cumsumvec[:-window_width]) / window_width
      
    q3005_cumsumvec         = np.cumsum(np.insert(q3005,0,0)) 
    q3005_avg               = (q3005_cumsumvec[window_width:] - q3005_cumsumvec[:-window_width]) / window_width 

    q11962_cumsumvec        = np.cumsum(np.insert(q11962,0,0)) 
    q11962_avg              = (q11962_cumsumvec[window_width:] - q11962_cumsumvec[:-window_width]) / window_width
     
    q24510_cumsumvec        = np.cumsum(np.insert(q24510,0,0)) 
    q24510_avg              = (q24510_cumsumvec[window_width:] - q24510_cumsumvec[:-window_width]) / window_width
     
        
    
    datablock=  np.column_stack((ip912_avg, ip3005_avg, ip11962_avg, ip24510_avg, q912_avg, q3005_avg, q11962_avg, q24510_avg))
    
    fileout=filename+'_Line'+str(fline)+'_block'+str(blocksize)
    
    if writedata:
        np.savez_compressed(fileout,datablock=datablock, gps=gps_avg, alt= alt_avg,\
                            lon=site_x, lat=site_y)
     
    comment=fileout
    return datablock, gps_avg, alt_avg, site_x, site_y, comment
        
    
def get_filelist(mystr='*',mypath='.'):
    """
    geberates filelist from path and wildcard
    author: vrath
    """
    import sys
    import os
    import fnmatch
    
    filelist = fnmatch.filter(os.listdir(mypath), mystr)
    return filelist
    
    

           
    
def prep_blocksdata(filename,blocksize,out='med',istart=1,siteshift=1, writedata=False):
    '''
     Reads data from file filename and averages into blocks of size blocksize. 
    '''
    
    if (blocksize % 2 == 0):
        hwidth=blocksize/2+1
        blocksize=blocksize+1
    else:
        hwidth =(blocksize-1)/2 +1

    data_obs= np.loadtxt(filename+'.dat', skiprows=1)
    fline= data_obs[0,0]
    work=data_obs[:,1:]
    sizework=np.shape(work)
    
    work_val = []

    for isite in np.arange(1,sizework[1]+1,siteshift):
       i1=np.maximum(0,isite-hwidth)
       i2=np.minimum(sizework[1],isite+hwidth)
       nds=i2-i1
 
   
       # var_block= np.nanvar(work[i1:i2,:],axis = 0)
#       #percentiles 15.9/84.1 (1s)  2.3/97.7 (2s)
#       pm_block= np.nanpercentile(work[i1:i2,:],15.9,axis = 0)
#       pp_block= np.nanpercentile(work[i1:i2,:],84.1,axis = 0)
      
       if out.lower[0:1] == 'block': 
         work_val.append(np.ravel(work[i1:i2,:],0))
       elif out.lower[0:1] == 'ave':
         val= np.nanmean(work[i1:i2,:],axis = 0)
         work_val.append(val)
       elif out.lower[0:1] == 'med':
         val= np.nanmedian(work[i1:i2,:],axis = 0)
         work_val.append(val)
    data_block = work_val
    
    fileout =filename+'_block_'+out+str(blocksize)+'_shift'+str(siteshift)
    if comment == None:
       comment = fileout
    if writedata:
        np.savez_compressed(fileout, data_block=data_block, comment=None)
     

    
    return data_block, comment
