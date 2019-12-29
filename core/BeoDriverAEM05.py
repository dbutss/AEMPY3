#! /usr/bin/python

import core1d as fd
import numpy as np

"""
PROGRAM beodriverfd

This programme is the equivalent of BeoDriverGTK.f90. 
"""
#----------------------------------------------------------------------
#  Uses AEM1D_FD to compute  the frequency-domain layered earth 
#  H field for a dipole of unit moment and current.
#
#                             INPUT
#                             -----
#        JS = station reference
#      FREQ - array of NFRQ frequencies
#     TXCLN - angle in radians that TX dipole makes with vertical (climb = +)
#     TXA90 - true for vertical co-planar briadside array
#     NSTAT - number of stations in survey line.
#        SZ - array of transmitter altitudes
#       ZRX - vertical offset of each receiver from transmitter  (below = +)
#       XRX - in-line horizontal offset of RX J;  (behind = +)
#       YRX - transverse horizontal offset of RX J (port = +)
#      NLYR - number of layers
#       RES - layer resistivities
#      REPS - array of relative dislectric constants
#      RMUX - mu(i) / mu(0)
#       THK - array of layer thicknesses
#     CALF, CTAU, CFREQ are the layered earth Cole-Cole parameters.
#
#                             OUTPUT
#                             ------
#     BFD(JF,JS,1) - the in-line component of the layered earth response at
#                    time JT, station JS. (nT)
#     BFD(JF,JS,2) - the transverse component
#     BFD(JF,JS,3) - the vertical component
#
#------------------------------------------------------------------------
#!  SIGN CONVENTION:
#  The normal layered earth field coordinate system used in this
#  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
#  positive to starboard, and Z (JC=3) positive down.
#------------------------------------------------------------------------
# LAST CHANGE:  11 July  2016   DK  
#------------------------------------------------------------------------
#default values for input parameter
aem_system = 'aem05'
print ('\n  AEM forward modeling for '+ aem_system+' system\n \n') 

alt = 60.
#nlyr = 5
nlyr = 5
#  generate synthetic data
znlyr = np.zeros(nlyr);iznlyr=znlyr.astype(int)
onlyr = np.ones(nlyr) ;ionlyr=onlyr.astype(int)
isactive = np.concatenate([ionlyr,iznlyr,iznlyr,iznlyr,iznlyr,iznlyr,iznlyr[:len(iznlyr)-1]])
isactive = np.array([isactive])
sizepar  = np.shape(isactive); #sizepar = sizepar[:]
sizeact  = np.shape(np.flatnonzero(isactive))
mpara=sizeact[0];

thk = np.array([10. , 15., 20., 25 ])   
res = np.array([100., 10., 100.,1. ,100.])
#thk = np.array([10.0, 5.0, 20.0, 25.0])
#thk = np.array([30.0, 20.0, 44.11])
#res = np.array([100.0, 10.0, 100.0, 1.0, 100.0]) 
#res = np.array([10.0, 1000.0, 10.0, 300.0])
reps = np.ones(nlyr)
rmu = np.ones(nlyr)
calf = np.ones(nlyr) 
ctau = np.zeros(nlyr)
cfreq = np.ones(nlyr)


model =np.ones((1,7*nlyr))
model[0,0*nlyr:1*nlyr] = res
model[0,1*nlyr:2*nlyr] = reps
model[0,2*nlyr:3*nlyr] = rmu
model[0,3*nlyr:4*nlyr] = calf
model[0,4*nlyr:5*nlyr] = ctau
model[0,5*nlyr:6*nlyr] = cfreq
model[0,6*nlyr:7*nlyr-1]= thk

calc_data = np.zeros(8)

calc_data = fd.aemfwd1d_aem05(mode,sz,nlyr,model)


calc_data, jac = fd.aemjac1d_aem05(js,sz,nlyr,model,mpara,isactive)

#print "Results from python wrapper %6.2f" % (calc_data)

calc_data
