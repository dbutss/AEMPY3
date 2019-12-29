#! /usr/bin/python
import sys 
import os
sys.path.insert(1,os.path.join(os.path.dirname(__file__), "/local/volker/ElectroMagnetics/AirborneEM/aem/AEMC3/modules"))

import core1d as td
import numpy as np
#import pylab
import matplotlib.pyplot as plt


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
#                    time JT, station JS. (nT)MAP
#     BFD(JF,JS,2) - the transverse component
#     BFD(JF,JS,3) - the vertical component
#
#------------------------------------------------------------------------
#!  SIGN CONVENTION:
#  The normal layered earth field coordinate system used in this
#  subroutine has X (JC=1) positive along the flight path, Y (JC=2)
#  positive to starboard, and Z (JC=3) positive down.
#------------------------------------------------------------------------
# LAST CHANGE:  28 August  2016   DK  
#------------------------------------------------------------------------

nan = float('NaN')
#default values for input parameter
mode= 0
sz = 100
#nlyr = 5
nlyr = 3
#thk = [30.0, 20.0, 44.11, nan]
#res = [10.0, 1000.0, 10.0, 300.0]; res = np.log10(res);
#thk = np.array([10.0, 5.0, 20.0, 25.0])
thk = np.array([15.0, 30.0, nan])
#res = np.array([100.0, 10.0, 100.0, 1.0, 100.0]) 
res = np.array([50.0, 300.0, 10.0]); 
# res = np.log10(res);
#res = np.array([50.0, 300.0, 10.0])
znlyr = np.zeros(nlyr)
onlyr = np.ones(nlyr)
reps = onlyr
rmu = onlyr
calf = onlyr
ctau = znlyr
cfreq = onlyr
calc_data = np.zeros(22)

m_true = np.concatenate([res,reps,rmu,calf,ctau,cfreq,thk])
m_true = np.array([m_true])

calc_data = td.aemfwd1d_aem05(mode,sz,nlyr,m_true[:,:-1])
print "Results from python wrapper" % (calc_data)



# plot data
plt.clf()

freq = [912.0, 3005.0, 11962.0, 24510.0];
derr=30.0

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.errorbar(freq, calc_data[0:4], yerr=derr, fmt='or', markersize=9, label ='In Phase')
plt.errorbar(freq, calc_data[4:8], yerr=derr, fmt='ob', markersize=9, label ='Quadrature')
plt.xscale('log')
plt.legend(loc='upper left', numpoints=1)
plt.xlabel('frequency (Hz)', fontsize=20)
plt.ylabel('data (ppm)', fontsize=20)
ax1.tick_params(labelsize=15)

plt.savefig('data_truemodel_#1', format='png')
#plt.grid(True)
plt.show()

