#!/usr/bin/env python2

"""
Created on Fri Oct 14 14:18:46 2016

@author: VR/DK/RD

"""
import sys
import os
import fnmatch
import mmap
import struct
#from mpl_toolkits.basemap import Basemap
#import mpl_toolkits.basemap.pyproj as pyproj
import pyproj
from sys import exit as error

#import math  as ma
#import numpy as np
#import scipy as sc       


def get_files(SearchString=None, SearchDirectory='.'):
       """ 
       FileList = get_files(Filterstring) produces a list
       of files frm a searchstring (allows wildcars)
             
       """
       FileList = fnmatch.filter(os.listdir(SearchDirectory), SearchString)
       
       return FileList
   
    


def process_aem05_data(InputFile=None,FmtCtrl=None, OutputFile=None, FColumns=None, Header=None):
    """
     transforms aem05 data to aempy-compatible format. 
     takes 3 arguments: 
    
     @1 is the input file name (with extendion)
     @2 is the subformat (A1 = A1 survey, TF = Overlap lines, 
        TB = Tellus Border, TE =Tellus)
     @3 is the output file name 
     
    """
    
#    import sys
#    from sys import exit as error
#    import os
#    import numpy as np
    #input Tellus project
    
    ierr =0 
    
    if InputFile==None:
        error('No Input file given!')
        
    if FmtCtrl==None:
        if FColumns==None:
            error('No Format given!')

    
    filext = '.dat'    
    if OutputFile==None:
        filen,filext0 = os.path.splitext(InputFile)
    else:
        filen,filext0 = os.path.splitext(OutputFile)
 
    OutAll=filen+'_All'+filext
    OutSub=filen+'_Sub'
           
    if Header==None:
        header='# aempy aem05 Format\n'
        
  

#    inputF = InputFile+'.dat'
    inputF = InputFile
    inFile = open(inputF,'r')
    # +'.dat'
    outFileAll = open(OutAll,'w')
    # output one line at a time, initially dummy val
    subF=OutSub
    outFileSub = open(subF,'w')
    
#    print(OutAll)
#    print(OutSub)
 
    if FmtCtrl==None:
        ucolumns=FColumns
        linecol=0
        
    elif  FmtCtrl == 'A1' or FmtCtrl == 'A2' or FmtCtrl == 'W1':
        '''
       Translation table from .xyz files to internal data format:
       
   0        0 LINE           08           -        -         Line number - LLLL.SR (L=line, S=segment, R=reflight)
          1 FLIGHT         06           -        -         Flight number    
           2  DATE           10           -        -         Date YYYYMMDD
           3  DAY            05           -        -         Day of year
           4  TIME           11           s        -         UTC seconds
           5  ITM-X          13           m        *         X coordinate, IRENET95 ITM
           6  ITM-Y          13           m        *         Y coordinate, IRENET95 ITM
   1        7  UTM-X          13           m        *         X coordinate, WGS-84 UTM 29N
   2        8  UTM-Y          13           m        *         Y coordinate, WGS-84 UTM 29N
   3        9  WGSHGT         13           m        *         GPS Elevation (above WGS-84 Ellipsoid)
   4       10 RADAR          11           m        *         Radar Altimeter Clearance above Terrai4n
       11 P09ppm         08          ppm       *         In-phase 912 Hz
          12 Q09ppm         08          ppm       *         Quadrature 912 Hz
          13 P3ppm          08          ppm       *         In-phase 3005 Hz
          14 Q3ppm          08          ppm       *         Quadrature 3005 Hz
          15 P12ppm         08          ppm       *         In-phase 11962 Hz
          16 Q12ppm         08          ppm       *         Quadrature 11962 Hz
          17 P25ppm         08          ppm       *         In-phase 24510 Hz
          18 Q25ppm         08          ppm       *         Quadrature 24510 Hz
          19 P09filt        08          ppm       *         Filtered in-phase 912 Hz
          20 Q09filt        08          ppm       *         Filtered quadrature 912 Hz
          21 P3filt         08          ppm       *         Filtered in-phase 3005 Hz
          22 Q3filt         08          ppm       *         Filtered quadrature 3005 Hz
          23 P12filt        08          ppm       *         Filtered in-phase 11962 Hz
          24 Q12filt        08          ppm       *         Filtered quadrature 11962 Hz
          25 P25filt        08          ppm       *         Filtered in-phase 24510 Hz
          26 Q25filt        08          ppm       *         Filtered quadrature 24510 Hz
   5       27 P09lev         08          ppm       *         Levelled and filtered in-phase 912 Hz
   9       28 Q09lev         08          ppm       *         Levelled and filtered quadrature 912 Hz
   6       29 P3lev          08          ppm       *         Levelled and filtered in-phase 3005 Hz
   10  30 Q3lev          08          ppm       *         Levelled and filtered quadrature 3005 Hz
   7       31 P12lev         08          ppm       *         Levelled and filtered in-phase 11962 Hz
   11  32 Q12lev         08          ppm       *         Levelled and filtered quadrature 11962 Hz
   8       33 P25lev         08          ppm       *         Levelled and filtered in-phase 24510 Hz
   12  34 Q25lev         08          ppm       *         Levelled and filtered quadrature 24510 Hz
       35 Radio_Flag     09           -        *         Radio call flag
   13  36 PLM_mV         10           mV       *         Power line monitor
       37 Res09          10         ohm-m      *         Apparent resistivity, half-space model, 912 Hz
          38 Res3           10         ohm-m      *         Apparent resistivity, half-space model, 3005 Hz
          39 Res12          10         ohm-m      *         Apparent resistivity, half-space model, 11962 Hz
          40 Res25          10         ohm-m      *         Apparent resistivity, half-space model, 24510 Hz
          41 Res09_MLEV     10         ohm-m      *         Microlevelled apparent resistivity, half-space model, 912 Hz
          42 Res3_MLEV      10         ohm-m      *         Microlevelled apparent resistivity, half-space model, 3005 Hz
          43 Res12_MLEV     10         ohm-m      *         Microlevelled apparent resistivity, half-space model, 11962 Hz
          44 Res25_MLEV     10         ohm-m      *         Microlevelled apparent resistivity, half-space model, 24510 Hz
          45 Depth09        10           m        *         Centroid depth 912Hz
          46 Depth3         10           m        *         Centroid depth 3005Hz
          47 Depth12        10           m        *         Centroid depth 11962Hz
          48 Depth25        10           m        *         Centroid depth 24510Hz
          49 ResSlice10     11         ohm-m      *         Resistivity depth slice at 10m
          50 ResSlice30     11         ohm-m      *         Resistivity depth slice at 30m
          51 ResSlice60     11         ohm-m      *         Resistivity depth slice at 60m
          52 ResSlice100    11         ohm-m      *         Resistivity depth slice at 100m
          53 ResSlice10_ML  11         ohm-m      *         Microlevelled resistivity depth slice at 10m
          54 ResSlice30_ML  14         ohm-m      *         Microlevelled resistivity depth slice at 30m
          55 ResSlice60_ML  14         ohm-m      *         Microlevelled resistivity depth slice at 60m
          56 ResSlice100_ML 14         ohm-m      *         Microlevelled resistivity depth slice at 100m
   
        
    DEM for plotting  can be derived from radar altimetry (c10) and GPS (c09) as: c09-c10
        
        '''
        ucolumns = [ 7, 8, 9, 10, 27, 29, 31, 33, 28, 30 ,32, 34, 36] # 
        linecol=0
    elif FmtCtrl == 'TF':
        error('Data Format'+FmtCtrl+' not yet adapted !')
        ucolumns = [7, 6, 10, 19, 20, 21, 22, 23, 24 ,25, 26]
        linecol=1
    elif FmtCtrl == 'TB':
        error('Data Format'+FmtCtrl+' not yet adapted !')
        ucolumns = [7, 6, 11, 19, 20, 21, 22, 23, 24 ,25, 26]                
        linecol=1
    elif FmtCtrl == 'T4':  
        error('Data Format'+FmtCtrl+' not yet adapted !')
        ucolumns = [7, 6, 10, 19, 20, 21, 22, 23, 24 ,25, 26]
        linecol=1
    
    else:
        ierr = 2        
        error('Data Format'+FmtCtrl+' not implemented!')

    
    outFileAll.write(header)

    with inFile as f:
        for line in f:
            if line[0] == '/':
                continue # skip header lines
            elif line[0:3] == 'TIE':
                outFileSub.close()
                lineN = int(float(line[4:].strip()))
                subF = OutSub+'_T' + str(lineN) + filext
                outFileSub = open(subF,'w')
                outFileSub.write(header)
                continue
            elif line[0:4] == 'LINE':
                outFileSub.close()
                # dummyVar...#    from sys import exit as error

                lineN = int(float(line[5:].strip()))
                subF = OutSub+'_L' + str(lineN) + filext
                outFileSub.close()
                outFileSub = open(subF,'w')
                outFileSub.write(header)
                continue
            # should be in data body by this point
            ln=filter(None,line.strip().split(' '))
            
            #print(ln)
            fline=''
            fline=fline + ' ' + str(int(float(ln[linecol])))
            try:
                for ind in ucolumns:
                    if ln[ind] == '*':
                        fline=fline + ' ' + str(float('NaN'))
                    else:
                        fline=fline + ' ' + str(float(ln[ind]))
            except IndexError:
                ierr = 1
                print ('Data error, record skipped in line: ', lineN)
                
            outFileAll.write(fline)
            outFileAll.write('\n')
            outFileSub.write(fline)
            outFileSub.write('\n')

    outFileSub.close()
    outFileAll.close()
    inFile.close()     
    
    return ierr

def process_genesis_data(InputFile=None,FmtCtrl=None,OutputFile=None,FColumns=None,Header=None):

    """
     tranforms CGG Genesis to aempy format 
     @1 is the file name  
     @2 is the subformat (default = NM)
     @3 is the output file name 
     
    """

#    import sys
#    import    os 
#    import numpy as np

    ierr = 0

    if InputFile==None:
        error('No Input file given!')
        
    if FmtCtrl==None:
        if FColumns==None:
            error('No Format given!')
            
            
    filext = '.dat'    
    if OutputFile==None:
        filen,filext0 = os.path.splitext(InputFile)
    else:
        filen,filext0 = os.path.splitext(OutputFile)
 
    OutAll=filen+'_All'+filext
    OutSub=filen+'_Sub'

            
    if Header==None:
        header='# aempy CGG Genesis Format\n'


    inputF = InputFile
    inFile = open(inputF,'r')

    outFileAll = open(OutAll,'w')
    # output one line at a time, initially dummy val
    subF=OutAll
    outFileSub = open(OutSub,'w')
    
        
    
    if FmtCtrl==None:
        ucolumns=FColumns
        linecol=0

    elif FmtCtrl == 'NM':
        ucolumns = [7,8,11,12,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82]
        linecol=0
   
    else:
        ierr= 2
        error('Data Format'+FmtCtrl+' not implemented!')
    
    


    
    outFileAll.write(header)
    with inputF as f:
        for line in f:
            if line[0] == '/':
                continue # skip header lines
            elif line[0:3] == 'Tie':
                outFileSub.close()
                #print line
                lineN = int(float(line[4:].strip()))
                subF = OutSub+'_T' + str(lineN) + filext
                outFileSub = open(subF,'w')
                outFileSub.write(header)
                continue
            elif line[0:4] == 'Line':
                outFileSub.close()
                # dummyVar...
                lineN = int(float(line[5:].strip()))
                subF = OutSub+'_L' + str(lineN) + filext
                outFileSub = open(subF,'w')
                outFileSub.write(header)

                continue
            
            # should be in data body by this point
            
            ln=filter(None,line.strip().split(' '))
     
            fline=''
            fline=fline + ' ' + str(int(float(ln[linecol])))
            # print 'fline', fline
            try:
                for ind in ucolumns:
                    if ln[ind] == '*':
                        fline=fline + ' ' + str(float('NaN'))
                    else:
                        fline=fline + ' ' + str(float(ln[ind]))
            except IndexError:
                ierr = 1
                print ('Data error, record skipped in line: ', lineN)
                
            outFileAll.write(fline)
            outFileAll.write('\n')
            outFileSub.write(fline)
            outFileSub.write('\n')
            
            
# now, to exchange dummyVar
    outFileSub.close()
    outFileAll.close()    
    inFile.close()    
    return ierr
    
def choose_data_rect(InputFile=None,Corners=None, OutputFile=None):
    """
     Chooses rectangular area frm aemc3-data format file. 
     takes 3 arguments: 
    
     @1 is the input file name 
     @2 are the left lower and right uper corners in m as [minX maxX minY maxY]
     @3 is the output file name 
     
    """
    
#    from sys import exit as error
#    import os
#    import numpy as np
    
    if Corners==None:
        error('No Rectangle given!')
    Corners=filter(None,Corners.strip().split(' '))
    
    filext = '.dat'    
    if OutputFile==None:
        filen,filext0 = os.path.splitext(InputFile)
    else:
        filen,filext0 = os.path.splitext(OutputFile)

    header = '# '+InputFile+' aempy format area: '+str(Corners)+'\n'
    print (header)
    
#    inputF = InputFile+'.dat'
    inputF = InputFile
    inFile = open(inputF,'r')

    outputF=filen+'_Rectangle'+filext
    outFile = open(outputF,'w')
    
    lInd=0
    with inFile as f:
        for line in f:
            if lInd < 1:
                if lInd == 2:
                    outFile.write(header)
                    lInd += 1
                    continue
                else:
                    outFile.write(line)
                    lInd += 1
                    continue
            ln=filter(None,line.strip().split(' '))
            
            #print(ln[1],Corners[2],Corners[3])
            if ln[1] > Corners[0] and ln[1] < Corners[1]:
                if ln[2] > Corners[2] and ln[2] < Corners[3]:
                    outFile.write(line)
                    continue
            lInd += 1
            
    
    outFile.close()
#    nL=sum(1 for lineC in open(outFile.name))-1
#    print 'number of lines: ', nL
    inFile.close()


def point_inside_polygon(x,y,poly):
    '''
    Determine if a point is inside a given polygon or not, where
    THE Polygon is given as a list of (x,y) pairs.
    Returns True  when point (x,y) ins inside polygon poly, False otherwise
    '''
    n = len(poly)
    inside =False
    
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside
    
def proj_to_line(x,y,line):
    '''
    Projects a point onto a line, where line is represented by two arbitrary 
    points. as an array
    '''
#    http://www.vcskicks.com/code-snippet/point-projection.php
#    private Point Project(Point line1, Point line2, Point toProject)
#{
#    double m = (double)(line2.Y - line1.Y) / (line2.X - line1.X);
#    double b = (double)line1.Y - (m * line1.X);
#
#    double x = (m * toProject.Y + toProject.X - m * b) / (m * m + 1);
#    double y = (m * m * toProject.Y + m * toProject.X + b) / (m * m + 1);
#
#    return new Point((int)x, (int)y);
#}
    x1=line[0,0]
    x2=line[1,0]
    y1=line[0,1]
    y2=line[1,1]
    m = (y2 - y1) / (x2 - x1);
    b = y1 - (m * x1);
#
    xn = (m * y + x - m * b) / (m * m + 1.);
    yn = (m * m * y + m * x + b) / (m * m + 1.);
#
    return xn, yn
#}


"""    
# This file is mostly a straight translation of
# GeographicLib/src/Geoid.cpp from C++ to Python
# by Kim Vandry <vandry@TZoNE.ORG>
#
# /**
# * \file Geoid.cpp
# * \brief Implementation for GeographicLib::Geoid class
# *
# * Copyright (c) Charles Karney (2009) <charles@karney.com>
# * and licensed under the LGPL.  For more information, see
# * http://geographiclib.sourceforge.net/
# **********************************************************************/
#
# Geoid height grade not supported

"""
class GeoidBadDataFile(Exception):
    pass

class GeoidHeight(object):
    """Calculate the height of the WGS84 geoid above the
    ellipsoid at any given latitude and longitude

    :param name: name to PGM file containing model info
    download from http://geographiclib.sourceforge.net/1.18/geoid.html
    """
    c0 = 240
    c3 = (
    (  9, -18, -88,    0,  96,   90,   0,   0, -60, -20),
    ( -9,  18,   8,    0, -96,   30,   0,   0,  60, -20),
    (  9, -88, -18,   90,  96,    0, -20, -60,   0,   0),
    (186, -42, -42, -150, -96, -150,  60,  60,  60,  60),
    ( 54, 162, -78,   30, -24,  -90, -60,  60, -60,  60),
    ( -9, -32,  18,   30,  24,    0,  20, -60,   0,   0),
    ( -9,   8,  18,   30, -96,    0, -20,  60,   0,   0),
    ( 54, -78, 162,  -90, -24,   30,  60, -60,  60, -60),
    (-54,  78,  78,   90, 144,   90, -60, -60, -60, -60),
    (  9,  -8, -18,  -30, -24,    0,  20,  60,   0,   0),
    ( -9,  18, -32,    0,  24,   30,   0,   0, -60,  20),
    (  9, -18,  -8,    0, -24,  -30,   0,   0,  60,  20),
    )

    c0n = 372
    c3n = (
    (  0, 0, -131, 0,  138,  144, 0,   0, -102, -31),
    (  0, 0,    7, 0, -138,   42, 0,   0,  102, -31),
    ( 62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31),
    (124, 0,  -62, 0,    0, -124, 0,   0,    0,  62),
    (124, 0,  -62, 0,    0, -124, 0,   0,    0,  62),
    ( 62, 0,  -31, 0,    0,  -62, 0,   0,    0,  31),
    (  0, 0,   45, 0, -183,   -9, 0,  93,   18,   0),
    (  0, 0,  216, 0,   33,   87, 0, -93,   12, -93),
    (  0, 0,  156, 0,  153,   99, 0, -93,  -12, -93),
    (  0, 0,  -45, 0,   -3,    9, 0,  93,  -18,   0),
    (  0, 0,  -55, 0,   48,   42, 0,   0,  -84,  31),
    (  0, 0,   -7, 0,  -48,  -42, 0,   0,   84,  31),
    )

    c0s = 372
    c3s = (
    ( 18,  -36, -122,   0,  120,  135, 0,   0,  -84, -31),
    (-18,   36,   -2,   0, -120,   51, 0,   0,   84, -31),
    ( 36, -165,  -27,  93,  147,   -9, 0, -93,   18,   0),
    (210,   45, -111, -93,  -57, -192, 0,  93,   12,  93),
    (162,  141,  -75, -93, -129, -180, 0,  93,  -12,  93),
    (-36,  -21,   27,  93,   39,    9, 0, -93,  -18,   0),
    (  0,    0,   62,   0,    0,   31, 0,   0,    0, -31),
    (  0,    0,  124,   0,    0,   62, 0,   0,    0, -62),
    (  0,    0,  124,   0,    0,   62, 0,   0,    0, -62),
    (  0,    0,   62,   0,    0,   31, 0,   0,    0, -31),
    (-18,   36,  -64,   0,   66,   51, 0,   0, -102,  31),
    ( 18,  -36,    2,   0,  -66,  -51, 0,   0,  102,  31),
    )

    def __init__(self, name="egm2008-1.pgm"):
        self.offset = None
        self.scale = None

        with open(name, "r") as f:
            line = f.readline()
            if line != "P5\012" and line != "P5\015\012":
                raise GeoidBadDataFile("No PGM header")
            headerlen = len(line)
            while True:
                line = f.readline()
                if len(line) == 0:
                    raise GeoidBadDataFile("EOF before end of file header")
                headerlen += len(line)
                if line.startswith('# Offset '):
                    try:
                        self.offset = int(line[9:])
                    except ValueError as e:
                        raise GeoidBadDataFile("Error reading offset", e)
                elif line.startswith('# Scale '):
                    try:
                        self.scale = float(line[8:])
                    except ValueError as e:
                        raise GeoidBadDataFile("Error reading scale", e)
                elif not line.startswith('#'):
                    try:
                        self.width, self.height = map(int, line.split())
                    except ValueError as e:
                        raise GeoidBadDataFile("Bad PGM width&height line", e)
                    break
            line = f.readline()
            headerlen += len(line)
            levels = int(line)
            if levels != 65535:
                raise GeoidBadDataFile("PGM file must have 65535 gray levels")
            if self.offset is None:
                raise GeoidBadDataFile("PGM file does not contain offset")
            if self.scale is None:
                raise GeoidBadDataFile("PGM file does not contain scale")

            if self.width < 2 or self.height < 2:
                raise GeoidBadDataFile("Raster size too small")

            fd = f.fileno()
            fullsize = os.fstat(fd).st_size

            if fullsize - headerlen != self.width * self.height * 2:
                raise GeoidBadDataFile("File has the wrong length")

            self.headerlen = headerlen
            self.raw = mmap.mmap(fd, fullsize, mmap.MAP_SHARED, mmap.PROT_READ)

        self.rlonres = self.width / 360.0
        self.rlatres = (self.height - 1) / 180.0
        self.ix = None
        self.iy = None

    def _rawval(self, ix, iy):
        if iy < 0:
            iy = -iy;
            ix += self.width/2;
        elif iy >= self.height:
            iy = 2 * (self.height - 1) - iy;
            ix += self.width/2;
        if ix < 0:
            ix += self.width;
        elif ix >= self.width:
            ix -= self.width

        return struct.unpack_from('>H', self.raw,
            (iy * self.width + ix) * 2 + self.headerlen)[0]

    def get(self, lat, lon, cubic=True):
        if lon < 0:
            lon += 360
        fy = (90 - lat) * self.rlatres
        fx = lon * self.rlonres
        iy = int(fy)
        ix = int(fx)
        fx -= ix
        fy -= iy
        if iy == self.height - 1:
            iy -= 1
        if ix != self.ix or iy != self.iy:
            self.ix = ix
            self.iy = iy
        if not cubic:
            self.v00 = self._rawval(ix, iy)
            self.v01 = self._rawval(ix+1, iy)
            self.v10 = self._rawval(ix, iy+1)
            self.v11 = self._rawval(ix+1, iy+1)
        else:
            v = (
            self._rawval(ix    , iy - 1),
            self._rawval(ix + 1, iy - 1),
            self._rawval(ix - 1, iy    ),
            self._rawval(ix    , iy    ),
            self._rawval(ix + 1, iy    ),
            self._rawval(ix + 2, iy    ),
            self._rawval(ix - 1, iy + 1),
            self._rawval(ix    , iy + 1),
            self._rawval(ix + 1, iy + 1),
            self._rawval(ix + 2, iy + 1),
            self._rawval(ix    , iy + 2),
            self._rawval(ix + 1, iy + 2)
            )
            
        if iy == 0:
            c3x = GeoidHeight.c3n
            c0x = GeoidHeight.c0n
        elif iy == self.height - 2:
            c3x = GeoidHeight.c3s
            c0x = GeoidHeight.c0s
        else:
            c3x = GeoidHeight.c3
            c0x = GeoidHeight.c0
            self.t = [
            sum([ v[j] * c3x[j][i] for j in range(12) ]) / float(c0x)
            for i in range(10)]
    
        if not cubic:
            a = (1 - fx) * self.v00 + fx * self.v01
            b = (1 - fx) * self.v10 + fx * self.v11
            h = (1 - fy) * a + fy * b
        else:
            h = (self.t[0] +
            fx * (self.t[1] + fx * (self.t[3] + fx * self.t[6])) +
            fy * (self.t[2] + fx * (self.t[4] + fx * self.t[7]) +
            fy * (self.t[5] + fx * self.t[8] + fy * self.t[9]))
        )
        return self.offset + self.scale * h
 
 
def utm2latlon(x,y):
   wgs84=pyproj.Proj("+init=EPSG:4326") #LatLon with WGS84 datum used by GPS units and Google Earth
   UTM29N=pyproj.Proj("+init=EPSG:32629") #UTM coords, zone 29N, WGS84 datum
   lon_wgs84,lat_wgs84 = pyproj.transform(UTM29N, wgs84, x, y)
   return lon_wgs84, lat_wgs84